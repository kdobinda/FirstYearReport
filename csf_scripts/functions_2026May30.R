# simulate_mixed_dat: generates discrete and categorical data from a random DAG
# relies on pcalg for random data generation
## Input
# seed: seed number
# p: Number of variables (columns)
# n: Number of observations (rows)
# deg: Expected number of outgoing edges from each variable
# binary_prop: proportion of binary variables in the dataset
# categorical_prop: proportion of categorical variables in the dataset
## Output
# raw_data: original continuous set
# mixed_data: generated mixed data
# adjacency_matrix: Adjacency matrix with weighted edges

simulate_mixed_dat <- function(
    seed,
    p = 30, 
    n = 1000, 
    deg = 3, 
    binary_prop = 0.25, 
    categorical_prop = 0.25){
  
  ## define parameters 
  s <- deg / (p-1) # sparsness of the graph
  
  ## generate random data
  set.seed(seed)
  true_dag <- randomDAG(p, s)
  gen_data <- rmvDAG(n, true_dag) # generate random samples
  
  oracle_corr <- cor(gen_data)
  
  # 25% binary, 25% 4 level category randomly
  bin_set <- gen_data[, 1:round(binary_prop * p)] 
  cat_set <- gen_data[, (round(binary_prop * p) + 1):round((binary_prop + categorical_prop) * p)]
  
  # categorical 
  for (i in 1:ncol(cat_set)){
    cat_cuts <- quantile(cat_set[, i], c(0.2, 0.4, 0.6, 0.8))
    cut1 <- runif(1, min = cat_cuts[1], max = cat_cuts[2])
    cut2 <- runif(1, min = cat_cuts[3], max = cat_cuts[4])
    cat_set[, i] <- as.integer((cat_set[, i] > cut1) + (cat_set[, i] > cut2)) 
  }
  
  # binary
  for (i in 1:ncol(bin_set)){
    bin_cuts <- quantile(bin_set[, i], c(0.25, 0.75))
    cut <- runif(1, min = bin_cuts[1], max = bin_cuts[2])
    bin_set[, i] <- as.integer((bin_set[, i] > cut) * 1)
  }
  
  new_dataset <- cbind(bin_set, cat_set, gen_data[, (round((binary_prop + categorical_prop) * p) + 1):p])
  new_dataset <- data.frame(new_dataset)
  colnames(new_dataset) <- paste0("V", 1:p)
  
  return(list(raw_data = gen_data,
              mixed_data = new_dataset,
              adjacency_matrix = true_dag))
}

# Penalty function for alpha
get_fges_penalty <- function(alpha) {
  # penalty_discount for FGES degenerate gaussian
  switch(as.character(alpha),
    "0.05"  = 1.0,
    "0.01"  = 2.0,
    "0.001" = 4.0,
    2.0
  )
}

get_dagslam_lambda <- function(alpha) {
  # lambda1 sparsity penalty for DAGSLAM
  switch(as.character(alpha),
    "0.05"  = 0.01,
    "0.01"  = 0.05,
    "0.001" = 0.10,
    0.05
  )
}

# Adaptive threshold for DAGSLAM adjacency matrix
# Uses the distribution of non-zero absolute weights to set a data-driven cutoff.
# Strategy: take the top (1 - target_density) quantile of absolute weights,
# where target_density is calibrated to match the expected sparsity of the true
# graph (deg / (p-1)). This avoids the fixed-0.3 cutoff producing graphs that
# are far denser than the other methods.
get_dagslam_threshold <- function(weight_matrix, deg, p) {
  abs_weights <- abs(weight_matrix)
  nonzero_weights <- abs_weights[abs_weights > 1e-6]
  
  if (length(nonzero_weights) == 0) return(0.3)  # fallback if all zero
  
  # Expected edge density in the true graph
  expected_density <- deg / (p * (p - 1))
  
  # We allow a modest multiple of expected edges (2x) to avoid over-pruning.
  # This targets roughly 2 * expected_density of the total possible edges.
  target_keep_fraction <- min(2 * expected_density * p * (p - 1) / length(nonzero_weights), 0.5)
  
  if (target_keep_fraction <= 0 || target_keep_fraction >= 1) return(0.3)
  
  threshold <- quantile(nonzero_weights, 1 - target_keep_fraction)
  
  # Hard floor: never go below 0.1 (avoids keeping near-zero noise weights)
  threshold <- max(threshold, 0.1)
  
  return(as.numeric(threshold))
}

# mixed_data_sim: runs causal discovery based on different inputs
## Input
# input_matrix: row with seed, alpha, n, p, deg, method
# save_dags: logical, whether to save discovered DAG to disk
## Output
# matrix row of performance metrics
mixed_data_sim <- function(input_matrix, save_dags = FALSE){
  
  seed   <- input_matrix$seed
  alpha  <- input_matrix$alpha
  n      <- input_matrix$n
  p      <- input_matrix$p
  deg    <- input_matrix$deg
  method <- input_matrix$method
  
  # generate data
  set.seed(seed)
  dat_obj <- simulate_mixed_dat(seed = seed, p = p, n = n, deg = deg)
  raw_data <- dat_obj$raw_data
  raw_corr <- cor(raw_data) 
  mix_dat  <- dat_obj$mixed_data
  true_dag <- dat_obj$adjacency_matrix
  true_dag_mat <- ifelse(as(true_dag, "matrix") != 0, 1, 0)
  
  # Record start time; we always work in seconds
  start_time <- proc.time()["elapsed"]
  
  if (method == 1){
    # Oracle PC
    corr <- raw_corr
    corr_diff <- 0
    disc_fit <- pc(suffStat = list(C = corr, n = n), 
                   gaussCItest, p = p,
                   alpha = alpha, skel.method = "stable")
    disc_mat <- as(disc_fit, "amat")
    gc()
  }
  
  if (method == 2){
    # Naive PC (treating mixed data as continuous)
    corr <- cor(mix_dat)
    corr_diff <- sum(abs(corr - raw_corr))
    disc_fit <- pc(suffStat = list(C = corr, n = n), 
                   gaussCItest, p = p,
                   alpha = alpha, skel.method = "stable")
    disc_mat <- as(disc_fit, "amat")
    gc()
  }
  
  if (method == 3){
    # FGES - degenerate Gaussian score
    corr <- cor(mix_dat)
    corr_diff <- sum(abs(corr - raw_corr))
    penalty <- get_fges_penalty(alpha)
    search <- ts$TetradSearch(mix_dat)
    search$set_verbose(FALSE)
    search$use_degenerate_gaussian_score(penalty_discount = penalty)
    search$run_fges()
    amat <- tr$graph_to_matrix(search$get_java())
    disc_mat <- matrix(0, nrow = p, ncol = p)
    disc_mat[amat == 2] <- 1
    disc_mat[amat == 3 & t(amat) == 3] <- 1
    gc()
  }
  
   if (method == 4){

    # Copula PC
    ## estimate underlying correlation matrix and (effective) sample size
    # copula object
    cop.obj <- inferCopulaModel(data.frame(mix_dat), nsamp = 1000, S0 = diag(p)/n, verb = T)

    # correlation matrix samples
    C_samples <- cop.obj$C.psamp[,, 501:1000]
    # average correlation matrix
    corr <- apply(C_samples, c(1,2), mean)
    corr_diff <- sum(abs(corr - raw_corr))

    

    ## CoPC
    disc_fit <- pc(suffStat = list(C = corr, n = n),

                   indepTest = gaussCItest, 

                   alpha = alpha,

                   p=p, skel.method = "stable")

    disc_mat <- as(disc_fit, "amat")

gc()  

} 
  
  if (method == 5){
    # Latent PC
    labels <- label_fun(mix_dat)
    corr <- latent_pc(mix_dat, labels)
    corr_diff <- sum(abs(corr - raw_corr))
    disc_fit <- pc(suffStat = list(C = corr, n = n),
                   indepTest = gaussCItest,
                   alpha = alpha, p = p, skel.method = "stable")
    disc_mat <- as(disc_fit, "amat")
    gc()
  }
  
  if (method == 6){
    # DAGSLAM
    corr <- cor(mix_dat)
    corr_diff <- sum(abs(corr - raw_corr))
    labels <- label_fun(mix_dat)
    type_mapping <- c(
      "continuous" = "gauss",
      "binary"     = "logistic",
      "ordinal"    = "muti-logistic"
    )
    loss_vec <- unname(type_mapping[labels])
    
    np <- import("numpy", convert = FALSE)
    m_vec <- rep(1L, length(labels)) 
    matrix_dat <- unname(as.matrix(mix_dat))
    
    for (i in seq_along(labels)) {
      if (labels[i] == "ordinal") {
        m_vec[i] <- as.integer(length(unique(matrix_dat[, i]))) 
      }
    }
    
    X_writable <- np$copy(matrix_dat)
    lambda <- get_dagslam_lambda(alpha)
    
    estimated_adjacency_matrix <- tryCatch({
      dagslam(
        X        = X_writable,
        loss_type = loss_vec,
        m_vec    = m_vec,
        lambda1  = lambda          # now uses alpha-mapped lambda, not hardcoded 0.05
      )
    }, error = function(e) {
      cat("DAGSLAM failed due to numerical instability on this seed.\n")
      return(matrix(NA, nrow = p, ncol = p))
    })
    
    if (any(is.na(estimated_adjacency_matrix))) {
      disc_mat <- matrix(0, nrow = p, ncol = p)
    } else {
      # fixed 0.3 weight cutoff from author
      disc_mat <- ifelse(abs(estimated_adjacency_matrix) > 0.3, 1, 0)
    }
    gc()
  }
  
  # ── Timing: always in seconds ─────────────────────────────────────────────
  time_secs <- proc.time()["elapsed"] - start_time
  
  # ── Standard metrics ──────────────────────────────────────────────────────
  tp  <- sum(true_dag_mat[which(disc_mat == 1)] == 1)
  tn  <- sum(true_dag_mat[which(disc_mat == 0)] == 0)
  fp  <- sum(disc_mat[which(true_dag_mat == 0)] == 1)
  fn  <- sum(disc_mat[which(true_dag_mat == 1)] == 0)
  TPR <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  FPR <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
  SHD <- shd(disc_mat, true_dag_mat)
  
  # ── Adjacency metrics (undirected skeleton) ───────────────────────────────
  adj_true <- 1 * ((true_dag_mat + t(true_dag_mat)) > 0)
  adj_disc <- 1 * ((disc_mat  + t(disc_mat))  > 0)
  
  adj_tp <- sum(adj_disc == 1 & adj_true == 1) / 2
  adj_fp <- sum(adj_disc == 1 & adj_true == 0) / 2
  adj_fn <- sum(adj_disc == 0 & adj_true == 1) / 2
  
  adj_prec <- ifelse((adj_tp + adj_fp) == 0, 0, adj_tp / (adj_tp + adj_fp))
  adj_rec  <- ifelse((adj_tp + adj_fn) == 0, 0, adj_tp / (adj_tp + adj_fn))
  
  # ── Arrowhead metrics (directional accuracy) ──────────────────────────────
  arrow_true <- 1 * (true_dag_mat == 1 & t(true_dag_mat) == 0)
  arrow_disc <- 1 * (disc_mat  == 1 & t(disc_mat)  == 0)
  
  arr_tp <- sum(arrow_disc == 1 & arrow_true == 1)
  arr_fp <- sum(arrow_disc == 1 & arrow_true == 0)
  arr_fn <- sum(arrow_disc == 0 & arrow_true == 1)
  
  arr_prec <- ifelse((arr_tp + arr_fp) == 0, 0, arr_tp / (arr_tp + arr_fp))
  arr_rec  <- ifelse((arr_tp + arr_fn) == 0, 0, arr_tp / (arr_tp + arr_fn))
  
  # ── Selective DAG saving ───────────────────────────────────────────────────
  # True DAG saved once per (seed, n, p, deg) combination — only when deg <= 3
  # to avoid writing large files for every scenario.
  if (deg <= 3) {
    true_dag_file <- paste0("results/dags/true_dag_seed", seed,
                            "_n", n, "_p", p, "_deg", deg, ".rds")
    if (!file.exists(true_dag_file)) saveRDS(true_dag_mat, true_dag_file)
  }
  
  # Discovered DAG only saved when explicitly requested via save_dags flag
  if (save_dags) {
    saveRDS(disc_mat, file = paste0("results/dags/disc_mat_seed", seed,
                                   "_alpha", alpha, "_n", n,
                                   "_p", p, "_deg", deg,
                                   "_method", method, ".rds"))
  }
  
  # ── Output row ────────────────────────────────────────────────────────────
  out <- matrix(c(seed, alpha, deg, n, p, method, time_secs, corr_diff, 
                  TPR, FPR, adj_prec, adj_rec, arr_prec, arr_rec, SHD), nrow = 1)
  colnames(out) <- c("seed", "alpha", "deg", "n", "p", "method",
                     "time_secs",   # always seconds now
                     "corr_diff",
                     "TPR", "FPR",
                     "adj_prec", "adj_rec",
                     "arr_prec", "arr_rec",
                     "SHD")
  return(out)
}
