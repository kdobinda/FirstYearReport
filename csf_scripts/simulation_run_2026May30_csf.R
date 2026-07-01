library(reticulate)
use_virtualenv("/mnt/iusers01/bk01-icvs/z55517kd/mixed_sim/.venv", required = TRUE)
py_config()
source_python('/mnt/iusers01/bk01-icvs/z55517kd/mixed_sim/python/import_jpype.py')
source_python('/mnt/iusers01/bk01-icvs/z55517kd/mixed_sim/python/DAGSLAM.py')
library(pcalg)
library(mvtnorm)
library(dplyr)
library(parallel)

# Source helper scripts explicitly from your home directory path
home_path <- file.path(Sys.getenv("HOME"), "mixed_sim")

source(file.path(home_path, 'gaussCItestLocal.R'))
source(file.path(home_path, 'inferCopulaModel.R'))  
source(file.path(home_path, "latent_pc.R"))
source(file.path(home_path, "utility.R"))
source(file.path(home_path, "functions_2026May30.R"))   # <-- updated functions

ts <- import("pytetrad.tools.TetradSearch")
tr <- import("pytetrad.tools.translate")

dir.create("results/dags", recursive = TRUE, showWarnings = FALSE)

# ── Job identity ──────────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
task_id   <- if (length(args) > 0) as.numeric(args[1]) else 1
my_p      <- if (length(args) > 1) as.numeric(args[2]) else 100
my_mgroup <- if (length(args) > 2) args[3] else "fast"
my_run_id <- if (length(args) > 3) args[4] else "default_run"

my_seed   <- 50 + task_id

cat("====================================================\n")
cat("Task ID:", task_id, "\n")
cat("Seed:", my_seed, "| p:", my_p, "| Method group:", my_mgroup, "\n")
cat("Run ID:", my_run_id, "\n")
cat("====================================================\n")

# ── Simulation parameters ─────────────────────────────────────────────────────
alpha <- c(0.001, 0.01, 0.05)
n     <- 5000
deg   <- c(1, 3, 5)   # updated: three sparsity levels as requested

# Method mapping
#   1 = Oracle PC
#   2 = Naive PC
#   3 = FGES (degenerate Gaussian)
#   4 = Copula PC  -- intentionally excluded from this run
#   5 = Latent PC
#   6 = DAGSLAM
#
# "fast"  group: methods 1, 2, 3, 5  (no Python/Java overhead for 1 & 2)
# "slow"  group: method 4 and 6 (DAGSLAM -- computationally expensive)
# "all"   group: 1, 2, 3, 4, 5, 6 (useful for small p validation runs)

method <- switch(my_mgroup,
  "fast" = c(1, 2, 3, 5),
  "slow" = c(4, 6),
  "all"  = c(1, 2, 3, 4, 5, 6),
  c(1, 2, 3, 5)   # fallback
)

input_matrix <- expand.grid(my_seed, alpha, n, my_p, deg, method)
colnames(input_matrix) <- c("seed", "alpha", "n", "p", "deg", "method")
# no copula for > 50 vars
input_matrix <- subset(input_matrix, !(p > 50 & method == 4))
colnames(input_matrix) <- c("seed", "alpha", "n", "p", "deg", "method")

cat("Total parameter combinations:", nrow(input_matrix), "\n")

# ── Parallelisation ───────────────────────────────────────────────────────────
# Use the CPUs allocated by SLURM. We cap at nrow(input_matrix) so we never
# request more workers than there are jobs.
#
# IMPORTANT: reticulate/JPype cannot be used inside forked workers because the
# JVM and Python interpreter are not fork-safe.  Methods 3 (FGES) and 6
# (DAGSLAM) depend on reticulate/JPype, so we split the work:
#   * Methods 1, 2, 5  --> safe to parallelise with mclapply (fork)
#   * Methods 3, 6     --> run sequentially to avoid JVM crashes
#
# If your mgroup is "slow" (method 6 only), everything runs sequentially.

n_cores_available <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat("CPUs available:", n_cores_available, "\n")

parallel_safe_methods <- c(1, 2, 4, 5)
needs_jvm_methods     <- c(3, 6)

rows_parallel   <- which(input_matrix$method %in% parallel_safe_methods)
rows_sequential <- which(input_matrix$method %in% needs_jvm_methods)

cat("Rows to parallelise:", length(rows_parallel),
    "| Rows sequential:", length(rows_sequential), "\n\n")

# save_dags = TRUE only for deg == 1 or deg == 3 runs, to keep disk usage sane
run_row <- function(i) {
  row <- input_matrix[i, ]
  save_flag <- (row$deg <= 3)
  tryCatch(
    mixed_data_sim(input_matrix = row, save_dags = save_flag),
    error = function(e) {
      cat("ERROR on row", i, ":", conditionMessage(e), "\n")
      out <- matrix(NA, nrow = 1, ncol = 15)
      colnames(out) <- c("seed", "alpha", "deg", "n", "p", "method",
                         "time_secs", "corr_diff",
                         "TPR", "FPR",
                         "adj_prec", "adj_rec",
                         "arr_prec", "arr_rec", "SHD")
      out[1, c("seed","alpha","deg","n","p","method")] <-
        c(row$seed, row$alpha, row$deg, row$n, row$p, row$method)
      return(out)
    }
  )
}

# Run parallel-safe methods in parallel
results_parallel <- list()
if (length(rows_parallel) > 0) {
  n_workers <- min(n_cores_available, length(rows_parallel))
  cat("Launching", n_workers, "parallel workers for methods", 
      paste(parallel_safe_methods, collapse = ","), "...\n")
  results_parallel <- mclapply(
    rows_parallel,
    run_row,
    mc.cores = n_workers,
    mc.preschedule = TRUE   # better load balance for variable-cost jobs
  )
}

# Run JVM-dependent methods sequentially
results_sequential <- list()
if (length(rows_sequential) > 0) {
  cat("Running", length(rows_sequential),
      "sequential rows for methods", paste(needs_jvm_methods, collapse = ","), "...\n")
  results_sequential <- lapply(seq_along(rows_sequential), function(j) {
    i <- rows_sequential[j]
    if (j %% 5 == 0) cat("  Sequential row", j, "of", length(rows_sequential), "\n")
    run_row(i)
  })
}

# Reassemble in original row order
all_results <- vector("list", nrow(input_matrix))
for (j in seq_along(rows_parallel))   all_results[[rows_parallel[j]]]   <- results_parallel[[j]]
for (j in seq_along(rows_sequential)) all_results[[rows_sequential[j]]] <- results_sequential[[j]]

out <- as.data.frame(do.call(rbind, all_results))

# ── Save results ──────────────────────────────────────────────────────────────
out_filename <- sprintf("results/sim_results_seed%d_p%d_%s_%s.csv",
                        my_seed, my_p, my_mgroup, my_run_id)
write.csv(out, file = out_filename, row.names = FALSE)
cat("Done. Results written to:", out_filename, "\n")
