# Define exact scratch workspace paths
scratch_results_dir <- file.path(Sys.getenv("HOME"), "scratch", "mixed_sim_runs", "run_May30_redo", "results")
output_master_file  <- file.path(Sys.getenv("HOME"), "mixed_sim", "results", "simulation_results_060226.csv")
# output_summary_file <- file.path(Sys.getenv("HOME"), "mixed_sim", "summary_metrics_table.csv")

csv_files <- list.files(path = scratch_results_dir, pattern = "\\.csv$", full.names = TRUE)

# Systematically import and merge all entries while preserving complete integrity
all_data_list <- lapply(seq_along(csv_files), function(i) {
  read.csv(csv_files[i], stringsAsFactors = FALSE)
})
master_df <- do.call(rbind, all_data_list)

# Drop any accidental rows filled entirely with NA from failed runs
master_df <- master_df[!is.na(master_df$seed), ]

# Standardize method identities into distinct labels for clear reading
method_mapping <- c(
  "1" = "Oracle PC",
  "2" = "Naive PC",
  "3" = "fGES (DG Score)",
  "4" = "Copula PC",
  "5" = "Latent PC",
  "6" = "DAGSLAM"
)
master_df$method_name <- method_mapping[as.character(master_df$method)]

# Save the absolute un-edited, combined raw dataset
write.csv(master_df, file = output_master_file, row.names = FALSE)
