
rm(list = ls())
suppressPackageStartupMessages({
  library(bridges)
  library(data.table)
  library(ape)
})

source("getters.R")

# 1) Read index from CLI or fallback to SLURM_ARRAY_TASK_ID or 1
args <- commandArgs(trailingOnly = TRUE)
idx <- if (length(args) >= 1) as.integer(args[1]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (is.na(idx) || idx < 1) stop("Invalid 'idx' provided. It must be a positive integer.", call. = FALSE)

# 2) Get samples and bounds-check
samples <- get_sample_names()
n <- length(samples)
if (n == 0) stop("No samples returned by get_sample_names().", call. = FALSE)
if (idx > n) stop(sprintf("idx=%d exceeds number of samples=%d.", idx, n), call. = FALSE)

sample_id <- samples[idx]

# 3) Ensure output dir exists
out_dir <- file.path("results", "bridges_trees")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 4) Fit + timing (with basic error capture)
s <- Sys.time()
fit <- bridges::fit(get_input(sample_id), alleles = "CN", k_jitter_fix = 2)
e <- Sys.time()

res <- list(
  sample_id = sample_id,
  tree      = fit$tree,
  time      = e - s
)

saveRDS(res, file.path(out_dir, paste0(sample_id, ".rds")))
cat(sprintf("[OK] %s -> %s (elapsed: %s)\n",
            sample_id,
            file.path(out_dir, paste0(sample_id, ".rds")),
            e - s))
