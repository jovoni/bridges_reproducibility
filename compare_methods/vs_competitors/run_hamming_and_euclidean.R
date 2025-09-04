
# Load required library
library(ape)
library(dplyr)
library(readr)

# -------- Function to compute pairwise distances -------- #
compute_pairwise_distance <- function(input_tsv, method = c("hamming", "euclidean")) {
  method <- match.arg(method)

  # Read input and order bins
  df <- readr::read_tsv(input_tsv, show_col_types = FALSE) %>%
    arrange(chr, pos)

  # Extract copy number matrix (rows: bins, columns: cells)
  cn_matrix <- as.matrix(df[ , !(names(df) %in% c("chr", "pos")) ])
  cn_matrix <- apply(cn_matrix, 2, as.numeric)

  cell_ids <- colnames(cn_matrix)
  n <- ncol(cn_matrix)
  D <- matrix(0, nrow = n, ncol = n, dimnames = list(cell_ids, cell_ids))

  # Compute distances
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (method == "hamming") {
        d <- sum(cn_matrix[, i] != cn_matrix[, j], na.rm = TRUE)
      } else {
        d <- sqrt(sum((cn_matrix[, i] - cn_matrix[, j])^2, na.rm = TRUE))
      }
      D[i, j] <- D[j, i] <- d
    }
  }

  return(D)
}

# -------- Main Script -------- #

# Command-line arguments: input_tsv outdir
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_tsv> <outdir>")
}

tsv_path <- args[1]
outdir <- args[2]

# Create output directory if missing
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Loop through methods
for (method in c("hamming", "euclidean")) {
  message("Processing method: ", method)
  s <- Sys.time()

  # Compute distance matrix and NJ tree
  D <- compute_pairwise_distance(tsv_path, method = method)
  tree <- ape::nj(D)

  e <- Sys.time()
  time_elapsed <- as.numeric(difftime(e, s, units = "secs"))
  message("Time for ", method, ": ", round(time_elapsed, 2), " seconds")

  # Create subdirectory for this method
  method_dir <- file.path(outdir, method)
  dir.create(method_dir, showWarnings = FALSE, recursive = TRUE)

  # Save outputs
  saveRDS(D, file = file.path(method_dir, "distance_matrix.rds"))
  saveRDS(tree, file = file.path(method_dir, "nj_tree.rds"))
  writeLines(as.character(time_elapsed), con = file.path(method_dir, "runtime_seconds.txt"))
}

message("Done.")
