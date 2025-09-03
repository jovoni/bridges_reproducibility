
rm(list = ls())
library(tidyverse)

cn_data = read.delim("data/ov2295_cell_cn.csv", sep = ",")
clone_df = read.delim("data/ov2295_clone_clusters.csv", sep = ",")

N = 900
N = min(N, nrow(clone_df))
cell_ids = sample(clone_df$cell_id, size = N)

cn_data = cn_data %>%
  dplyr::filter(cell_id %in% cell_ids) %>%
  dplyr::mutate(CN = state) %>%
  dplyr::select(cell_id, sample_id, library_id, start, end, chr, CN)

bridges_fit = bridges::fit(
  data = cn_data,
  alleles = c("CN"),
  k_jitter_fix = 2
)

dir.create("results", recursive = T)
saveRDS(bridges_fit, "results/bridges_fit.rds")

bridges_fit_wo_reconstruction = bridges::fit(
  data = cn_data,
  alleles = c("CN"),
  k_jitter_fix = 2, 
  avoid_reconstruction = T
)

dir.create("results", recursive = T)
saveRDS(bridges_fit_wo_reconstruction, "results/bridges_fit_wo_reconstruction.rds")
