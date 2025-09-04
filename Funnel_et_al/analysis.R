
rm(list = ls())
library(data.table)
source("utils.R")
library(tidyverse)

sample_names = get_sample_names()
sample_id = sample_names[10]
MODE = "fitch"

df_clonal_discordance = lapply(sample_names, function(sample_id) {
  print(sample_id)
  if (file.exists(paste0("results/bridges_trees/", sample_id, ".rds"))) {
    clonal_discordance(sample_id, mode = MODE)
  }
}) %>% do.call(rbind, .)
saveRDS(df_clonal_discordance, "results/metrics/clonal_discordance.rds")

sample_id = "OV2295"
dissimilarities_df = lapply(sample_names, function(sample_id) {
  print(sample_id)
  compute_sibling_similarities(sample_id)
}) %>% do.call(rbind, .)
saveRDS(dissimilarities_df, "results/metrics/dissimilarities.rds")
