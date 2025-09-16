
rm(list = ls())
library(data.table)
source("utils.R")
library(tidyverse)

sample_names = get_sample_names()
sample_id = sample_names[2]
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

# RF distances
df_v_sitka = lapply(sample_names, function(sample_id) {
  print(sample_id)
  compute_rf_distance_w_sitka(sample_id)  
}) %>% do.call("bind_rows", .)

df_v_lazac = lapply(sample_names, function(sample_id) {
  print(sample_id)
  compute_rf_distance_w_lazac(sample_id)
}) %>% do.call("bind_rows", .)

saveRDS(dplyr::bind_rows(df_v_sitka, df_v_lazac), "results/metrics/RF_against.rds")
