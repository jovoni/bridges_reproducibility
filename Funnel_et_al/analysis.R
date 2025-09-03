
rm(list = ls())
library(data.table)
source("utils.R")
library(tidyverse)

sample_names = get_sample_names()

sample_id = sample_names[31]
sample_id = "SA1035"

MODE = "fitch"

df = lapply(sample_names, function(sample_id) {
  print(sample_id)
  if (file.exists(paste0("results/bridges_trees/", sample_id, ".rds"))) {
    
    clonal_discordance(sample_id, mode = MODE)
    
    # fit = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))
    # 
    # dplyr::tibble(
    #   sample_id = sample_id,
    #   bridges = bridges_score(fit, sample_id, mode = MODE),
    #   sitka = sitka_score(sample_id, mode = MODE)
    # )
  }
}) %>% do.call(rbind, .)

df %>% 
  ggplot(mapping = aes(x = sample_id, y = value, fill = method)) +
  geom_col(position = "dodge") +
  theme_bw()

df %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::filter(value == min(value)) %>% 
  pull(method) %>% 
  table()



