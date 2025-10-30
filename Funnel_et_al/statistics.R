
rm(list = ls())
library(kableExtra)
library(knitr)
library(tidyverse)
source("getters.R")

df = lapply(get_sample_names(), function(s) {
  x = readRDS(file.path("results/bridges_trees/", paste0(s, ".rds")))
  dplyr::tibble(sample = s, ncells = x$tree$tip.label %>% length())
}) %>% do.call("bind_rows", .)


# Clonal discordance
df_clonal_discordance = readRDS("results/metrics/clonal_discordance.rds")
df_clonal_discordance %>%
  tidyr::pivot_wider(names_from = method, values_from = value) %>%
  kbl(format = "latex",
      booktabs = TRUE,
      col.names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
      escape = FALSE)

df_clonal_discordance %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(method != "lazac") %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(value == min(value)) %>%
  dplyr::group_by(Dataset, method) %>%
  dplyr::summarise(n = n())

df_clonal_discordance %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::filter(Dataset == "Funnel", method != "sitka") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(value == min(value)) %>%
  dplyr::group_by(Dataset, method) %>%
  dplyr::summarise(n = n())

# Dissimilarities
df = readRDS("results/metrics/dissimilarities.rds")

df %>%
  dplyr::mutate(name = factor(name, levels = c("bridges", "sitka", "lazac"))) %>%
  na.omit() %>%
  dplyr::group_by(name, sample_id) %>%
  dplyr::summarise(m = mean(diss)) %>%
  dplyr::mutate(m = formatC(m, format = "e", digits = 2)) %>%
  tidyr::pivot_wider(names_from = name, values_from = m) %>%
  kbl(format = "latex",
      booktabs = TRUE,
      col.names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
      escape = FALSE)

df %>%
  na.omit() %>%
  dplyr::group_by(name, sample_id) %>%
  dplyr::summarise(m = mean(diss)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(name != "lazac") %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(m == min(m)) %>%
  dplyr::group_by(Dataset, name) %>%
  dplyr::summarise(n = n())

df %>%
  na.omit() %>%
  dplyr::group_by(name, sample_id) %>%
  dplyr::summarise(m = mean(diss)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::filter(Dataset == "Funnel", name != "sitka") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(m == min(m)) %>%
  dplyr::group_by(Dataset, name) %>%
  dplyr::summarise(n = n())


# Youden
df = readRDS("results/metrics/youden.rds")

df %>%
  dplyr::mutate(mehtod = factor(mehtod, levels = c("bridges", "sitka", "lazac"))) %>%
  dplyr::mutate(Jouden = formatC(Jouden, format = "e", digits = 2)) %>%
  tidyr::pivot_wider(names_from = mehtod, values_from = Jouden) %>%
  kbl(format = "latex",
      booktabs = TRUE,
      col.names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
      escape = FALSE)

colnames(df)
df %>%
  na.omit() %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::filter(mehtod != "lazac") %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(Jouden == max(Jouden)) %>%
  dplyr::group_by(Dataset, mehtod) %>%
  dplyr::summarise(n = n())

df %>%
  na.omit() %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::filter(Dataset == "Funnel", mehtod != "sitka") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::filter(Jouden == max(Jouden)) %>%
  dplyr::group_by(Dataset, mehtod) %>%
  dplyr::summarise(n = n())
