rm(list = ls())
require(tidyverse)
library(ape)
library(ggplot2)
library(tidyverse)
source("../utils_comparisons.R")

# source("../../functions_for_tree_algorithm/compute_distance_matrix.R")
# source("../../functions_for_tree_algorithm/distance_functions.R")
# source("../../functions_for_tree_algorithm/main.R")
# source("../../functions_for_tree_algorithm/preprocess_input.R")
# source("../../functions_for_tree_algorithm/plot.R")

params_df = read.delim("../data/param_grid.csv", sep = ",")
DATA_DIR = "../data/"
RES_DIR = "results/"
dir.create(RES_DIR, recursive = T)
i = 1

message("Reading params")
params_df = read.delim(file.path(DATA_DIR, "param_grid.csv"), sep = ",")
message("Starting inference")

RES = parallel::mclapply(1:nrow(params_df), function(i) {
#RES = parallel::mclapply(1:N, function(i) {
  print(i)
  sim = params_df[i,]
  sim_id = sim$sim_id

  data = readRDS(file.path(DATA_DIR, sim_id, "simulation.RDS"))
  true_tree = data$tree

  # Bridges
  bridges_tree = readRDS(file.path(RES_DIR, sim_id, "bridges_tree.RDS"))
  time_df = readRDS(file.path(RES_DIR, sim_id, "time_tibble.RDS"))

  # DICE
  root_tree = ape::read.tree(file.path(RES_DIR, sim_id,"dice", "standard_root_balME_tree.nwk"))
  root_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id,"dice", "dice_time.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "root", seconds = root_time))

  # MEDICC
  medicc_tree = ape::read.tree(file.path(RES_DIR, sim_id, "medicc2", "medicc_input_final_tree.new"))
  medicc_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id, "medicc2", "medicc_time.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "medicc", seconds = medicc_time))

  # Hamming
  hamming_tree = readRDS(file.path(RES_DIR, sim_id, "hamming/nj_tree.rds"))
  hamming_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id, "hamming/runtime_seconds.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "hamming", seconds = hamming_time))

  # Euclidean
  euc_tree = readRDS(file.path(RES_DIR, sim_id, "euclidean/nj_tree.rds"))
  euc_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id, "euclidean/runtime_seconds.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "euclidean", seconds = hamming_time))

  trees = list(
    "bridges" = bridges_tree,
    "dice" = root_tree,
    "medicc" = medicc_tree,
    "hamming" = hamming_tree,
    "euclidean" = euc_tree
  )

  time_df = time_df %>%
    dplyr::mutate(algorithm = method) %>%
    dplyr::select(algorithm, seconds) %>%
    tidyr::pivot_longer(!algorithm, names_to = "metric")

  sim_df = lapply(1:length(trees), function(j) {
    get_tree_metrics(true_tree, trees[[j]]) %>%
      dplyr::mutate(algorithm = names(trees)[j])
  }) %>% do.call("bind_rows", .) %>% dplyr::bind_rows(time_df)

  cbind(sim_df, sim)
}, mc.cores = 1) %>% do.call("bind_rows", .)

saveRDS(RES, file.path(RES_DIR, "results_summary.RDS"))
