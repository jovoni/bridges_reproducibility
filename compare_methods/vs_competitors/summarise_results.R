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
i = 2

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
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "dice", seconds = root_time))

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
  
  # Sitka
  sitka_tree = ape::read.tree(file = file.path(RES_DIR, sim_id, "sitka/tree.newick"))
  sitka_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id, "sitka/time.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "sitka", seconds = sitka_time))
  
  # Lazac
  lazac_tree = ape::read.tree(file = file.path(RES_DIR, sim_id, "lazac/_tree.newick"))
  lazac_time = as.numeric(unlist(strsplit(read_file(file.path(RES_DIR, sim_id, "lazac/time.txt")), "\n")))
  time_df = dplyr::bind_rows(time_df, dplyr::tibble(method = "lazac", seconds = lazac_time))

  trees = list(
    "bridges" = bridges_tree,
    "dice" = root_tree,
    "medicc" = medicc_tree,
    "hamming" = hamming_tree,
    "euclidean" = euc_tree,
    "sitka" = sitka_tree,
    "lazac" = lazac_tree
  )

  time_df = time_df %>%
    dplyr::mutate(algorithm = method) %>%
    dplyr::select(algorithm, seconds) %>%
    tidyr::pivot_longer(!algorithm, names_to = "metric")
  
  get_tree_metrics = function(tree1, tree2) {
    common_taxa <- intersect(tree1$tip.label, tree2$tip.label)
    tree1 <- ape::drop.tip(tree1, setdiff(tree1$tip.label, common_taxa))
    tree2 <- ape::drop.tip(tree2, setdiff(tree2$tip.label, common_taxa))
    
    rf_dist <- phangorn::RF.dist(tree1, tree2)
    rf_normalized <- phangorn::RF.dist(tree1, tree2, normalize = T)
    
    GRF = TreeDist::TreeDistance(tree1, tree2)
    MutualClusteringInf = TreeDist::MutualClusteringInfo(tree1, tree2, normalize = T)
    
    if (tree1$Nnode >= 477) {
      quartet_divergence = NA
    } else {
      
      calculate_q = function(tree1, tree2) {
        tryCatch({
          sq_status <- Quartet::QuartetStatus(tree1, tree2)
          quartet_divergence = Quartet::QuartetDivergence(sq_status, similarity = F)
          return(quartet_divergence)
        }, error = function(e) {
          dist_mat <- ape::cophenetic.phylo(tree2)
          tree2_rebuilt <- ape::nj(dist_mat)  # neighbor-joining reconstruction
          tree2_rebuilt$tip.label <- rownames(dist_mat)
          sq_status <- Quartet::QuartetStatus(tree1, tree2_rebuilt)
          quartet_divergence = Quartet::QuartetDivergence(sq_status, similarity = F)
          return(quartet_divergence)
        })
      }
      
      quartet_divergence = calculate_q(tree1, tree2)
    }
    
    names = c("RF distance", "Quartet divergence", "RF normalized", "Generalized RF", "Mutual Clustering Info")
    values = c(rf_dist, quartet_divergence, rf_normalized, GRF, MutualClusteringInf)
    dplyr::tibble(metric = names, value = values)
  }
  
  sim_df = lapply(1:length(trees), function(j) {
    get_tree_metrics(true_tree, trees[[j]]) %>%
      dplyr::mutate(algorithm = names(trees)[j])
  }) %>% do.call("bind_rows", .) %>% dplyr::bind_rows(time_df)
  
  #sim_df %>% tidyr::pivot_wider(values_from = value, names_from = algorithm) %>% view()

  cbind(sim_df, sim)
}, mc.cores = 1) %>% do.call("bind_rows", .)

saveRDS(RES, file.path(RES_DIR, "results_summary.RDS"))
