
suppressPackageStartupMessages({
  library(bridges)
  library(data.table)
  library(ape)
})

rm(list = ls())
source("utils.R")

samples <- get_sample_names()
sample_id <- "DG1134"

data = get_input(sample_id)

length(unique(data$cell_id))

s <- Sys.time()
fit <- bridges::fit(data, k_jitter_fix = 0, alleles = "CN")
e <- Sys.time()
e-s

fit0 = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))

# get_labels(sample_id) %>% 
#   dplyr::rename(name=cell_id, branch_color = clone_id) %>% 
#   write_delim("results/mapping.txt", delim = "\t")
# 
# library(RColorBrewer)
# mapping <- get_labels(sample_id) %>%
#   # pick a palette (expand if needed)
#   mutate(
#     branch_color = as.factor(clone_id),
#     branch_color = setNames(
#       brewer.pal(max(3, length(unique(clone_id))), "Set3"),
#       levels(as.factor(clone_id))
#     )[branch_color]
#   ) %>%
#   rename(name = cell_id) %>% 
#   dplyr::select(name, branch_color)
#write_delim(mapping, "results/mapping.txt", delim = "\t")


bridges_score(fit, sample_id, mode = "fitch")
bridges_score(fit0, sample_id, mode = "fitch")
sitka_score(sample_id, mode = "fitch")

sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))

bridges::plot_heatmap(data, tree = fit$tree, annotations = get_labels(sample_id))



labels_df = get_labels(sample_id)
sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
labels_df_sitka = labels_df %>% dplyr::mutate(cell_id = paste0("cell_", cell_id))
keep <- labels_df_sitka$cell_id
tree_pruned <- keep.tip(sitka_tree, labels_df_sitka$cell_id)
labels_pruned <- labels_df_sitka %>%
  dplyr::filter(cell_id %in% keep)
bridges::plot_heatmap(data %>% dplyr::mutate(cell_id = paste0("cell_", cell_id)), tree = tree_pruned, annotations = labels_pruned)
