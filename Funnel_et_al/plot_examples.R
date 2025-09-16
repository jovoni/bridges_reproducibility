
rm(list = ls())
source("getters.R")
library(ape)
library(ggtree)
library(dplyr)
library(tibble)

# ---- Load tree & labels ----
df_clonal_discordance = readRDS("results/metrics/clonal_discordance.rds")

sample_ids = names(table(df_clonal_discordance$sample_id)[table(df_clonal_discordance$sample_id) == 3])

lapply(sample_ids, function(x) {
  s = sample_id = x
  
  obj <- readRDS(file.path("results/bridges_trees", paste0(s, ".rds")))
  tree <- obj$tree
  
  labels_df <- get_labels(s) %>%
    distinct(cell_id, .keep_all = TRUE)
  
  # Keep only overlapping tips
  keep <- intersect(labels_df$cell_id, tree$tip.label)
  if (length(keep) == 0L) stop("No overlapping cells between labels and tree tips.")
  
  tree <- ape::keep.tip(tree, keep)
  labels_df <- labels_df %>% filter(cell_id %in% keep)
  clonal_disc = df_clonal_discordance %>% dplyr::filter(method == "bridges", sample_id == s) %>% 
    pull(value)
  
  # (Optional) drop edge lengths if they cause warnings
  tree$edge.length <- NULL
  
  # ---- Build heatmap data (rows must be tips; columns are annotations) ----
  # Here we only use clone_id. Add more columns if you have them.
  hm <- labels_df %>%
    select(cell_id, clone_id) %>%
    mutate(clone_id = as.factor(clone_id)) %>%
    column_to_rownames("cell_id")
  
  # Order rows exactly like tip labels (CRUCIAL for gheatmap)
  hm <- hm[tree$tip.label, , drop = FALSE]
  
  # ---- Circular tree + heatmap ----
  circ <- ggtree(tree, layout = "circular")
  p1 <- gheatmap(
    circ, hm,
    offset = 0.5, width = 0.1,
    colnames = TRUE, colnames_angle = 1, colnames_offset_y = 0.25, legend_title = "Clone", font.size = 0
  ) + ggtitle("BRIDGES", subtitle = paste0("Clonal discordance = ", clonal_disc))
  print(p1)
  
  
  # Sitka
  labels_df = get_labels(sample_id)
  sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
  labels_df = labels_df %>% dplyr::mutate(cell_id = paste0("cell_", cell_id))
  #keep <- labels_df_sitka$cell_id
  keep = intersect(labels_df$cell_id, sitka_tree$tip.label)
  tree <- keep.tip(sitka_tree, keep)
  clonal_disc = df_clonal_discordance %>% dplyr::filter(method == "sitka", sample_id == s) %>% 
    pull(value)
  
  # (Optional) drop edge lengths if they cause warnings
  tree$edge.length <- NULL
  
  # ---- Build heatmap data (rows must be tips; columns are annotations) ----
  # Here we only use clone_id. Add more columns if you have them.
  hm <- labels_df %>%
    select(cell_id, clone_id) %>%
    mutate(clone_id = as.factor(clone_id)) %>%
    column_to_rownames("cell_id")
  
  # Order rows exactly like tip labels (CRUCIAL for gheatmap)
  hm <- hm[tree$tip.label, , drop = FALSE]
  
  # ---- Circular tree + heatmap ----
  circ <- ggtree(tree, layout = "circular")
  p2 <- gheatmap(
    circ, hm,
    offset = 0.5, width = 0.1,
    colnames = TRUE, colnames_angle = 1, colnames_offset_y = 0.25, legend_title = "Clone", font.size = 0
  ) + ggtitle("Sitka", subtitle = paste0("Clonal discordance = ", clonal_disc))
  print(p2)
  
  
  # Lazac
  labels_df = get_labels(sample_id)
  newick_text = read_file(paste0("signatures_dataset/lazac_trees/",sample_id,"_hscn_tree.newick"))
  newick_text = paste0(newick_text, ";")
  tree = ape::read.tree(text = newick_text)
  keep = intersect(labels_df$cell_id, tree$tip.label)
  tree <- keep.tip(tree, keep)
  #Xdiss = get_X_for_dissimilarity(sample_id)
  labels_df <- labels_df %>%
    dplyr::filter(cell_id %in% keep)
  
  tree$edge.length <- NULL
  
  # ---- Build heatmap data (rows must be tips; columns are annotations) ----
  # Here we only use clone_id. Add more columns if you have them.
  hm <- labels_df %>%
    select(cell_id, clone_id) %>%
    mutate(clone_id = as.factor(clone_id)) %>%
    column_to_rownames("cell_id")
  
  # Order rows exactly like tip labels (CRUCIAL for gheatmap)
  hm <- hm[tree$tip.label, , drop = FALSE]
  
  clonal_disc = df_clonal_discordance %>% dplyr::filter(method == "lazac", sample_id == s) %>% 
    pull(value)
  
  # ---- Circular tree + heatmap ----
  circ <- ggtree(tree, layout = "circular")
  p3 <- gheatmap(
    circ, hm,
    offset = 0.5, width = 0.1,
    colnames = TRUE, colnames_angle = 1, colnames_offset_y = 0.25, legend_title = "Clone", font.size = 0
  ) + ggtitle("Lazac", subtitle = paste0("Clonal discordance = ", clonal_disc))
  print(p3)
  
  # Plot also Sibling dissimilarity
  df = readRDS("results/metrics/dissimilarities.rds")
  dissimilarity_density = df %>% 
    dplyr::filter(sample_id == s) %>% 
    ggplot(mapping = aes(x = diss, col = name)) +
    geom_density() +
    theme_bw() +
    scale_color_bmj() +
    labs(x = "Sibling Dissimilarity", y = "Density", col = "Algorithm")
  
  design = "
  AABBCC
  AABBCC
  AABBCC
  #DDDD#
  "
  
  p = p1 + p2 + p3 + 
    dissimilarity_density +
    plot_layout(guides = "collect", design = design) +
    plot_annotation(tag_levels = "A")
  ggsave(paste0("plot/samples/", x, ".pdf"), plot = p, width = 15, height = 8)
})

#sample_id = s = get_sample_names()[12]
