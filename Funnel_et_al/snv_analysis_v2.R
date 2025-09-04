
rm(list = ls())
library(data.table)
library(tidyverse)
library(ape)
library(phangorn)

snv_annotations = fread("signatures_dataset/DLP/SNV/snv_annotations.csv.gz")
snv_all = fread("signatures_dataset/DLP/SNV/snv_counts.csv.gz")

s = "SA604"



# Your inputs:
# tree: phylo object
# snv_df: dataframe with columns cell_id, chr, start, alt_count

# First, create SNV presence/absence matrix
# snv_matrix <- snv_df %>%
#   dplyr::filter(alt_count > 0) %>%  # Only consider SNVs with reads supporting alternative
#   dplyr::mutate(snv_id = paste(chr, start, sep = "_")) %>%
#   select(cell_id, snv_id) %>%
#   distinct() %>%
#   mutate(present = 1) %>%
#   pivot_wider(names_from = snv_id, values_from = present, values_fill = 0)
# 
# # Convert to named list for easier lookup
# snv_support <- snv_matrix %>%
#   pivot_longer(-cell_id, names_to = "snv_id", values_to = "present") %>%
#   filter(present == 1) %>%
#   group_by(snv_id) %>%
#   summarise(supporting_cells = list(cell_id), .groups = "drop") %>%
#   deframe()


# Function to get all descendant tips (cells) for a node
get_descendant_tips <- function(tree, node) {
  if (node <= length(tree$tip.label)) {
    return(node)  # It's already a tip
  }
  descendants <- Descendants(tree, node, type = "tips")[[1]]
  return(descendants)
}

# Function to select subtrees based on size criteria  
select_subtrees <- function(tree, min_pct = 10, max_pct = 25) {
  total_tips <- length(tree$tip.label)
  min_size <- ceiling(total_tips * min_pct / 100)
  max_size <- floor(total_tips * max_pct / 100)
  
  selected_nodes <- c()
  
  # Get all internal nodes (breadth-first order)
  internal_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
  
  for (node in internal_nodes) {
    descendant_tips <- get_descendant_tips(tree, node)
    subtree_size <- length(descendant_tips)
    
    if (subtree_size >= min_size && subtree_size <= max_size) {
      # Check if not contained in previously selected subtrees
      is_contained <- any(sapply(selected_nodes, function(selected_node) {
        selected_descendants <- get_descendant_tips(tree, selected_node)
        all(descendant_tips %in% selected_descendants)
      }))
      
      if (!is_contained) {
        selected_nodes <- c(selected_nodes, node)
      }
    }
  }
  
  return(selected_nodes)
}

# Function to count SNVs supporting a subtree
count_supporting_snvs <- function(subtree_cells, snv_support_list) {
  subtree_cell_names <- tree$tip.label[subtree_cells]
  
  supported_count <- 0
  for (supporting_cells in snv_support_list) {
    # SNV supports subtree if ALL supporting cells are in the subtree
    if (all(supporting_cells %in% subtree_cell_names)) {
      supported_count <- supported_count + 1
    }
  }
  return(supported_count)
}


count_supporting_snvs_fast <- function(subtree_idx, A_list, N) {
  # N = length(tree$tip.label), passed in to avoid re-computing
  in_subtree <- logical(N)
  in_subtree[subtree_idx] <- TRUE
  
  s <- 0L
  for (A in A_list) {
    ok <- TRUE
    # early exit on first outside cell
    for (a in A) { if (!in_subtree[a]) { ok <- FALSE; break } }
    if (ok) s <- s + 1L
  }
  s
}


# Permutation test function
permutation_test <- function(subtree_cells, tree, snv_support_list, n_permutations = 500) {
  original_support <- count_supporting_snvs(subtree_cells, snv_support_list)
  
  #more_extreme <- 0
  all_cells <- 1:length(tree$tip.label)
  subtree_size <- length(subtree_cells)
  
  tip_index <- setNames(seq_along(tree$tip.label), tree$tip.label)
  # Convert each SNV’s supporting cells to integer indices
  A_list <- lapply(snv_support_list, function(cells) unname(tip_index[cells]))
  # (Optional) remove NAs if any unknown names slipped in:
  A_list <- lapply(A_list, function(v) v[!is.na(v)])
  N = length(tree$tip.label)
  
  more_extreme = lapply(1:n_permutations, function(i) {
    # Randomly sample cells of same size
    permuted_cells <- sample(all_cells, subtree_size)
    permuted_support <- count_supporting_snvs_fast(permuted_cells, A_list, N)
    return(permuted_support >= original_support)
  })
  
  more_extreme = sum(unlist(more_extreme))
  
  p_value <- more_extreme / n_permutations
  return(list(p_value = p_value, original_support = original_support))
}

# Main validation function for a single sample
validate_single_sample <- function(tree, snv_df, min_pct = 5, max_pct = 15, n_permutations = 500) {
  
  # Prepare SNV support data
  snv_support_list <- snv_df %>%
    filter(alt_count > 0) %>%  # Only consider SNVs with supporting reads
    mutate(snv_id = paste(chr, start, sep = "_")) %>%
    group_by(snv_id) %>%
    summarise(supporting_cells = list(unique(cell_id)), .groups = "drop") %>%
    pull(supporting_cells)
  
  cat("Found", length(snv_support_list), "SNVs with supporting reads\n")
  
  # Select subtrees for analysis
  
  clades = find_min_clades(tree, min_frac = 0.05)
  selected_subtrees = clades$nodes
  
  #selected_subtrees <- select_subtrees(tree, min_pct, max_pct)
  cat("Selected", length(selected_subtrees), "subtrees for analysis\n")
  
  # Analyze each subtree
  results <- data.frame()
  
  for (i in seq_along(selected_subtrees)) {
    node <- selected_subtrees[i]
    subtree_tips <- get_descendant_tips(tree, node)
    
    cat("Analyzing subtree", i, "with", length(subtree_tips), "cells\n")
    
    # Perform permutation test
    test_result <- permutation_test(subtree_cells = subtree_tips, 
                                    tree = tree, 
                                    snv_support_list = snv_support_list, 
                                    n_permutations = n_permutations)
    
    results <- rbind(results, data.frame(
      subtree_id = i,
      node = node,
      subtree_size = length(subtree_tips),
      subtree_percentage = round(length(subtree_tips) / length(tree$tip.label) * 100, 1),
      supporting_snvs = test_result$original_support,
      p_value = test_result$p_value,
      significant = test_result$p_value < 0.05
    ))
  }
  
  return(results)
}

find_min_clades <- function(tree, min_frac = 0.10, min_size = NULL) {
  tree <- reorder.phylo(tree, "postorder")
  N <- length(tree$tip.label)
  if (is.null(min_size)) min_size <- ceiling(min_frac * N)
  
  parent <- tree$edge[,1]
  child  <- tree$edge[,2]
  children_list <- split(child, parent)
  
  total_nodes <- N + tree$Nnode
  # subtree tip counts
  cnt <- integer(total_nodes)
  cnt[1:N] <- 1L
  # postorder accumulation: children appear before parent in "postorder"
  for (i in seq_len(nrow(tree$edge))) {
    cnt[parent[i]] <- cnt[parent[i]] + cnt[child[i]]
  }
  
  # inclusion-minimal nodes whose subtree size >= min_size
  elig <- which(cnt >= min_size)
  sel <- integer(0)
  for (v in elig) {
    ch <- children_list[[as.character(v)]]
    # internal nodes only qualify if all children are below threshold
    # (so you can't push the cut further down)
    if (length(ch) == 0L) next
    if (all(cnt[ch] < min_size)) sel <- c(sel, v)
  }
  
  # build tip groups for selected nodes
  clades <- lapply(sel, function(v) tree$tip.label[Descendants(tree, v, "tips")[[1]]])
  names(clades) <- paste0("clade_", seq_along(clades))
  
  # tip → clade mapping (each tip belongs to exactly one selected clade)
  tip2clade <- setNames(rep(NA_character_, N), tree$tip.label)
  for (nm in names(clades)) tip2clade[clades[[nm]]] <- nm
  
  list(
    min_size      = min_size,
    nodes         = sel,            # node indices (N+1 .. N+Nnode)
    clades        = clades,         # named list of tip.label vectors
    tip_clusters  = tip2clade       # named vector, length N
  )
}

SNV_analyis = function(s, snv_all) {
  # Bridges first
  tree = readRDS(paste0("results/bridges_trees/",s,".rds"))$tree
  snv_df = snv_all[snv_all$patient == s,] %>% 
    dplyr::rename(alt_count = alt_counts)
  keep = intersect(snv_df$cell_id, tree$tip.label)
  tree <- keep.tip(tree, keep)
  snv_df = snv_df %>% dplyr::filter(cell_id %in% keep)
  bridges_res = validate_single_sample(tree, snv_df, min_pct = 10, max_pct = 100)
  
  
  
  # Sitka second
  tree = ape::read.tree(file = paste0("signatures_dataset/sitka_trees/",s,"-cn-tree.newick"))
  snv_df = snv_all[snv_all$patient == s,] %>% 
    dplyr::rename(alt_count = alt_counts)
  snv_df$cell_id = paste0("cell_", snv_df$cell_id)
  keep = intersect(snv_df$cell_id, tree$tip.label)
  tree <- keep.tip(tree, keep)
  snv_df = snv_df %>% dplyr::filter(cell_id %in% keep)
  sitka_res = validate_single_sample(tree, snv_df)
}
