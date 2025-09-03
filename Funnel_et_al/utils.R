
source("getters.R")
library(ape)
library(phangorn)

# Cost matrix for Sankoff parsimony
.make_cost <- function(n_states, mode = c("fitch", "wagner")) {
  mode <- match.arg(mode)
  if (mode == "fitch") {
    M <- matrix(1, n_states, n_states); diag(M) <- 0
  } else {
    i <- seq_len(n_states); M <- abs(outer(i, i, "-"))
  }
  M
}

clonal_discordance_from_df <- function(tree, labels_df,
                                       mode = c("fitch", "wagner"),
                                       tie_break = c("first", "random"),
                                       seed = 1) {
  
  stopifnot(inherits(tree, "phylo"))
  stopifnot(all(c("cell_id", "clone_id") %in% names(labels_df)))
  mode <- match.arg(mode)
  tie_break <- match.arg(tie_break)
  
  # Keep only labels for tips present in the tree, aligned to tree$tip.label
  lab_map <- setNames(labels_df$clone_id, labels_df$cell_id)
  missing_tips <- setdiff(tree$tip.label, names(lab_map))
  if (length(missing_tips)) {
    stop("These tips have no clone_id in labels_df: ",
         paste(head(missing_tips, 10), collapse = ", "),
         if (length(missing_tips) > 10) " ...")
  }
  tip_labels_chr <- lab_map[tree$tip.label]
  
  # Encode as factor states (levels preserve unique order of appearance)
  states <- factor(tip_labels_chr, levels = unique(tip_labels_chr))
  state_levels <- levels(states)
  K <- length(state_levels)
  
  # Build phyDat with one character site named "clone"
  dat <- phyDat(matrix(states, ncol = 1,
                       dimnames = list(tree$tip.label, "clone")),
                type = "USER", levels = state_levels)
  
  # Sankoff parsimony score with chosen cost matrix
  cost <- .make_cost(K, mode)
  score <- as.numeric(parsimony(tree, dat, method = "sankoff", cost = cost))
  score
}


clonal_discordance = function(sample_id, mode) {
  fit = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))
  labels_df = get_labels(sample_id)
  keep = intersect(labels_df$cell_id, fit$tree$tip.label)
  tree_pruned <- keep.tip(fit$tree, keep)
  #Xdiss = get_X_for_dissimilarity(sample_id)
  labels_pruned <- labels_df %>%
    dplyr::filter(cell_id %in% keep)
  d1 = dplyr::tibble(sample_id = sample_id, value = clonal_discordance_from_df(tree = tree_pruned, labels_df, mode), method = "bridges")
  
  labels_df = get_labels(sample_id)
  sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
  labels_df_sitka = labels_df %>% dplyr::mutate(cell_id = paste0("cell_", cell_id))
  keep = intersect(labels_df_sitka$cell_id, sitka_tree$tip.label)
  tree_pruned <- keep.tip(sitka_tree, keep)
  labels_pruned <- labels_df_sitka %>%
    dplyr::filter(cell_id %in% keep)

  d2 = dplyr::tibble(sample_id = sample_id, value = clonal_discordance_from_df(tree = tree_pruned, labels_pruned, mode), method = "sitka")
  dplyr::bind_rows(d1, d2)
}

bridges_score = function(sample_id, mode) {
  fit = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))
  labels_df = get_labels(sample_id)
  keep = intersect(labels_df$cell_id, fit$tree$tip.label)
  tree_pruned <- keep.tip(fit$tree, keep)
  Xdiss = get_X_for_dissimilarity(sample_id)
  labels_pruned <- labels_df %>%
    dplyr::filter(cell_id %in% keep)
  
  list(
    clonal_discordance = clonal_discordance_from_df(tree = tree_pruned, labels_df, mode),
    clone_f1 = clone_f1_against_best_clade(tree_pruned, labels_pruned),
    sibling_diss = sibling_dissimilarity(tree_pruned, Xdiss)
  )
}

sitka_score = function(sample_id, mode) {
  labels_df = get_labels(sample_id)
  sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
  labels_df_sitka = labels_df %>% dplyr::mutate(cell_id = paste0("cell_", cell_id))
  #keep <- labels_df_sitka$cell_id
  keep = intersect(labels_df_sitka$cell_id, sitka_tree$tip.label)
  tree_pruned <- keep.tip(sitka_tree, keep)
  labels_pruned <- labels_df_sitka %>%
    dplyr::filter(cell_id %in% keep)
  Xdiss = get_X_for_dissimilarity(sample_id)
  rownames(Xdiss) = paste0("cell_", rownames(Xdiss))
  Xdiss = Xdiss[rownames(Xdiss) %in% keep,]
  
  list(
    clonal_discordance = clonal_discordance_from_df(tree = tree_pruned, labels_pruned, mode),
    clone_f1 = clone_f1_against_best_clade(tree_pruned, labels_pruned),
    sibling_diss = sibling_dissimilarity(tree_pruned, Xdiss)
  )
}

# library(ape)
# library(dplyr)
# 
# # ---- (A) Mean sibling dissimilarity -----------------------------------------
# # X: matrix/data.frame of features per tip (rows = tips, cols = features)
# #     Works best with discrete features; for numeric we use exact match unless tol>0.
# # tol: numeric tolerance for treating two numeric entries as "equal".
# # na.rm: ignore NA positions in the per-pair comparison.
# sibling_dissimilarity <- function(tree, X, tol = 0, na.rm = TRUE) {
#   # keep common tips and order rows as tree tips
#   tips <- tree$tip.label
#   X <- as.data.frame(X, check.names = FALSE)
#   common <- intersect(tips, rownames(X))
#   if (length(common) < 2L) stop("Not enough overlapping tips and X rows.")
#   tree <- keep.tip(tree, common)
#   X <- X[tree$tip.label, , drop = FALSE]
#   
#   # per-entry equality function (handles numeric (with tol) and non-numeric)
#   eq_entry <- function(a, b) {
#     if (is.numeric(a) && is.numeric(b)) {
#       abs(a - b) <= tol
#     } else {
#       a == b
#     }
#   }
#   
#   # vectorized normalized Hamming distance between two rows
#   nhamming <- function(i, j) {
#     a <- as.list(X[i, , drop = TRUE])
#     b <- as.list(X[j, , drop = TRUE])
#     same <- mapply(eq_entry, a, b)
#     if (na.rm) {
#       keep <- !(is.na(a) | is.na(b))
#       if (!any(keep)) return(NA_real_)
#       same <- same[keep]
#     }
#     mean(!same)
#   }
#   
#   # build sibling pairs: tips sharing the same immediate parent
#   Ntip <- length(tree$tip.label)
#   parents <- tree$edge[, 1]; children <- tree$edge[, 2]
#   # children of each internal node
#   child_list <- split(children, parents)
#   # for each parent, take its direct tip-children; all unordered pairs
#   pairs <- list()
#   for (p in names(child_list)) {
#     kids <- child_list[[p]]
#     tip_kids <- kids[kids <= Ntip]
#     if (length(tip_kids) >= 2L) {
#       comb <- t(combn(tip_kids, 2L))
#       pairs[[length(pairs) + 1L]] <- comb
#     }
#   }
#   if (length(pairs) == 0L) stop("No sibling tip pairs found (tree may be too unresolved).")
#   pairs <- do.call(rbind, pairs)
#   
#   # map node indices (1..Ntip) to row indices in X
#   idx_to_row <- seq_len(Ntip)  # already ordered by tree$tip.label
#   dists <- apply(pairs, 1L, function(pr) nhamming(idx_to_row[pr[1]], idx_to_row[pr[2]]))
#   
#   #mean(dists, na.rm = TRUE)
#   dists
# }
# 
# # ---- (B) Clone F1 against best clade ----------------------------------------
# # labels_df: data.frame(cell_id, clone_id)
# # Returns per-clone F1 plus macro and size-weighted averages.
# clone_f1_against_best_clade <- function(tree, labels_df) {
#   labs <- labels_df %>% distinct(cell_id, clone_id)
#   common <- intersect(tree$tip.label, labs$cell_id)
#   if (length(common) < 2L) stop("Not enough overlapping tips and labels.")
#   tree <- keep.tip(tree, common)
#   labs <- labs %>% filter(cell_id %in% common)
#   
#   # all clades = descendant tip sets for each internal node
#   Ntip <- length(tree$tip.label)
#   internals <- (Ntip + 1):(Ntip + tree$Nnode)
#   clades <- lapply(internals, function(nd) {
#     tips_idx <- Descendants(tree, nd, type = "tips")[[1]]
#     tree$tip.label[tips_idx]
#   })
#   # remove duplicates (can happen in polytomies)
#   key <- vapply(clades, function(x) paste(sort(x), collapse = "|"), "")
#   clades <- clades[!duplicated(key)]
#   
#   # helper: F1 between two tip sets
#   f1_set <- function(A, B) {
#     TP <- length(intersect(A, B))
#     if (TP == 0) return(0)
#     FP <- length(setdiff(B, A))
#     FN <- length(setdiff(A, B))
#     2 * TP / (2 * TP + FP + FN)
#   }
#   
#   # compute per-clone best F1
#   clones <- split(labs$cell_id, labs$clone_id)
#   per_clone <- lapply(names(clones), function(cl) {
#     A <- clones[[cl]]
#     best <- max(vapply(clades, function(B) f1_set(A, B), numeric(1)))
#     data.frame(clone_id = cl, size = length(A), best_F1 = best, row.names = NULL)
#   }) %>% bind_rows()
#   
#   macro_F1 <- mean(per_clone$best_F1)
#   weighted_F1 <- with(per_clone, sum(best_F1 * size) / sum(size))
#   
#   list(per_clone = per_clone,
#        macro_F1 = macro_F1,
#        size_weighted_F1 = weighted_F1)
# }
# 
# 
# library(ape)
# library(phangorn)
# library(dplyr)
# 
# sample_id = "SA1292"
# 
# Xdiss = get_X_for_dissimilarity(sample_id, allele = "CN")
# 
# # Bridges
# 
# tree = fit$tree
# labels_df = get_labels(sample_id)
# 
# bridges_score(fit, sample_id, "fitch")
# clone_f1_against_best_clade(tree, labels_df)
# sibling_diss = sibling_dissimilarity(tree, Xdiss)
# hist(sibling_diss)
# mean(sibling_diss)
# 
# sitka_tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
# labels_df_sitka = labels_df %>% dplyr::mutate(cell_id = paste0("cell_", cell_id))
# tree = sitka_tree
# labels_df = labels_df_sitka
# 
# sitka_score(sample_id, mode = "fitch")
# rownames(Xdiss) = paste0("cell_", rownames(Xdiss))
# clone_f1_against_best_clade(tree, labels_df)
# 
# sibling_diss = sibling_dissimilarity(tree, Xdiss)
# hist(sibling_diss)
# mean(sibling_diss)
