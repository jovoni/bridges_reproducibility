
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

# 1.) Clonal discordance ####
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
  
  # Lazac
  labels_df = get_labels(sample_id)
  
  if (file.exists(paste0("signatures_dataset/lazac_trees/",sample_id,"_hscn_tree.newick"))) {
    newick_text = read_file(paste0("signatures_dataset/lazac_trees/",sample_id,"_hscn_tree.newick"))
    newick_text = paste0(newick_text, ";")
    tree = ape::read.tree(text = newick_text)
    keep = intersect(labels_df$cell_id, tree$tip.label)
    tree_pruned <- keep.tip(tree, keep)
    #Xdiss = get_X_for_dissimilarity(sample_id)
    labels_pruned <- labels_df %>%
      dplyr::filter(cell_id %in% keep)
    d3 = dplyr::tibble(sample_id = sample_id, value = clonal_discordance_from_df(tree = tree_pruned, labels_pruned, mode), method = "lazac")
    
    return(dplyr::bind_rows(d1, d2, d3))
  }
  
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
# 2.) Mean sibling dissimilarity ####
# X: matrix/data.frame of features per tip (rows = tips, cols = features)
#     Works best with discrete features; for numeric we use exact match unless tol>0.
# tol: numeric tolerance for treating two numeric entries as "equal".
# na.rm: ignore NA positions in the per-pair comparison.
# sibling_dissimilarity_discrete_rcpp.R
# Fast sibling dissimilarity for discrete features using Rcpp in a single R file.
# Dependencies: Rcpp (for compilation), ape (for keep.tip).
# Usage:
#   source("sibling_dissimilarity_discrete_rcpp.R")
#   d <- sibling_dissimilarity_discrete_rcpp(tree, X, na.rm = TRUE)

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required. Please install it first.")
}

# Compile the C++ kernel when this file is sourced
Rcpp::cppFunction(code = '
#include <Rcpp.h>
using namespace Rcpp;

// Compute normalized Hamming distance for each pair of rows in X.
// X: IntegerMatrix (rows = tips in tree order, cols = features), can contain NA_INTEGER
// pairs: IntegerMatrix with 2 columns, 1-based indices into rows of X
// na_rm: if true, positions with NA in either row are skipped; if all skipped -> NA_real_
//        if false, any NA in compared positions makes result NA_real_
// [[Rcpp::export]]
NumericVector pair_hamming_discrete_cpp(IntegerMatrix X, IntegerMatrix pairs, bool na_rm = true) {
  const int P = pairs.nrow();
  const int M = X.ncol();
  NumericVector out(P);

  for (int k = 0; k < P; ++k) {
    int i = pairs(k, 0) - 1; // convert to 0-based
    int j = pairs(k, 1) - 1;

    int keep = 0, diff = 0;
    bool any_na = false;

    for (int m = 0; m < M; ++m) {
      int a = X(i, m);
      int b = X(j, m);

      if (a == NA_INTEGER || b == NA_INTEGER) {
        if (!na_rm) { any_na = true; break; }
        continue;
      }
      ++keep;
      if (a != b) ++diff;
    }

    if (!na_rm && any_na) {
      out[k] = NA_REAL;
    } else if (keep == 0) {
      out[k] = NA_REAL;
    } else {
      out[k] = static_cast<double>(diff) / static_cast<double>(keep);
    }
  }

  return out;
}
')

#' Sibling dissimilarity for discrete features (fast, Rcpp-backed)
#'
#' Computes normalized Hamming distances between all sibling tip pairs
#' (tips that share the same immediate parent) in a possibly non-binary tree.
#' Works with **discrete** features only; columns are coerced to integer codes once.
#'
#' @param tree an \code{ape::phylo} object
#' @param X data.frame or matrix with rownames = tip labels, cols = discrete features
#' @param na.rm logical; if TRUE (default), positions with NA in either row are ignored.
#'   If all positions are NA, returns NA for that pair. If FALSE, any NA yields NA.
#' @return numeric vector of sibling distances (one per sibling pair).
#'         Attributes \code{pairs_idx} (integer matrix of row indices used)
#'         and \code{pairs_labels} (character matrix of tip labels) are attached for convenience.
#' @examples
#' # d <- sibling_dissimilarity_discrete_rcpp(tree, X)
sibling_dissimilarity_discrete_rcpp <- function(tree, X, na.rm = TRUE) {
  # 1) Align X to tree tips
  tips <- tree$tip.label
  X <- as.data.frame(X, check.names = FALSE)
  if (is.null(rownames(X))) stop("X must have rownames matching tree tip labels.")
  common <- intersect(tips, rownames(X))
  if (length(common) < 2L) stop("Not enough overlapping tips and X rows.")
  tree <- ape::keep.tip(tree, common)
  X <- X[tree$tip.label, , drop = FALSE]
  
  # 2) Coerce every column to integer codes (preserving NA)
  Xint_list <- lapply(X, function(col) {
    if (is.integer(col)) return(col)
    if (is.factor(col)) return(as.integer(col))
    # character/logical/other → factor → integer codes; exclude = NULL preserves NA
    as.integer(factor(col, exclude = NULL))
  })
  Xint <- as.matrix(as.data.frame(Xint_list, check.names = FALSE))
  storage.mode(Xint) <- "integer"
  
  # 3) Build sibling pairs (all unordered pairs of direct tip-children per parent)
  Ntip <- length(tree$tip.label)
  parents  <- tree$edge[, 1]
  children <- tree$edge[, 2]
  child_list <- split(children, parents)
  
  pair_list <- vector("list", length(child_list))
  n_added <- 0L
  it <- 0L
  for (kids in child_list) {
    it <- it + 1L
    tip_kids <- kids[kids <= Ntip]
    if (length(tip_kids) >= 2L) {
      n_added <- n_added + 1L
      pair_list[[n_added]] <- t(combn(tip_kids, 2L))
    }
  }
  if (n_added == 0L) stop("No sibling tip pairs found (tree may be too unresolved).")
  if (n_added < length(pair_list)) pair_list <- pair_list[seq_len(n_added)]
  pairs_idx <- do.call(rbind, pair_list)
  mode(pairs_idx) <- "integer"
  
  # 4) Compute distances with the C++ kernel
  d <- pair_hamming_discrete_cpp(Xint, pairs_idx, na_rm = isTRUE(na.rm))
  
  # Attach labels for convenience
  pairs_labels <- cbind(tree$tip.label[pairs_idx[, 1]], tree$tip.label[pairs_idx[, 2]])
  colnames(pairs_labels) <- c("tip1", "tip2")
  attr(d, "pairs_idx") <- pairs_idx
  attr(d, "pairs_labels") <- pairs_labels
  d
}

compute_sibling_similarities = function(sample_id) {
  # Get Xdiss
  Xdiss = get_X_for_dissimilarity(sample_id)
  
  # Read bridges fit
  fit = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))
  tree = fit$tree
  keep = intersect(rownames(Xdiss), tree$tip.label)
  tree_pruned <- keep.tip(tree, keep)
  bridges_diss = sibling_dissimilarity_discrete_rcpp(tree_pruned, Xdiss)
  
  # Read lazac fit
  if (file.exists(paste0("signatures_dataset/lazac_trees/",sample_id,"_hscn_tree.newick"))) {
    newick_text = read_file(paste0("signatures_dataset/lazac_trees/",sample_id,"_hscn_tree.newick"))
    newick_text = paste0(newick_text, ";")
    tree = ape::read.tree(text = newick_text)
    keep = intersect(rownames(Xdiss), tree$tip.label)
    tree_pruned <- keep.tip(tree, keep)
    lazac_diss = sibling_dissimilarity_discrete_rcpp(tree_pruned, Xdiss)  
  } else {
    lazac_diss = NULL
  }
  
  # Read sitka tree
  tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
  rownames(Xdiss) = paste0("cell_", rownames(Xdiss))
  keep = intersect(rownames(Xdiss), tree$tip.label)
  tree_pruned <- keep.tip(tree, keep)
  sitka_diss = sibling_dissimilarity_discrete_rcpp(tree_pruned, Xdiss)
  
  dplyr::bind_rows(
    dplyr::tibble(name = "sitka", diss = sitka_diss),
    dplyr::tibble(name = "bridges", diss = bridges_diss),
  ) %>% dplyr::mutate(sample_id = sample_id)
}

# 3.) Clone F1 against best clade ####
# # labels_df: data.frame(cell_id, clone_id)
# # Returns per-clone F1 plus macro and size-weighted averages.
clone_f1_against_best_clade <- function(tree, labels_df) {
  labs <- labels_df %>% distinct(cell_id, clone_id)
  common <- intersect(tree$tip.label, labs$cell_id)
  if (length(common) < 2L) stop("Not enough overlapping tips and labels.")
  tree <- keep.tip(tree, common)
  labs <- labs %>% filter(cell_id %in% common)

  # all clades = descendant tip sets for each internal node
  Ntip <- length(tree$tip.label)
  internals <- (Ntip + 1):(Ntip + tree$Nnode)
  clades <- lapply(internals, function(nd) {
    tips_idx <- Descendants(tree, nd, type = "tips")[[1]]
    tree$tip.label[tips_idx]
  })
  # remove duplicates (can happen in polytomies)
  key <- vapply(clades, function(x) paste(sort(x), collapse = "|"), "")
  clades <- clades[!duplicated(key)]

  # helper: F1 between two tip sets
  f1_set <- function(A, B) {
    TP <- length(intersect(A, B))
    if (TP == 0) return(0)
    FP <- length(setdiff(B, A))
    FN <- length(setdiff(A, B))
    2 * TP / (2 * TP + FP + FN)
  }

  # compute per-clone best F1
  clones <- split(labs$cell_id, labs$clone_id)
  per_clone <- lapply(names(clones), function(cl) {
    A <- clones[[cl]]
    best <- max(vapply(clades, function(B) f1_set(A, B), numeric(1)))
    data.frame(clone_id = cl, size = length(A), best_F1 = best, row.names = NULL)
  }) %>% bind_rows()

  macro_F1 <- mean(per_clone$best_F1)
  weighted_F1 <- with(per_clone, sum(best_F1 * size) / sum(size))

  list(per_clone = per_clone,
       macro_F1 = macro_F1,
       size_weighted_F1 = weighted_F1)
}

compute_best_f1 = function(sample_id) {
  # Get labels
  labels = get_labels(sample_id)
  
  # Read bridges fit
  fit = readRDS(paste0("results/bridges_trees/", sample_id, ".rds"))
  tree = fit$tree
  keep = intersect(labels$cell_id, tree$tip.label)
  tree_pruned <- keep.tip(tree, keep)
  labels_pruned = labels %>% dplyr::filter(cell_id %in% keep)
  f1_bridges = clone_f1_against_best_clade(tree_pruned, labels_pruned)
  
  # Read sitka tree
  tree = ape::read.tree(paste0("signatures_dataset/sitka_trees/",sample_id,"-cn-tree.newick"))
  labels$cell_id = paste0("cell_", labels$cell_id)
  keep = intersect(labels$cell_id, tree$tip.label)
  tree_pruned <- keep.tip(tree, keep)
  labels_pruned = labels %>% dplyr::filter(cell_id %in% keep)
  f1_sitka = clone_f1_against_best_clade(tree_pruned, labels_pruned)
  
  
  list(
    f1_bridges = f1_bridges,
    f1_sitka = f1_sitka
  )
}