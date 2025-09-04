
rm(list = ls())
library(data.table)
library(tidyverse)
library(ape)
library(phangorn)
source("getters.R")

snv_annotations = fread("signatures_dataset/DLP/SNV/snv_annotations.csv.gz")
snv_all = fread("signatures_dataset/DLP/SNV/snv_counts.csv.gz")

# ---- Packages ----
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Matrix)     # for sparse matrices and fast crossprod
  library(data.table) # for fast SNV -> matrix conversion
})

# ---- 1) SNVs (long) -> binary Y (cells x loci) ----
snv_to_binary_Y <- function(snv_df,
                            cells_in_order,
                            min_alt = 1L,
                            min_cells_with_alt = 2L,
                            min_cells_without_alt = 2L) {
  dt <- as.data.table(snv_df)
  setnames(dt, tolower(names(dt)))
  req <- c("cell_id","chromosome","position","alt_count")
  if (!all(req %in% names(dt))) stop("snv_df must have columns: ", paste(req, collapse=", "))
  
  # only cells that are in the tree
  dt <- dt[cell_id %in% cells_in_order]
  
  # define locus id
  dt[, locus := paste0(chromosome, ":", position)]
  # binary presence by threshold
  dt[, present := as.integer(!is.na(alt_count) & alt_count >= min_alt)]
  dt <- dt[present == 1L, .(cell_id, locus)]   # keep positives only
  
  # filter loci by prevalence (optional but recommended)
  # count positives per locus; we need total cells C to evaluate negatives
  C <- length(cells_in_order)
  if (C < 2L) stop("Need at least 2 cells after intersecting with the tree.")
  
  pos_counts <- dt[, .N, by = locus]
  keep_loci <- pos_counts[N >= min_cells_with_alt & (C - N) >= min_cells_without_alt, locus]
  if (length(keep_loci) == 0L) stop("No loci passed prevalence filters.")
  
  dt <- dt[locus %in% keep_loci]
  
  # build sparse Y: rows = cells_in_order, cols = loci (alphabetical, stable)
  loci <- sort(unique(dt$locus))
  row_index <- match(dt$cell_id, cells_in_order)
  col_index <- match(dt$locus, loci)
  
  Y <- sparseMatrix(i = row_index,
                    j = col_index,
                    x = 1L,
                    dims = c(C, length(loci)),
                    dimnames = list(cells_in_order, loci))
  
  Y
}

# ---- 2) Tree -> clade-membership G (cells x K) ----
build_clade_matrix <- function(tree) {
  if (!is.rooted(tree)) tree <- midpoint(tree)
  # tips order is the row order
  cells <- tree$tip.label
  innodes <- (Ntip(tree) + 1L):(Ntip(tree) + tree$Nnode)
  D <- phangorn::Descendants(tree, innodes, type = "tips")
  
  # build sparse matrix: one column per internal node; ones at descendant tip rows
  C <- length(cells); K <- length(D)
  if (K == 0L) stop("Tree has no internal nodes.")
  
  # pack i,j for all columns
  i <- integer(0); j <- integer(0)
  for (k in seq_len(K)) {
    tips_k <- as.integer(D[[k]])
    i <- c(i, tips_k)
    j <- c(j, rep.int(k, length(tips_k)))
  }
  G <- sparseMatrix(i = i, j = j, x = 1L,
                    dims = c(C, K),
                    dimnames = list(cells, paste0("node", innodes)))
  list(G = G, tree = tree) # return possibly re-rooted tree too
}

# ---- 3) Youden J from Y and tree ----
youden_from_tree_and_Y <- function(tree, Y, tie_break = c("first","smallest","random")) {
  tie_break <- match.arg(tie_break)
  
  # keep only cells present in both Y and tree, and align orders
  common <- intersect(tree$tip.label, rownames(Y))
  tree <- keep.tip(tree, common)
  # after keep.tip, tip order changed; reindex Y rows to that order
  Y <- Y[tree$tip.label, , drop = FALSE]
  
  G_pack <- build_clade_matrix(tree)
  G <- G_pack$G  # cells x K
  C <- nrow(G); L <- ncol(Y); K <- ncol(G)
  if (nrow(Y) != C) stop("Row mismatch after alignment.")
  
  # Precompute counts
  s1 <- Matrix::colSums(G)   # predicted positives per clade
  r1 <- Matrix::colSums(Y)   # positives per locus
  
  # Overlaps for all clades x all loci: V[k,l] = sum_c G[c,k]*Y[c,l]
  V <- Matrix::crossprod(G, Y)  # K x L (sparse-friendly)
  
  # Score to maximize accuracy: 2*TP - PredPos
  # recycle s1 (K vector) across columns
  score <- 2 * V - s1
  
  # argmax per locus with tie-breaking
  # get column-wise max indices
  vmax <- apply(score, 2L, max)
  # logical KxL of maxima (sparse-friendly way: compare to vmax per column)
  # we'll just loop columns to honor tie rules deterministically
  pick <- integer(L)
  for (l in seq_len(L)) {
    col <- score[, l]
    wmax <- which(col == vmax[l])
    if (length(wmax) == 1L) {
      pick[l] <- wmax
    } else {
      if (tie_break == "first") {
        pick[l] <- wmax[1L]
      } else if (tie_break == "smallest") {
        pick[l] <- wmax[which.min(s1[wmax])]
      } else {
        pick[l] <- sample(wmax, 1L)
      }
    }
  }
  
  # Extract per-locus TP, then FP,FN,TN
  tp <- as.numeric(V[cbind(pick, seq_len(L))])
  fp <- as.numeric(s1[pick] - tp)
  fn <- as.numeric(r1        - tp)
  tn <- as.numeric(C - (tp + fp + fn))
  
  # Aggregate confusion
  H <- c(TN = sum(tn), FP = sum(fp), FN = sum(fn), TP = sum(tp))
  sens <- if ((H["TP"] + H["FN"]) > 0) H["TP"] / (H["TP"] + H["FN"]) else 1
  spec <- if ((H["TN"] + H["FP"]) > 0) H["TN"] / (H["TN"] + H["FP"]) else 1
  J <- sens + spec - 1
  
  # Delta-method CI using per-locus fractions
  Cn <- C
  P <- cbind(tn, fp, fn, tp) / Cn   # L x 4
  pbar <- colMeans(P)
  den_spec <- pbar[1] + pbar[2]
  den_sens <- pbar[4] + pbar[3]
  eps <- .Machine$double.eps
  
  gvec <- c(  pbar[2] / (den_spec + eps)^2,     # d/d TN
              -pbar[1] / (den_spec + eps)^2,     # d/d FP
              -pbar[4] / (den_sens + eps)^2,     # d/d FN
              pbar[3] / (den_sens + eps)^2)     # d/d TP
  
  S <- stats::cov(P) # 4x4
  varJ <- as.numeric(t(gvec) %*% S %*% gvec) / nrow(P)
  seJ  <- sqrt(max(varJ, 0))
  ci95 <- c(lower = J - 1.96 * seJ, upper = J + 1.96 * seJ)
  
  list(J = J, ci = ci95, H = H,
       picks = data.frame(locus = colnames(Y), picked_clade = colnames(G)[pick],
                          tp = tp, fp = fp, fn = fn, tn = tn, row.names = NULL))
}

# ---- 4) End-to-end convenience wrapper ----
tree_gof_youden <- function(tree,
                            snv_df,
                            min_alt = 1L,
                            min_cells_with_alt = 2L,
                            min_cells_without_alt = 2L,
                            tie_break = c("first","smallest","random")) {
  if (!is.rooted(tree)) tree <- midpoint(tree)
  # keep overlap cells first, then build Y in tree order
  common_cells <- intersect(tree$tip.label, unique(snv_df$cell_id))
  tree <- keep.tip(tree, common_cells)
  cells <- tree$tip.label
  
  Y <- snv_to_binary_Y(
    snv_df,
    cells_in_order = cells,
    min_alt = min_alt,
    min_cells_with_alt = min_cells_with_alt,
    min_cells_without_alt = min_cells_without_alt
  )
  
  youden_from_tree_and_Y(tree, Y, tie_break = match.arg(tie_break))
}

get_youden_index = function(s) {
  print(s)
  # Bridges first
  tree = readRDS(paste0("results/bridges_trees/",s,".rds"))$tree
  snv_df = snv_all[snv_all$patient == s,] %>% 
    dplyr::rename(alt_count = alt_counts)
  keep = intersect(snv_df$cell_id, tree$tip.label)
  tree <- keep.tip(tree, keep)
  snv_df = snv_df %>% dplyr::filter(cell_id %in% keep)
  snv_df = snv_df %>% dplyr::select(cell_id, chr, start, alt_count) %>% 
    dplyr::rename(chromosome = chr, position = start)
  
  res_bridges <- tree_gof_youden(
    tree,
    snv_df,
    min_alt = 2L,                 # presence if alt_count >= 1
    min_cells_with_alt = 3L,      # filter very rare loci
    min_cells_without_alt = 3L,   # filter ubiquitous loci
    tie_break = "smallest"         # prefer smallest clade when tied
  )
  
  # Lazac
  if (file.exists(paste0("signatures_dataset/lazac_trees/",s,"_hscn_tree.newick"))) {
    newick_text = read_file(paste0("signatures_dataset/lazac_trees/",s,"_hscn_tree.newick"))
    newick_text = paste0(newick_text, ";")
    tree = ape::read.tree(text = newick_text)
    snv_df = snv_all[snv_all$patient == s,] %>% 
      dplyr::rename(alt_count = alt_counts)
    keep = intersect(snv_df$cell_id, tree$tip.label)
    tree <- keep.tip(tree, keep)
    snv_df = snv_df %>% dplyr::filter(cell_id %in% keep)
    snv_df = snv_df %>% dplyr::select(cell_id, chr, start, alt_count) %>% 
      dplyr::rename(chromosome = chr, position = start)
    
    res_lazac <- tree_gof_youden(
      tree,
      snv_df,
      min_alt = 2L,                 # presence if alt_count >= 1
      min_cells_with_alt = 3L,      # filter very rare loci
      min_cells_without_alt = 3L,   # filter ubiquitous loci
      tie_break = "smallest"         # prefer smallest clade when tied
    )
  } else {
    res_lazac = list(J = NA)
  }
  
  # Sitka second
  tree = ape::read.tree(file = paste0("signatures_dataset/sitka_trees/",s,"-cn-tree.newick"))
  snv_df = snv_all[snv_all$patient == s,] %>% 
    dplyr::rename(alt_count = alt_counts)
  snv_df$cell_id = paste0("cell_", snv_df$cell_id)
  keep = intersect(snv_df$cell_id, tree$tip.label)
  tree <- keep.tip(tree, keep)
  snv_df = snv_df %>% dplyr::filter(cell_id %in% keep)
  snv_df = snv_df %>% dplyr::select(cell_id, chr, start, alt_count) %>% 
    dplyr::rename(chromosome = chr, position = start)
  
  res_sitka <- tree_gof_youden(
    tree,
    snv_df,
    min_alt = 2L,                 # presence if alt_count >= 1
    min_cells_with_alt = 3L,      # filter very rare loci
    min_cells_without_alt = 3L,   # filter ubiquitous loci
    tie_break = "smallest"         # prefer smallest clade when tied
  )

  dplyr::bind_rows(
    dplyr::tibble(Jouden = res_bridges$J, mehtod = "bridges"),
    dplyr::tibble(Jouden = res_sitka$J, mehtod = "sitka"),
    dplyr::tibble(Jouden = res_lazac$J, mehtod = "lazac")
  ) %>% dplyr::mutate(sample_id = s)
}

samples = get_sample_names()
youden_df = lapply(samples, function(s) get_youden_index(s))
youden_df = do.call("rbind", youden_df)

saveRDS(youden_df, "results/metrics/youden.rds")
