
snv_df = fread("signatures_dataset/DLP/SNV/snv_counts.csv.gz") %>% 
  dplyr::select(cell_id, chr, start, alt_counts, patient)

s = "SA1035"

tree = readRDS("results/bridges_trees/SA1035.rds")$tree
snv_dt = snv_df[snv_df$patient == s,] %>% 
  dplyr::rename(chrom = chr, pos = start, alt_count = alt_counts)

library(ape)
library(data.table)

# Inputs per sample:
# tree: ape::phylo (rooted; tips are cell IDs)
# snv_dt: data.table with columns: cell_id, chrom, pos, alt_count [, total_count, etc.]
# params:
alpha <- 0.1; beta <- 0.25; B <- 500; alt_min <- 2L

# --- 1) Harmonize
tips <- tree$tip.label
snv_dt <- snv_dt[cell_id %in% tips]
tips <- intersect(tips, unique(snv_dt$cell_id))
tree <- keep.tip(tree, tips)
setkey(snv_dt, cell_id)

tip_index <- setNames(seq_along(tips), tips)

# --- 2) Build per-locus ALT cell sets
snv_dt[, locus := paste(chrom, pos, sep=":")]
snv_dt[, alt_flag := alt_count >= alt_min]
A_list <- snv_dt[alt_flag == TRUE, .(cells = list(unname(tip_index[cell_id]))), by = locus]$cells
# Optional: filter loci with too few/many ALT cells if desired

# --- 3) BFS subtree selection
# Ensure tree is rooted
if (!is.rooted(tree)) tree <- phangorn::midpoint(tree)

N <- length(tree$tip.label)
# Fast leaf counts per node
desc_tips <- Descendants(tree, 1:Nnode(tree)+N, type="tips")
sizes <- vapply(desc_tips, length, integer(1))
node_order <- seq_along(desc_tips)  # BFS from root is fine since Descendants is order by node

selected <- logical(length(desc_tips))
subtree_masks <- list()

covered_tips <- rep(FALSE, N)  # to enforce non-nesting
for (i in node_order) {
  n_i <- sizes[i]
  if (n_i >= alpha*N && n_i <= beta*N) {
    mask <- rep(FALSE, N)
    mask[desc_tips[[i]]] <- TRUE
    # check non-nesting: skip if fully contained in an already selected subtree
    if (!all(mask & covered_tips == mask)) {
      selected[i] <- TRUE
      subtree_masks[[length(subtree_masks)+1]] <- mask
      # mark coverage to prevent choosing descendants
      covered_tips <- covered_tips | mask
    }
  }
}

# --- 4) Observed support counts
supports_obs <- integer(length(subtree_masks))
not_masks <- lapply(subtree_masks, function(m) which(!m))

for (s in seq_along(subtree_masks)) {
  nm <- not_masks[[s]]
  # An SNV supports S if none of its ALT cells fall outside S
  cnt <- 0L
  for (A in A_list) {
    if (!any(A %in% nm)) cnt <- cnt + 1L
  }
  supports_obs[s] <- cnt
}

# --- 5) Permutation test
set.seed(1)
perm_mat <- replicate(B, sample.int(N), simplify = FALSE)  # list of length B of integer permutations

supports_perm_ge <- integer(length(subtree_masks))  # counts of perm >= obs

# Precompute A indices as a matrix (rows = A sets, cols = max length)
max_A_len <- max(lengths(A_list))
A_indices_matrix <- matrix(NA_integer_, nrow = length(A_list), ncol = max_A_len)
A_lengths <- lengths(A_list)
for (i in seq_along(A_list)) {
  A_indices_matrix[i, 1:A_lengths[i]] <- A_list[[i]]
}

# Main optimized loop
for (b in seq_len(B)) {
  print(b)
  if (b %% 100 == 0) print(b)
  perm <- perm_mat[[b]]
  
  for (s in seq_along(subtree_masks)) {
    nm_perm_set <- perm[not_masks[[s]]]
    nm_perm_lookup <- logical(max(perm))  # Assuming perm values are reasonable
    nm_perm_lookup[nm_perm_set] <- TRUE
    
    # Vectorized check across all A sets
    cnt <- 0L
    for (i in seq_along(A_list)) {
      A_len <- A_lengths[i]
      A_indices <- A_indices_matrix[i, 1:A_len]
      perm_A_indices <- perm[A_indices]
      
      # Check if any permuted A index is in nm_perm_set
      if (!any(nm_perm_lookup[perm_A_indices])) {
        cnt <- cnt + 1L
      }
    }
    
    if (cnt >= supports_obs[s]) supports_perm_ge[s] <- supports_perm_ge[s] + 1L
  }
}


pvals <- (1 + supports_perm_ge) / (1 + B)
pvals_adj <- p.adjust(pvals, method = "BH")

results <- data.table(
  subtree_id = seq_along(subtree_masks),
  size = vapply(subtree_masks, sum, integer(1)),
  support_obs = supports_obs,
  pval = pvals,
  qval = pvals_adj
)

# 'results' now tells you which subtrees are significantly supported by SNVs

