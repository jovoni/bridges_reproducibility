
rm(list = ls())
library(tidyverse)

cn_data = read.delim("data/ov2295_cell_cn.csv", sep = ",")
clone_df = read.delim("data/ov2295_clone_clusters.csv", sep = ",")
cell_ids = unique(clone_df$cell_id)

cn_data = cn_data %>%
  dplyr::filter(cell_id %in% cell_ids) %>%
  dplyr::mutate(CN = state) %>%
  dplyr::select(cell_id, sample_id, library_id, start, end, chr, CN)

CN_matrix <- cn_data %>%
  arrange(cell_id, chr, start) %>%  # Sort once
  group_by(cell_id) %>%
  mutate(position = row_number()) %>%  # Create position index
  ungroup() %>%
  select(cell_id, position, CN) %>%
  pivot_wider(names_from = cell_id, values_from = CN) %>%
  select(-position) %>%
  as.matrix() %>%
  t()

dim(CN_matrix)

cn_clustering <- function(
    M,
    method = c("gmm", "hdbscan"),
    gmm_k = 20,
    hdbscan_minPts = 10,
    min_cells = 50,
    seed = 123,
    umap_args = list()
) {
  method <- match.arg(method)
  stopifnot(is.matrix(M) || is.data.frame(M))
  M <- as.matrix(M)
  if (!is.numeric(M)) stop("M must be numeric.")
  
  # Dependencies
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Please install 'uwot'.")
  if (method == "gmm" && !requireNamespace("mclust", quietly = TRUE))
    stop("Please install 'mclust' for GMM.")
  if (method == "hdbscan" && !requireNamespace("dbscan", quietly = TRUE))
    stop("Please install 'dbscan' for HDBSCAN.")
  
  # UMAP
  set.seed(seed)
  # umap_defaults <- list(n_components = 2, n_neighbors = 15, min_dist = 0.1,
  #                       metric = "euclidean", verbose = FALSE, ret_model = FALSE)
  # umap_params <- utils::modifyList(umap_defaults, umap_args, keep.null = TRUE)
  # emb <- do.call(uwot::umap, c(list(X = M), umap_params))
  # embedding <- data.frame(UMAP_1 = emb[,1], UMAP_2 = emb[,2], row.names = rownames(M))
  
  library(umap)
  emb = umap::umap(M)
  embedding <- data.frame(UMAP_1 = emb$layout[,1], UMAP_2 = emb$layout[,2], row.names = rownames(M))
  
  # Clustering
  if (method == "gmm") {
    # Over-specify K at gmm_k; let mclust pick covariance type via BIC (G fixed).
    library(mclust)
    fit <- mclust::Mclust(data = embedding, G = gmm_k)
    cl_raw <- fit$classification
    # mclust labels clusters 1..K; no noise class
    noise_label <- NULL
  } else {
    # HDBSCAN; cluster 0 is noise/outliers -> drop
    fit <- dbscan::hdbscan(embedding, minPts = hdbscan_minPts)
    cl_raw <- fit$cluster
    noise_label <- 0L
  }
  
  # Remove noise/outlier cluster (HDBSCAN) and small clusters
  cl <- cl_raw
  if (!is.null(noise_label)) cl[cl == noise_label] <- NA_integer_
  
  # Recode to consecutive integers for kept clusters after size filtering
  sizes <- table(cl, useNA = "no")
  keep_ids <- names(sizes)[sizes >= min_cells]
  cl[!(cl %in% as.integer(keep_ids))] <- NA_integer_
  
  # Compact the labels (1..K_kept) but preserve mapping
  kept_levels <- sort(unique(stats::na.omit(cl)))
  recode_map <- setNames(seq_along(kept_levels), kept_levels)
  cl_recode <- ifelse(is.na(cl), NA_integer_, recode_map[as.character(cl)])
  cluster_factor <- factor(cl_recode, levels = seq_along(kept_levels))
  
  kept_clusters_df <- data.frame(
    cluster_id = seq_along(kept_levels),
    size = as.integer(table(cluster_factor)),
    row.names = NULL
  )
  
  # Median CN per cluster
  if (length(kept_levels) > 0) {
    by_cluster <- split(data.frame(M, check.names = FALSE), cluster_factor)
    # remove NA group if present
    by_cluster <- by_cluster[!is.na(names(by_cluster))]
    med_list <- lapply(by_cluster, function(df) {
      apply(as.matrix(df), 2, stats::median, na.rm = TRUE)
    })
    cluster_medians <- do.call(rbind, med_list)
    rownames(cluster_medians) <- paste0("cluster_", seq_len(nrow(cluster_medians)))
  } else {
    cluster_medians <- matrix(numeric(0), nrow = 0, ncol = ncol(M),
                              dimnames = list(NULL, colnames(M)))
  }
  
  list(
    embedding = embedding,
    cluster = cluster_factor,          # NA for noise/small clusters
    kept_clusters = kept_clusters_df,  # cluster_id, size
    cluster_medians = cluster_medians
  )
}

res_gmm <- cn_clustering(M = CN_matrix, method = "gmm", gmm_k = 20, min_cells = 50, seed = 1234)
res_gmm$embedding$cluster = res_gmm$cluster
res_gmm$embedding$cell_id = rownames(res_gmm$embedding)

saveRDS(res_gmm$embedding %>% dplyr::left_join(clone_df), "results/embedding.RDS")

res_gmm$embedding %>% 
  dplyr::left_join(clone_df) %>% 
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2, col = clone_id)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = clone_colors)

