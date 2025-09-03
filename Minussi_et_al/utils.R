
SAMPLE_NAMES = c("TN1", "TN2")

run_umap <- function(popseg_long,
                     ploidy_VAL = NULL,
                     umap_dist = "manhattan",
                     umap_min_dist = 0,
                     umap_spread = 1,
                     umap_n_neighbors = 40,
                     mc.cores = 10,
                     seed = 55,
                     round = FALSE) {

  if (round == TRUE) {
    popseg_long <- ploidy_scale(ploidy_VAL = ploidy_VAL, popseg_long, round = round)
  }

  message("Constructing UMAP embedding.")
  set.seed(seed)
  dat_umap <- uwot::umap(popseg_long,
                         metric = umap_dist,
                         min_dist = umap_min_dist,
                         n_neighbors = umap_n_neighbors,
                         spread = umap_spread,
                         n_components = 2,
                         n_thread = mc.cores, )

  umap_df <- as.data.frame(dat_umap)

  rownames(umap_df) <- rownames(popseg_long)
  umap_df$cell <- rownames(umap_df)

  return(umap_df)
}

run_clustering <- function(umap_df,
                           k_snn_major = 35,
                           k_snn_minor = 17) {

  library(scran)
  # building a snn graph for superclones
  message("Building SNN graph.")
  g_major <- scran::buildSNNGraph(umap_df[,c(1:2)], k = k_snn_major, transposed = T)

  g_clusters <- igraph::membership(igraph::components(g_major))
  g_clusters <- paste0("s", g_clusters)

  # Clustering
  message("Running hdbscan.")
  subclones <- dbscan::hdbscan(umap_df[,c(1:2)],
                               minPts = k_snn_minor)
  umap_df$subclones <- paste0("c",subclones$cluster)

  # for hdb
  # adding the ones classified as outliers to the closest cluster possible according to euclidean distance
  dist_umap <- dist(umap_df[,c(1:2)]) %>% as.matrix() %>% as.data.frame() %>%
    rownames_to_column("cell2") %>%
    gather(key = "cell1",
           value = "dist",
           -cell2) %>%
    dplyr::filter(cell1 != cell2)

  dist_min <- dist_umap %>%
    right_join(umap_df %>% dplyr::select(cell, subclones), by = c("cell2" = "cell")) %>%
    filter(subclones != "c0") %>%
    group_by(cell1) %>%
    slice_min(dist) %>%
    ungroup()


  for (i in 1:nrow(umap_df)) {

    if(umap_df$subclones[i] == "c0") {
      cellname <- rownames(umap_df)[i]
      closest_cell <- filter(dist_min, cell1 == rownames(umap_df)[i])$cell2
      closest_cell_cluster <- filter(umap_df, cell == closest_cell)$subclones
      umap_df$subclones[i] <- closest_cell_cluster

    }

  }

  cl_df <- tibble::tibble(
    superclones = g_clusters,
    subclones =  umap_df$subclones,
    cells = rownames(umap_df)
  )

  cl_df <- cl_df %>%
    arrange(superclones, subclones)

  # calculating number of cells in every cluster
  freq_df <- janitor::tabyl(umap_df$subclones)
  names(freq_df)[1] <- "cluster"
  print(freq_df)

  classification <- cl_df

  message("Done.")

  return(classification)

}
