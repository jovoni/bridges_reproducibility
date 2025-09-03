
rm(list = ls())
library(tidyverse)
library(bridges)
source("utils.R")

# Fit data
for (sample_id in SAMPLE_NAMES) {
  cna_data = read.delim(paste0("data/", sample_id, ".tsv"), sep = "\t")

  data = cna_data %>%
    dplyr::mutate(cell_id = sample_id, chr = str_replace(chrom, "chr", ""), A = cn_a, B = cn_b) %>%
    dplyr::select(cell_id, chr, A, B, start, end)

  res = bridges::fit(data = data, alleles = c("A", "B"))

  popseg_long = lapply(unique(cna_data$sample_id), function(cid) {
    df = cna_data %>% dplyr::filter(sample_id == cid) %>%
      dplyr::arrange(chrom, start)

    c(df$cn_a, df$cn_b)
    #df$cn_a + df$cn_b
  }) %>% do.call(rbind, .)
  rownames(popseg_long) = unique(cna_data$sample_id)

  umap_df = run_umap(popseg_long)

  if (sample_id == "TN2") {
    classif = run_clustering(umap_df, k_snn_major = 63, k_snn_minor = 17)
  } else {
    classif = run_clustering(umap_df, k_snn_major = 45, k_snn_minor = 17)
  }
  classif$cell_id = classif$cells
  classif$cells = NULL

  dir.create("results", recursive = T)

  saveRDS(res, file = paste0("results/", "briges_fit_", sample_id, ".rds"))
  saveRDS(classif, file = paste0("results/", "classification_", sample_id, ".rds"))
  saveRDS(umap_df, file = paste0("results/", "umap_df_", sample_id, ".rds"))
}


#bridges::plot_heatmap(data, res$tree, to_plot = c("A", "B"), use_raster = F, annotations = classif, ladderize = T)




