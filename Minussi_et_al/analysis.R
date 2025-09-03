
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

  res = readRDS(paste0("results/", "briges_fit_", sample_id, ".rds"))
  classif = readRDS(paste0("results/", "classification_", sample_id, ".rds"))
  umap_df = readRDS(paste0("results/", "umap_df_", sample_id, ".rds"))

  hm = bridges::plot_heatmap(data, res$tree, to_plot = c("A", "B"), use_raster = F, annotations = classif, ladderize = T)

  umap = umap_df %>%
    dplyr::mutate(cell_id = cell) %>%
    dplyr::left_join(classif, by = "cell_id") %>%
    ggplot(mapping = aes(x = V1, y = V2, col = subclones, shape = superclones)) +
    geom_point() +
    theme_bw() +
    labs(x = "UMAP 1", y = "UMAP 2")

  ggsave(paste0("plot/", sample_id, "_umap.pdf"), width = 6, height = 6, plot = umap)
  pdf(paste0("plot/", sample_id, "_gm.pdf"), width = 16, height = 10)
  print(hm)
  dev.off()
}




