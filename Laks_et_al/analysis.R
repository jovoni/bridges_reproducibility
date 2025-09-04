
rm(list = ls())
library(tidyverse)

clone_colors = c(
  "A" = '#a6cee3',
  "B" = '#1f78b4',
  "C" = '#b2df8a',
  "D" = '#33a02c',
  "E" = '#fb9a99',
  "F" = '#e31a1c',
  "G" = '#fdbf6f',
  "H" = '#ff7f00',
  "I" = '#cab2d6'
)

bridges_fit = readRDS("results/bridges_fit.rds")
#bridges_fit = readRDS("results/bridges_fit_wo_reconstruction.rds")
clone_df = read.delim("data/ov2295_clone_clusters.csv", sep = ",")
clone_snvs_df <- read_csv("data/ov2295_clone_snvs.csv")
clone_snvs_df = clone_snvs_df %>%
  dplyr::group_by(clone_id) %>%
  dplyr::summarise(s = sum(is_present))

library(ape)
tree = bridges_fit$tree
leaves_lenghts_df <- data.frame(
  cell_id = tree$tip.label,
  length = node.depth.edgelength(tree)[1:Ntip(tree)]
) %>% dplyr::left_join(clone_df, by = "cell_id") %>%
  dplyr::left_join(clone_snvs_df, by = "clone_id")

p_SNV_CNA_corr = leaves_lenghts_df %>%
  dplyr::group_by(clone_id) %>%
  dplyr::mutate(avg_length = mean(length), sd_length = sd(length)) %>%
  dplyr::select(clone_id, s, avg_length, sd_length) %>%
  dplyr::distinct() %>%
  ggplot(mapping = aes(x = s, y = avg_length)) +
  geom_smooth(method = "lm", col = "black") +
  ggpubr::stat_cor(p.accuracy = 0.001) +
  geom_point(mapping = aes(col = clone_id), size = 3) +
  theme_bw() +
  scale_color_manual(values = clone_colors) +
  labs(x = "SNVs root-to-clade length",
       y = "Average bridges SCNA root-to-leaf length",
       col = "Laks et al. clone id")

p_SNV_CNA_corr
ggsave("plot/SNV_CNA_correlation.pdf", width = 9, height = 7, plot = p_SNV_CNA_corr)

cn_data = read.delim("data/ov2295_cell_cn.csv", sep = ",")
clone_df = read.delim("data/ov2295_clone_clusters.csv", sep = ",")

N = 900
N = min(N, nrow(clone_df))
cell_ids = sample(clone_df$cell_id, size = N)

cn_data = cn_data %>%
  dplyr::filter(cell_id %in% cell_ids) %>%
  dplyr::mutate(CN = state) %>%
  dplyr::select(cell_id, sample_id, library_id, start, end, chr, CN)

hm = bridges::plot_heatmap(cn_data, tree = bridges_fit$tree, to_plot = "CN", use_raster = F, ladderize = T, annotations = clone_df)

png("plot/heatmap.png", width = 16, height = 10, units = "in", res = 300)
print(hm)
dev.off()


#bridges_fit = readRDS("results/bridges_fit.rds")
bridges_fit = readRDS("results/bridges_fit_wo_reconstruction.rds")
clone_df = read.delim("data/ov2295_clone_clusters.csv", sep = ",")
clone_snvs_df <- read_csv("data/ov2295_clone_snvs.csv")
clone_snvs_df = clone_snvs_df %>%
  dplyr::group_by(clone_id) %>%
  dplyr::summarise(s = sum(is_present))

library(ape)
tree = bridges_fit$tree
leaves_lenghts_df <- data.frame(
  cell_id = tree$tip.label,
  length = node.depth.edgelength(tree)[1:Ntip(tree)]
) %>% dplyr::left_join(clone_df, by = "cell_id") %>%
  dplyr::left_join(clone_snvs_df, by = "clone_id")

leaves_lenghts_df %>%
  dplyr::group_by(clone_id) %>%
  dplyr::mutate(avg_length = mean(length), sd_length = sd(length)) %>%
  dplyr::select(clone_id, s, avg_length, sd_length) %>%
  dplyr::distinct() %>%
  ggplot(mapping = aes(x = s, y = avg_length)) +
  geom_smooth(method = "lm", col = "black") +
  ggpubr::stat_cor(p.accuracy = 0.001) +
  geom_point(mapping = aes(col = clone_id), size = 3) +
  theme_bw() +
  scale_color_manual(values = clone_colors) +
  labs(x = "SNVs root-to-clade length",
       y = "Average NJ root-to-leaf length",
       col = "Laks et al. clone id")

ggsave("plot/SNV_NJ_correlation.pdf", width = 9, height = 7, plot = last_plot())
