
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(patchwork)
source("utils_plot.R")
library(ggsci)

dir.create("plot")
a = alpha_quantile = .25

res = readRDS("results/results_summary.RDS")
res = res %>%
  dplyr::group_by(sim_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n == max(n)) %>% 
  dplyr::mutate(algorithm = ifelse(algorithm == "root", "dice", algorithm))

res$algorithm = factor(res$algorithm, levels = c("bridges", "medicc", "dice", "hamming", "euclidean"))

p_time = res %>% 
  dplyr::filter(metric == "seconds") %>% 
  ggplot(mapping = aes(x = as.factor(ncells), y = value, col = algorithm)) +
  geom_jitter() +
  #geom_line() +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  scale_colour_nejm() +
  labs(x = "Number of cells", y = "Time (s)", col = "Algorithm")
p_time

p_RF_normed = plot_normed_boxplots(res, "RF distance", split_by = "BFB") +
  scale_colour_nejm() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  labs(y = "Norm'ed RF distance", col = "Algorithm")
p_RF_normed


final_plot = patchwork::free(p_RF_normed) + patchwork::free(p_time) +
  plot_layout(design = "AAABB") +
  plot_annotation(tag_levels = "A")
ggsave("../../figures/sim_performances.pdf", width = 8, height = 3.5, units = "in")

stop()

plot_trend(res, metric_name = "RF normalized", logy = F)
plot_trend(res, metric_name = "Quartet divergence", logy = F)

plot_normed_boxplots(res, "RF normalized", split_by = "BFB")
plot_normed_boxplots(res, "RF normalized")

#res = dplyr::bind_rows(res, readRDS("../method_details/results/tree_algorithms_comparisons/results_summary.RDS"))
res$bfb_prop = res$bfb_rate / (res$normal_rate + res$del_rate + res$amp_rate + res$bfb_rate)
res$bfb_prop %>% table()
res %>% dplyr::group_by(algorithm) %>% dplyr::summarise(n = n()) %>%
  dplyr::mutate(f = n / max(n))

res = res %>% dplyr::filter(algorithm %in% c("bridges", "medalt", "medicc", "nj", "root"))

summary_bar_plot = res %>%
  dplyr::filter(metric != "RF distance") %>%
  dplyr::group_by(algorithm, metric) %>%
  na.omit() %>%
  dplyr::summarise(m = median(value), qlow=quantile(value, alpha_quantile), qhigh=quantile(value, 1-alpha_quantile)) %>%
  na.omit() %>%
  ggplot(mapping = aes(x=algorithm, y=m, fill=algorithm, ymin=qlow, ymax=qhigh)) +
  geom_col(position = "dodge") +
  geom_errorbar() +
  ggh4x::facet_nested(~"Metric"+metric, scales = "free_y", independent = "y") +
  theme_bw() +
  labs(x = "", y="Value", col = "Method", fill="Method") +
  theme(legend.position = "bottom")
summary_bar_plot


general_performance_plot = res %>%
  dplyr::filter(metric != "RF distance", normal_rate != 1) %>%
  dplyr::group_by(algorithm, ncells, bfb_prop, metric, normal_rate) %>%
  dplyr::summarise(y = mean(value), yl = quantile(value, a, na.rm = TRUE), yh = quantile(value, 1-a, na.rm = TRUE)) %>%
  ggplot(mapping = aes(x = ncells, y=y, ymax=yh, ymin=yl, col=algorithm)) +
  geom_pointrange() +
  geom_line() +
  #scale_x_continuous(transform = "log10") +
  ggh4x::facet_nested("Metric"+metric~"BFB prop"+bfb_prop, scales = "free_y") +
  theme_bw() +
  labs(x = "Number of leaves", y="Value", col = "Method") +
  theme(legend.position = "bottom")
general_performance_plot

ggsave("plot/summary_bar_plot.png", summary_bar_plot, width = 8, height = 5, units = "in", dpi = 300)
ggsave("plot/general_performance_plot.png", general_performance_plot, width = 10, height = 7, units = "in", dpi = 300)


res %>%
  dplyr::filter(metric != "RF distance") %>%
  #dplyr::filter(algorithm != "medalt") %>%
  dplyr::mutate(BFB_level = ifelse(bfb_prop == 0, "No BFB", ifelse(bfb_prop >= .5 , "High", "Low"))) %>%
  dplyr::mutate(Ncells = ifelse(ncells <= 100, "<= 100", "> 100")) %>%
  dplyr::group_by(metric, sim_id) %>%
  dplyr::mutate(normed_value = value / max(value)) %>%
  ggplot(mapping = aes(x = algorithm, y=normed_value)) +
  geom_boxplot() +
  geom_jitter(alpha = .5, mapping = aes(col = algorithm)) +
  ggh4x::facet_nested(Ncells~metric+BFB_level) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_blank())

ggsave("plot/relative_boxplots.png", width = 10, height = 7, units = "in", dpi = 300)
