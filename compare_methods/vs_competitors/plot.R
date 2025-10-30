
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
source("utils_plot.R")
library(ggsci)

MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)

NAME_MAPPING = list(
  "bridges" = "BRIDGES",
  "dice" = "DICE",
  "medicc" = "MEDICC2",
  "hamming" = "Hamming",
  "lazac" = "Lazac",
  "euclidean" = "Euclidean",
  "sitka" = "Sitka"
)

dir.create("plot")
a = alpha_quantile = .25

res = readRDS("results/results_summary.RDS")
res = res %>%
  dplyr::group_by(sim_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n == max(n)) %>%
  dplyr::mutate(algorithm = ifelse(algorithm == "root", "dice", algorithm))
res$algorithm = unlist(NAME_MAPPING[res$algorithm])
res$algorithm = factor(res$algorithm, levels = c("BRIDGES", "MEDICC2", "DICE", "Hamming", "Lazac", "Euclidean", "Sitka"))

# res %>%
#   dplyr::filter(metric == "seconds") %>%
#   dplyr::group_by(algorithm, ncells) %>%
#   dplyr::summarise(m = mean(value)) %>%
#   dplyr::group_by(ncells) %>%
#   dplyr::mutate(f = m / m[algorithm == "bridges"]) %>%
#   dplyr::select(algorithm, ncells, f) %>%
#   tidyr::pivot_wider(values_from = f, names_from = algorithm) %>% view()

plots_time = time_plots(res)
plots_RFnorm = metric_plot(res, "RF normalized")
plots_RF = metric_plot(res, "RF distance")
plots_GRF = metric_plot(res, "Generalized RF")
plots_MCI = metric_plot(res, "Mutual Clustering Info")
plots_QD = metric_plot(res %>% dplyr::filter(ncells < 1000), "Quartet divergence")

final_plot = patchwork::free(plots_RFnorm$p_absolute) +
  patchwork::free(plots_time$p_time + theme(legend.position = "bottom")) +
  plot_layout(design = "AAAABB") +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    plot.tag = element_text(face = 'bold', size = 15),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

# export at 589 pt width (~8.15 in) and 11 pt text
ggsave(
  "../../figures/sim_performances.pdf",
  plot = final_plot,
  width = 8.15, height = 4.4, units = "in",  # adjust height to fit layout
  dpi = 72,  # ensures 1 pt ≈ 1 LaTeX pt
  useDingbats = FALSE
)

ggsave(
  "../../figures/supp_sim_speedup.pdf",
  plot = plots_time$p_time_frac,
  width = 8.15, height = 6.0, units = "in",
  dpi = 72, useDingbats = FALSE
)

# ggsave("../../figures/sim_performances.pdf", width = 11.7, height = 6.3, units = "in", plot = final_plot, dpi = 600)
# ggsave("../../figures/supp_sim_speedup.pdf", width = 8, height = 5.5, units = "in", plot = plots_time$p_time_frac)


library(ggplot2)
library(patchwork)

# 589 pt -> inches (LaTeX pt ~ 72.27 pt/in)
target_in <- 589 / 72.27  # ~8.15 in

# Patchwork design: 3 rows of A, 2 rows of B
design <- "
A
A
A
B
B
"

# Common theme (11 pt text at final size)
common_theme <- theme(
  text = element_text(size = 11),
  legend.title = element_text(face = "bold"),
  plot.tag = element_text(face = "bold", size = 15),
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  panel.grid.minor = element_blank()
)

# Helper to save with true point-for-point sizing
save_pdf <- function(file, plot, width_in = target_in, height_in = target_in){
  ggsave(
    filename = file, plot = plot,
    width = width_in, height = height_in, units = "in",
    dpi = 72, device = cairo_pdf, limitsize = FALSE
  )
}

# --- RF ---
p <- free(plots_RF$p_absolute) + free(plots_RF$p_norm + coord_flip()) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  common_theme

# Choose a sensible height (vertical stack: 3:2 rows) – tweak if needed
save_pdf("../../figures/supp_sim_RF.pdf", p, width_in = target_in, height_in = 9.2)

# --- QD ---
p <- free(plots_QD$p_absolute) + free(plots_QD$p_norm + coord_flip()) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  common_theme

save_pdf("../../figures/supp_sim_QD.pdf", p, width_in = target_in, height_in = 9.2)

# --- MCI ---
p <- free(plots_MCI$p_absolute) + free(plots_MCI$p_norm + coord_flip()) +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A") &
  common_theme

save_pdf("../../figures/supp_sim_MCI.pdf", p, width_in = target_in, height_in = 9.2)
