
library(ggsci)
library(tidyverse)
library(patchwork)

# --- target width: 589 pt -> inches (LaTeX pt ≈ 72.27 pt/in) ---
target_in <- 589 / 72.27  # ≈ 8.15 in
# Keep similar aspect ratio to your 15x9 in original -> height ≈ 9/15 * 8.15
target_h  <- target_in * 1.2

# --- load data ---
df_clonal_discordance <- readRDS("results/metrics/clonal_discordance.rds")
dissmilarities        <- readRDS("results/metrics/dissimilarities.rds")
youden_df             <- readRDS("results/metrics/youden.rds")
rf_df                 <- readRDS("results/metrics/RF_against.rds")

# Optional: fix common typos in columns
if ("Jouden" %in% names(youden_df)) youden_df <- youden_df %>% dplyr::rename(Youden = Jouden)
if ("mehtod" %in% names(youden_df)) youden_df <- youden_df %>% dplyr::rename(method = mehtod)

# --- common theme: 11 pt text at final size ---
common_theme <- theme(
  text = element_text(size = 11),
  legend.title = element_text(face = "bold"),
  plot.tag = element_text(face = "bold", size = 15),
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  panel.grid.minor = element_blank()
)

# --- patchwork design (two rows): "AABB" on row 1; "#CC#" on row 2 ---
design <- "
AABB
#CC#
"

design <- "
A
B
C
"

# (Optional) RF Normalized boxplot (not saved below, but keep if needed)
rf_df %>%
  group_by(sample_id) %>%
  mutate(n = n(),
         Dataset = ifelse(n == 2, "Funnel", "Funnel+")) %>%
  filter(metric == "RF normalized") %>%
  ggplot(aes(x = against, y = value, fill = against)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_bmj()

# --- Clonal discordance ---
p_clonal_discordance <- df_clonal_discordance %>%
  group_by(sample_id) %>%
  mutate(n = n(),
         Dataset = ifelse(n == 3, "Funnel", "Funnel+"),
         Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot(aes(x = sample_id, y = value, fill = method)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_fill_bmj() +
  labs(x = "Sample", y = "Clonal discordance", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# --- Mean sibling dissimilarity ---
p_dissmilarities <- dissmilarities %>%
  na.omit() %>%
  group_by(name, sample_id) %>%
  summarise(m = mean(diss), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(n = n(),
         Dataset = ifelse(n == 3, "Funnel", "Funnel+"),
         Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot(aes(x = sample_id, y = m, fill = name)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_fill_bmj() +
  labs(x = "Sample", y = "Mean Sibling Dissimilarity", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# --- Youden index ---
p_jouden <- youden_df %>%
  na.omit() %>%
  group_by(sample_id) %>%
  mutate(n = n(),
         Dataset = ifelse(n == 3, "Funnel", "Funnel+"),
         Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot(aes(x = sample_id, y = Youden, fill = method)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_fill_bmj() +
  labs(x = "Sample", y = "Youden Index", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# --- Combine & export at 589 pt width with true 11 pt text ---
p <- p_clonal_discordance + p_dissmilarities + p_jouden +
  plot_layout(design = design, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") &
  common_theme

ggsave(
  "../figures/funnel_cohort_v_sitka.pdf",
  plot = p,
  width = target_in, height = target_h, units = "in",
  dpi = 72, device = cairo_pdf, limitsize = FALSE
)
