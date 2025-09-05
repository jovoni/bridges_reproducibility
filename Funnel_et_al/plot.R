
library(ggsci)

df_clonal_discordance = readRDS("results/metrics/clonal_discordance.rds")
dissmilarities = readRDS("results/metrics/dissimilarities.rds")
youden_df = readRDS("results/metrics/youden.rds")

p_clonal_discordance = df_clonal_discordance %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::mutate(Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot(mapping = aes(x = sample_id, y = value, fill = method)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_fill_bmj() +
  labs(x = "Sample", y = "Clonal discordance", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p_dissmilarities = dissmilarities %>%
  na.omit() %>%
  dplyr::group_by(name, sample_id) %>%
  dplyr::summarise(m = mean(diss)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::mutate(Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot2::ggplot(mapping = aes(x = sample_id, y = m, fill = name)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_bmj() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  labs(x = "Sample", y = "Mean Sibling Dissimilarity", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p_jouden = youden_df %>%
  na.omit() %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(Dataset = ifelse(n == 3, "Funnel", "Funnel+")) %>%
  dplyr::mutate(Dataset = factor(Dataset, levels = c("Funnel", "Funnel+"))) %>%
  ggplot2::ggplot(mapping = aes(x = sample_id, y = Jouden, fill = mehtod)) +
  geom_col(position = "dodge") +
  theme_bw() +
  facet_grid(~Dataset, scales = "free_x", space = "free_x") +
  scale_fill_bmj() +
  labs(x = "Sample", y = "Jouden Index", fill = "Algorithm") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

p = p_clonal_discordance + p_dissmilarities + p_jouden +
  plot_layout(design = "AABB\n#CC#", guides = "collect") +
  plot_annotation(tag_levels = c("A")) & theme(legend.position = "bottom")

ggsave("../figures/funnel_cohort_v_sitka.pdf", width = 12, height = 6, units = "in")
