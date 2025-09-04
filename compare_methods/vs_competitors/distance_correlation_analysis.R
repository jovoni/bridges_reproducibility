.libPaths("/home/santacg/R_env/lib/R/library")
rm(list = ls())
require(tidyverse)
library(ape)
library(ggplot2)
library(corrplot)
#library(GGally)

params_df = read.delim("../data/param_grid.csv", sep = ",")
DATA_DIR = "../data/"
RES_DIR = "results/"
dir.create(RES_DIR, recursive = T)

read_medicc_D = function(sim_id) {
  df <- read.delim(file.path("results",sim_id,"medicc2","medicc_input_pairwise_distances.tsv"), 
                   sep = "\t", row.names = 1, check.names = FALSE)
  mat = as.matrix(df)
  if ("diploid" %in% colnames(mat)) {
    mat = mat[rownames(mat) != "diploid", colnames(mat) != "diploid"]
  }
}

# Function to extract upper triangular values (excluding diagonal)
extract_upper_tri <- function(matrix) {
  upper_indices <- upper.tri(matrix, diag = FALSE)
  return(matrix[upper_indices])
}

read_phy = function(path) {
  lines <- readLines(path)
  
  # Parse the number of taxa
  n_taxa <- as.integer(lines[1])
  
  # Read the distance matrix into a data frame
  df <- read.table(text = lines[-1], stringsAsFactors = FALSE)
  
  # Extract labels and matrix
  labels <- df[[1]]
  mat <- as.matrix(df[, -1])
  mode(mat) <- "numeric"
  rownames(mat) <- labels
  colnames(mat) <- labels
  
  # Convert to dist object
  #dist_mat <- as.dist(mat)
  mat
}

get_correlation_results = function(sim_idx) {
  sim = params_df[sim_idx, ]
  sim_id = sim$sim_id
  
  # Bridges
  bridges_D = readRDS(file.path(RES_DIR, sim_id, "bridges_D.RDS"))
  cell_order = colnames(bridges_D)
  
  # # Medalt
  # medalt_D = readRDS(file.path(RES_DIR, sim_id, "medalt_D.RDS"))
  # medalt_D = medalt_D[cell_order, cell_order]
  # medalt_D[medalt_D==1e+06] = 0
  
  # Root
  root_D = read_phy(file.path(RES_DIR, sim_id, "dice/standard_root_balME_dm.phy"))
  #root_D = readRDS(file.path(RES_DIR, sim_id, "dice/standard_root_balME_dm.phy"))
  root_D = root_D[cell_order, cell_order]
  
  # Medicc2
  medicc_D = read_medicc_D(sim_id)
  medicc_D = medicc_D[cell_order, cell_order]
  
  # Extract upper triangular values for each matrix
  bridges_upper <- extract_upper_tri(bridges_D)
  #medalt_upper <- extract_upper_tri(medalt_D)
  root_upper <- extract_upper_tri(root_D)
  medicc_upper <- extract_upper_tri(medicc_D)
  
  # Create a data frame with all upper triangular values
  correlation_data <- tibble(
    Bridges = bridges_upper,
    #Medalt = medalt_upper,
    Root = root_upper,
    Medicc2 = medicc_upper
  )
  
  # Calculate correlation matrix
  cor_matrix <- cor(correlation_data, use = "complete.obs")
  
  # Create correlation results tibble
  correlation_results <- tibble(
    Method1 = character(),
    Method2 = character(),
    Correlation = numeric(),
    P_value = numeric(),
    RMSE = numeric()
  )
  
  # Calculate pairwise correlations with p-values and RMSE
  methods <- colnames(correlation_data)
  for (i in 1:(length(methods)-1)) {
    for (j in (i+1):length(methods)) {
      cor_test <- cor.test(correlation_data[[methods[i]]], 
                           correlation_data[[methods[j]]], 
                           use = "complete.obs")
      
      # Calculate RMSE
      valid_indices <- complete.cases(correlation_data[[methods[i]]], correlation_data[[methods[j]]])
      rmse_value <- sqrt(mean((correlation_data[[methods[i]]][valid_indices] - 
                                 correlation_data[[methods[j]]][valid_indices])^2))
      
      correlation_results <- correlation_results %>%
        add_row(
          Method1 = methods[i],
          Method2 = methods[j],
          Correlation = cor_test$estimate,
          P_value = cor_test$p.value,
          RMSE = rmse_value
        )
    }
  }
  
  # Add significance stars and round values
  correlation_results <- correlation_results %>%
    mutate(
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        P_value < 0.1 ~ ".",
        TRUE ~ ""
      ),
      Correlation_rounded = round(Correlation, 3),
      P_value_rounded = round(P_value, 4),
      RMSE_rounded = round(RMSE, 4)
    )
  
  correlation_results
}

get_rmse_and_cor_plots = function(correlation_results) {
  # Create a separate RMSE visualization
  rmse_plot <- correlation_results %>%
    ggplot(aes(x = paste(Method1, "vs", Method2), y = RMSE)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = round(RMSE, 3)), vjust = -0.3, size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = "Root Mean Square Error (RMSE) Between Distance Methods",
      x = "Method Comparison",
      y = "RMSE",
      caption = "Lower RMSE indicates better agreement between methods"
    )
  
  
  # Create a separate Cor visualization
  cor_plot <- correlation_results %>%
    ggplot(aes(x = paste(Method1, "vs", Method2), y = Correlation)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = round(Correlation, 3)), vjust = -0.3, size = 3) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = "Correlation Between Distance Methods",
      x = "Method Comparison",
      y = "Spearman Correlation"
    ) +
    ylim(c(NA, 1))
  
  list(rmse_plot=rmse_plot, cor_plot=cor_plot)
}

sim_idx = 50
params_df[sim_idx,]
plots = get_rmse_and_cor_plots(get_correlation_results(sim_idx))
plots$rmse_plot
plots$cor_plot

#N = nrow(params_df)
N = 60
df = lapply(1:N, function(sim_idx) {
  print(sim_idx)
  get_correlation_results(sim_idx) %>% 
    dplyr::bind_cols(params_df[sim_idx,])
}) %>% do.call("bind_rows", .)
df$BFB_prop = df$bfb_rate / (df$bfb_rate + df$amp_rate + df$del_rate + df$normal_rate)

cols = c(
  "Bridges vs Medalt" = "gray20",
  "Bridges vs Root"= "gray20",
  "Bridges vs Medicc2"= "indianred",
  "Medalt vs Root"= "gray20",
  "Medalt vs Medicc2" = "gray20",
  "Root vs Medicc2"= "gray20"
)

df %>%
  dplyr::filter(ncells > 100) %>% 
  ggplot(aes(x = paste(Method1, "vs", Method2), y = RMSE, col = paste(Method1, "vs", Method2))) +
  geom_violin() +
  theme_bw() +
  scale_color_manual(values = cols) +
  facet_wrap(~paste0(ncells, " cells")) +
  scale_y_continuous(transform = "log10") +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Root Mean Square Error (RMSE) Between Distance Methods",
    x = "Method Comparison",
    y = "RMSE",
    caption = "Lower RMSE indicates better agreement between methods"
  )

df %>%
  dplyr::filter(ncells > 100) %>% 
  ggplot(aes(x = paste(Method1, "vs", Method2), y = RMSE, col = paste(Method1, "vs", Method2))) +
  geom_violin() +
  theme_bw() +
  scale_color_manual(values = cols) +
  facet_wrap(~paste0("BFB rate = ",BFB_prop)) +
  scale_y_continuous(transform = "log10") +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Root Mean Square Error (RMSE) Between Distance Methods",
    x = "Method Comparison",
    y = "RMSE",
    caption = "Lower RMSE indicates better agreement between methods"
  )


df %>%
  dplyr::filter(ncells > 100) %>% 
  ggplot(aes(x = paste(Method1, "vs", Method2), y = Correlation, col = paste(Method1, "vs", Method2))) +
  geom_violin() +
  theme_bw() +
  scale_color_manual(values = cols) +
  facet_wrap(~paste0(ncells, " cells")) +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Correlation Between Distance Methods",
    x = "Method Comparison",
    y = "Spearman Correlation"
  )

df %>%
  ggplot(aes(x = paste(Method1, "vs", Method2), y = Correlation, col = paste(Method1, "vs", Method2))) +
  geom_violin() +
  theme_bw() +
  facet_wrap(~paste0("BFB rate = ", BFB_prop)) +
  coord_flip() +
  scale_color_manual(values = cols) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    title = "Correlation Between Distance Methods",
    x = "Method Comparison",
    y = "Spearman Correlation"
  )

