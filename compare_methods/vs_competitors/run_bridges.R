
rm(list = ls())
require(tidyverse)
library(ape)
library(ggplot2)
library(tidyverse)
library(bridges)
source("../utils_comparisons.R")
source("utils_distances.R")

# Get the task ID from SLURM environment variable
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

params_df = read.delim("../data/param_grid.csv", sep = ",")
DATA_DIR = "../data/"
RES_DIR = "results/"
dir.create(RES_DIR, recursive = T)


# Process only the row corresponding to this task ID
i <- task_id
print(paste("Processing task", i, "of", nrow(params_df)))

sim_id = params_df[i,]$sim_id

# Read input matrix
input = readRDS(file.path(DATA_DIR, sim_id, "simulation.RDS"))

s = Sys.time()
res = bridges::fit(input$cna_data, alleles = c("A", "B"), 
                   k_jitter_fix = 0, bfb_penalty = 0)
e = Sys.time()
bridges_time = as.numeric(e - s, unit = "secs")

bridges_tree = res$tree
bridges_D = res$D

time_tibble = dplyr::tibble(
  method = c("bridges"),
  seconds = c(bridges_time)
)

dir.create(file.path(RES_DIR, sim_id), recursive = T)
saveRDS(time_tibble, file.path(RES_DIR, sim_id, "time_tibble.RDS"))
saveRDS(bridges_tree, file.path(RES_DIR, sim_id, "bridges_tree.RDS"))
saveRDS(bridges_D, file.path(RES_DIR, sim_id, "bridges_D.RDS"))

print(paste("Completed task", i, "for sim_id:", sim_id))