
source("getters.R")
library(ape)
library(ggtree)
library(dplyr)
library(tibble)

print(s)

# ---- Load tree & labels ----
s = get_sample_names()[12]
obj <- readRDS(file.path("results/bridges_trees", paste0(s, ".rds")))
tree <- obj$tree

labels_df <- get_labels(s) %>%
  distinct(cell_id, .keep_all = TRUE)

# Keep only overlapping tips
keep <- intersect(labels_df$cell_id, tree$tip.label)
if (length(keep) == 0L) stop("No overlapping cells between labels and tree tips.")

tree <- ape::keep.tip(tree, keep)
labels_df <- labels_df %>% filter(cell_id %in% keep)

# (Optional) drop edge lengths if they cause warnings
tree$edge.length <- NULL

# ---- Build heatmap data (rows must be tips; columns are annotations) ----
# Here we only use clone_id. Add more columns if you have them.
hm <- labels_df %>%
  select(cell_id, clone_id) %>%
  mutate(clone_id = as.factor(clone_id)) %>%
  column_to_rownames("cell_id")

# Order rows exactly like tip labels (CRUCIAL for gheatmap)
hm <- hm[tree$tip.label, , drop = FALSE]

# ---- Circular tree + heatmap ----
circ <- ggtree(tree, layout = "circular")
p1 <- gheatmap(
  circ, hm,
  offset = 0.5, width = 0.1,
  colnames = TRUE, colnames_angle = 1, colnames_offset_y = 0.25, legend_title = "Clone", font.size = 0
) + ggtitle("BRIDGES")
print(p1)
