
library(ape)
library(phangorn)

get_sample_names = function() {
  stringr::str_replace_all(list.files("signatures_dataset/DLP/CNA/persample/"), "_hscn.csv.gz", "")
}

get_X_for_dissimilarity = function(sample_id, allele = "CN") {
  data = get_input(sample_id)
  split_by_chr <- split(data, data$chr)
  split_by_chr <- split_by_chr[names(split_by_chr)]
  
  worker <- function(df_chr) {
    bridges:::tibble_to_matrix(df_chr, value_column = allele)
  }
  
  all_input_Xs <- lapply(split_by_chr, worker)
  X = all_input_Xs %>% do.call("cbind", .)
  X
}

get_labels = function(sample_id) {
  labels_df = fread(paste0("signatures_dataset/clones_trees/",sample_id,"_clones.tsv"))
  labels_df
}

get_input = function(sample_id) {
  df = fread(paste0("signatures_dataset/DLP/CNA/persample/",sample_id,"_hscn.csv.gz"))
  df = df %>% 
    dplyr::select(cell_id, chr, start, end, state, Maj, Min) %>% 
    dplyr::rename(A = Maj, B = Min, CN = state)  
}