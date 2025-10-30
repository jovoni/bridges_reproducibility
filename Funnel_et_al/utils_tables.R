
# Method 1: Custom function approach
library(dplyr)

# Function to bold minimum values in each row
bold_min_latex <- function(df, numeric_cols = NULL) {

  if (is.null(numeric_cols)) {
    # Automatically detect numeric columns (excluding first column assumed to be ID)
    numeric_cols <- names(df)[sapply(df, is.numeric)]
  }

  # Create a copy for modification
  df_latex <- df

  # Process each row
  for (i in 1:nrow(df)) {
    # Get numeric values for this row (excluding NAs)
    row_values <- as.numeric(df[i, numeric_cols])
    valid_values <- row_values[!is.na(row_values)]

    if (length(valid_values) > 0) {
      min_val <- min(valid_values, na.rm = TRUE)

      # Find which column(s) have the minimum value and bold them
      for (col in numeric_cols) {
        if (!is.na(df[i, col]) && df[i, col] == min_val) {
          df_latex[i, col] <- paste0("\\textbf{", df[i, col], "}")
        }
      }
    }
  }

  # Convert NA to dashes
  df_latex[is.na(df_latex)] <- "-"

  return(df_latex)
}

# Function to generate LaTeX table
generate_latex_table <- function(df, caption = "", label = "",
                                 col_names = NULL, table_spec = NULL) {

  if (is.null(col_names)) {
    col_names <- names(df)
  }

  if (is.null(table_spec)) {
    # Auto-generate table specification
    n_cols <- ncol(df)
    table_spec <- paste0("l|", paste(rep("r", n_cols - 1), collapse = ""))
  }

  # Start building the LaTeX code
  latex_code <- paste0(
    "\\begin{table}[ht]\n",
    "    \\centering\n",
    "    \\begin{tabular}{", table_spec, "}\n",
    "    \\toprule\n"
  )

  # Add header
  header <- paste(col_names, collapse = " & ")
  latex_code <- paste0(latex_code, "      ", header, " \\\\\n")
  latex_code <- paste0(latex_code, "      \\midrule\n")

  # Add data rows
  for (i in 1:nrow(df)) {
    row_data <- paste(df[i, ], collapse = " & ")
    latex_code <- paste0(latex_code, "    ", row_data, " \\\\\n")
  }

  # Close table
  latex_code <- paste0(
    latex_code,
    "    \\bottomrule\n",
    "    \\end{tabular}\n"
  )

  if (caption != "") {
    latex_code <- paste0(latex_code, "    \\caption{", caption, "}\n")
  }

  if (label != "") {
    latex_code <- paste0(latex_code, "    \\label{", label, "}\n")
  }

  latex_code <- paste0(latex_code, "\\end{table}")

  return(latex_code)
}

# Generate the complete LaTeX table
latex_table <- generate_latex_table(
  df_with_bold,
  caption = "Sample data with minimum values highlighted in bold.",
  label = "tab:sample_data",
  col_names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
  table_spec = "l|rrr"
)

# Print the result
cat(latex_table)

# =====================================================
# Method 2: Using kableExtra (more powerful)
# =====================================================

library(kableExtra)
library(knitr)

# Function using kableExtra
create_latex_table_kable <- function(df, numeric_cols = NULL) {

  if (is.null(numeric_cols)) {
    numeric_cols <- names(df)[sapply(df, is.numeric)]
  }

  # Create formatting matrix
  bold_matrix <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df))
  colnames(bold_matrix) <- names(df)

  # Find minimum values for each row
  for (i in 1:nrow(df)) {
    row_values <- as.numeric(df[i, numeric_cols])
    valid_values <- row_values[!is.na(row_values)]

    if (length(valid_values) > 0) {
      min_val <- min(valid_values, na.rm = TRUE)

      for (col in numeric_cols) {
        if (!is.na(df[i, col]) && df[i, col] == min_val) {
          bold_matrix[i, col] <- TRUE
        }
      }
    }
  }

  # Create the table
  df_display <- df
  df_display[is.na(df_display)] <- "-"

  kable_table <- df_display %>%
    kbl(format = "latex",
        booktabs = TRUE,
        col.names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
        escape = FALSE) %>%
    kable_styling(latex_options = c("hold_position")) %>%
    column_spec(1, border_right = TRUE)

  # Apply bold formatting
  for (i in 1:nrow(bold_matrix)) {
    for (j in 1:ncol(bold_matrix)) {
      if (bold_matrix[i, j]) {
        kable_table <- kable_table %>%
          cell_spec(i, j, bold = TRUE)
      }
    }
  }

  return(kable_table)
}

# Alternative kableExtra approach
kable_result <- df %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), "-", as.character(.)))) %>%
  kbl(format = "latex",
      booktabs = TRUE,
      col.names = c("Sample ID", "Bridges", "Sitka", "Lazac"),
      escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  column_spec(1, border_right = TRUE) %>%
  # You would need to manually specify which cells to bold
  # or write a more complex function to automate this
  cell_spec(1, 3, bold = TRUE) %>%  # DG1134 - Sitka (22)
  cell_spec(2, 2, bold = TRUE) %>%  # DG1197 - Bridges (2)
  cell_spec(3, 2, bold = TRUE)      # OV2295 - Bridges (57)
# ... continue for all minimum values

# Print kable result
print(kable_result)

# =====================================================
# Method 3: Using xtable with custom formatting
# =====================================================

library(xtable)

# Function using xtable
create_latex_xtable <- function(df, numeric_cols = NULL) {

  # Process data same as before
  df_processed <- bold_min_latex(df, numeric_cols)
  df_processed[is.na(df_processed)] <- "-"

  # Create xtable
  xt <- xtable(df_processed,
               caption = "Sample data with minimum values in bold",
               label = "tab:sample")

  # Print with custom options
  latex_output <- print(xt,
                        type = "latex",
                        include.rownames = FALSE,
                        sanitize.text.function = identity,  # Don't escape LaTeX commands
                        hline.after = c(-1, 0, nrow(df_processed)),
                        add.to.row = list(
                          pos = list(0),
                          command = "\\midrule\n"
                        ),
                        print.results = FALSE)

  # Replace default rules with booktabs
  latex_output <- gsub("\\\\hline", "\\\\toprule", latex_output, fixed = TRUE)
  latex_output <- gsub("\\\\toprule\n\\\\toprule", "\\\\toprule", latex_output, fixed = TRUE)
  latex_output <- gsub("\\\\hline\n\\\\end\\{tabular\\}", "\\\\bottomrule\n\\\\end{tabular}", latex_output)

  return(latex_output)
}

# Usage example
xtable_result <- create_latex_xtable(df, c("bridges", "sitka", "lazac"))
cat(xtable_result)

# =====================================================
# Method 4: Most flexible custom solution
# =====================================================

create_complete_latex_solution <- function(df, id_col = 1,
                                           caption = "", label = "",
                                           col_names = NULL) {

  # Identify numeric columns
  numeric_cols <- names(df)[-id_col]
  numeric_cols <- numeric_cols[sapply(df[numeric_cols], is.numeric)]

  # Create result dataframe
  result_df <- df

  # Bold minimum values
  for (i in 1:nrow(df)) {
    row_values <- as.numeric(df[i, numeric_cols])
    valid_values <- row_values[!is.na(row_values)]

    if (length(valid_values) > 0) {
      min_val <- min(valid_values)

      for (col in numeric_cols) {
        if (!is.na(df[i, col]) && df[i, col] == min_val) {
          result_df[i, col] <- paste0("\\textbf{", df[i, col], "}")
        }
      }
    }
  }

  # Handle missing values
  result_df[is.na(result_df)] <- "-"

  # Generate LaTeX
  if (is.null(col_names)) {
    col_names <- names(result_df)
  }

  # Auto table spec
  table_spec <- paste0("l|", paste(rep("r", length(numeric_cols)), collapse = ""))

  latex_lines <- c(
    "\\begin{table}[ht]",
    "    \\centering",
    paste0("    \\begin{tabular}{", table_spec, "}"),
    "    \\toprule",
    paste0("      ", paste(col_names, collapse = " & "), " \\\\"),
    "      \\midrule"
  )

  # Add data rows
  for (i in 1:nrow(result_df)) {
    row_line <- paste0("    ", paste(result_df[i, ], collapse = " & "), " \\\\")
    latex_lines <- c(latex_lines, row_line)
  }

  # Close table
  latex_lines <- c(
    latex_lines,
    "    \\bottomrule",
    "    \\end{tabular}"
  )

  if (caption != "") {
    latex_lines <- c(latex_lines, paste0("    \\caption{", caption, "}"))
  }

  if (label != "") {
    latex_lines <- c(latex_lines, paste0("    \\label{", label, "}"))
  }

  latex_lines <- c(latex_lines, "\\end{table}")

  return(paste(latex_lines, collapse = "\n"))
}

# Final usage
final_table <- create_complete_latex_solution(
  df,
  caption = "Data comparison with minimum values highlighted in bold.",
  label = "tab:comparison",
  col_names = c("Sample ID", "Bridges", "Sitka", "Lazac")
)

cat(final_table)
