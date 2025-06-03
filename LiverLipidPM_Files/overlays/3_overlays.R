# Load the CSV file
converted_gene_list <- read.csv("converted_gene_lists.csv")

# Define unique hex colors for each column
colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33F5", "#F5FF33", 
            "#33FFF5", "#F533FF", "#57FF33", "#5733FF", "#FF5733",
            "#FFB533", "#B533FF", "#33FFB5", "#B5FF33", "#FF33B5")

# Get column names
column_names <- names(converted_gene_list)

# Process each column
for (i in 1:length(column_names)) {
  col_name <- column_names[i]
  
  # Get values from the column, removing NAs
  values <- converted_gene_list[[col_name]]
  values <- values[!is.na(values)]
  
  # Create data frame with name and color columns
  output_df <- data.frame(
    name = values,
    color = rep(colors[i], length(values))
  )
  
  # Create filename
  filename <- paste0("overlays/", col_name, ".txt")
  
  # Write comment and header first
  cat("# NAME=", col_name, "\n", file = filename)
  
  # Write the data with header
  write.table(output_df, file = filename, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", append = TRUE)}

# Print confirmation
cat("Files created for", length(column_names), "columns in the overlays folder\n")
