# Convert rat genes to human orthologs and then Ensembl IDs to HGNC symbols

# Load the source CSV table
data <- read.csv("data.csv", stringsAsFactors = FALSE)

# Install packages if needed 
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("org.Hs.eg.db", "org.Rn.eg.db", "biomaRt"))

library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(biomaRt)

# Get all unique IDs from the dataframe
all_ids <- unique(unlist(data, use.names = FALSE))
all_ids <- all_ids[!is.na(all_ids)]

# Function to identify rat vs human Ensembl IDs
# Rat IDs typically start with ENSRNOG, human with ENSG
identify_species <- function(ids) {
  rat_pattern <- "^ENSRNOG"
  human_pattern <- "^ENSG"
  
  species <- rep("unknown", length(ids))
  species[grepl(rat_pattern, ids)] <- "rat"
  species[grepl(human_pattern, ids)] <- "human"
  
  return(species)
}

species_info <- identify_species(all_ids)
rat_ids <- all_ids[species_info == "rat"]
human_ids <- all_ids[species_info == "human"]

cat("Found", length(rat_ids), "rat IDs and", length(human_ids), "human IDs\n")

# Step 1: Convert rat Ensembl IDs to human orthologs using biomaRt
rat_to_human_mapping <- data.frame()

if(length(rat_ids) > 0) {
  cat("Converting rat genes to human orthologs...\n")
  
  # Try different Ensembl mirrors if main site is down
  mirrors <- c("https://www.ensembl.org",
               "https://uswest.ensembl.org", 
               "https://useast.ensembl.org",
               "https://asia.ensembl.org")
  
  rat_mart <- NULL
  for(mirror in mirrors) {
    tryCatch({
      cat("Trying mirror:", mirror, "\n")
      rat_mart <- useEnsembl("ensembl", 
                             dataset = "rnorvegicus_gene_ensembl",
                             mirror = mirror)
      break
    }, error = function(e) {
      cat("Mirror failed, trying next...\n")
    })
  }
  
  if(!is.null(rat_mart)) {
    # Get human orthologs for rat genes
    rat_to_human <- getBM(
      attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
      filters = "ensembl_gene_id",
      values = rat_ids,
      mart = rat_mart
    )
    
    # Remove empty orthologs
    rat_to_human <- rat_to_human[rat_to_human$hsapiens_homolog_ensembl_gene != "", ]
    rat_to_human_mapping <- rat_to_human
    
    cat("Found", nrow(rat_to_human_mapping), "rat-to-human ortholog mappings\n")
  } else {
    cat("All Ensembl mirrors failed. Using org.Rn.eg.db alternative method...\n")
    
    # Alternative: Use org.Rn.eg.db to get rat gene symbols, then map to human
    tryCatch({
      # Get rat gene symbols
      rat_symbols <- AnnotationDbi::select(org.Rn.eg.db,
                                           keys = rat_ids,
                                           columns = "SYMBOL",
                                           keytype = "ENSEMBL")
      
      # Get human genes with same symbols (approximate ortholog mapping)
      if(nrow(rat_symbols) > 0) {
        rat_symbols <- rat_symbols[!is.na(rat_symbols$SYMBOL), ]
        human_orthologs <- AnnotationDbi::select(org.Hs.eg.db,
                                                 keys = rat_symbols$SYMBOL,
                                                 columns = "ENSEMBL",
                                                 keytype = "SYMBOL")
        
        # Create mapping dataframe
        if(nrow(human_orthologs) > 0) {
          human_orthologs <- human_orthologs[!is.na(human_orthologs$ENSEMBL), ]
          # Match rat IDs to human IDs via symbols
          merged <- merge(rat_symbols, human_orthologs, by = "SYMBOL")
          rat_to_human_mapping <- data.frame(
            ensembl_gene_id = merged$ENSEMBL.x,
            hsapiens_homolog_ensembl_gene = merged$ENSEMBL.y
          )
          cat("Found", nrow(rat_to_human_mapping), "rat-to-human mappings via gene symbols\n")
        }
      }
    }, error = function(e) {
      cat("Alternative method also failed. Proceeding with human genes only.\n")
    })
  }
}

# Step 2: Get all human Ensembl IDs (original + converted from rat)
all_human_ids <- c(human_ids)
if(nrow(rat_to_human_mapping) > 0) {
  all_human_ids <- c(all_human_ids, rat_to_human_mapping$hsapiens_homolog_ensembl_gene)
}
all_human_ids <- unique(all_human_ids)

# Step 3: Convert human Ensembl IDs to HGNC symbols
cat("Converting human Ensembl IDs to HGNC symbols...\n")

ensembl_to_symbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = all_human_ids,
                                           columns = "SYMBOL", 
                                           keytype = "ENSEMBL")

# Create lookup tables
human_lookup <- setNames(ensembl_to_symbol$SYMBOL, ensembl_to_symbol$ENSEMBL)

# Create rat-to-symbol lookup (rat -> human -> symbol)
rat_lookup <- character()
if(nrow(rat_to_human_mapping) > 0) {
  for(i in 1:nrow(rat_to_human_mapping)) {
    rat_id <- rat_to_human_mapping$ensembl_gene_id[i]
    human_id <- rat_to_human_mapping$hsapiens_homolog_ensembl_gene[i]
    if(human_id %in% names(human_lookup) && !is.na(human_lookup[human_id])) {
      rat_lookup[rat_id] <- human_lookup[human_id]
    }
  }
}

# Step 4: Function to convert any ID (rat or human) to HGNC symbol
convert_to_hgnc <- function(x) {
  result <- character(length(x))
  
  for(i in seq_along(x)) {
    if(is.na(x[i])) {
      result[i] <- NA
    } else if(grepl("^ENSRNOG", x[i])) {
      # Rat gene
      if(x[i] %in% names(rat_lookup)) {
        result[i] <- rat_lookup[x[i]]
      } else {
        result[i] <- x[i]  # Keep original if no conversion
      }
    } else if(grepl("^ENSG", x[i])) {
      # Human gene
      if(x[i] %in% names(human_lookup) && !is.na(human_lookup[x[i]])) {
        result[i] <- human_lookup[x[i]]
      } else {
        result[i] <- x[i]  # Keep original if no conversion
      }
    } else {
      result[i] <- x[i]  # Keep as is if not recognized pattern
    }
  }
  
  return(result)
}

# Step 5: Apply conversion to all columns
data_converted <- data.frame(lapply(data, convert_to_hgnc), stringsAsFactors = FALSE)

# Results summary
cat("\n=== Conversion Summary ===\n")
cat("Original rat IDs:", length(rat_ids), "\n")
cat("Rat orthologs found:", length(rat_lookup), "\n")
cat("Rat conversion rate:", round(length(rat_lookup)/length(rat_ids)*100, 1), "%\n")
cat("Original human IDs:", length(human_ids), "\n")
cat("Human symbols found:", sum(!is.na(human_lookup)), "\n")
cat("Human conversion rate:", round(sum(!is.na(human_lookup))/length(human_ids)*100, 1), "%\n")

# Show unmapped genes
unmapped_rat <- rat_ids[!rat_ids %in% names(rat_lookup)]
unmapped_human <- human_ids[!human_ids %in% names(human_lookup) | is.na(human_lookup[human_ids])]

if(length(unmapped_rat) > 0) {
  cat("\nUnmapped rat genes (", length(unmapped_rat), "):\n")
  cat(paste(head(unmapped_rat, 10), collapse = ", "))
  if(length(unmapped_rat) > 10) cat(", ...")
  cat("\n")
}

if(length(unmapped_human) > 0) {
  cat("\nUnmapped human genes (", length(unmapped_human), "):\n")
  cat(paste(head(unmapped_human, 10), collapse = ", "))
  if(length(unmapped_human) > 10) cat(", ...")
  cat("\n")
}

# Count how many cells in final dataframe still contain Ensembl IDs
ensembl_pattern <- "^ENS[RG]"
remaining_ensembl <- sum(grepl(ensembl_pattern, unlist(data_converted), perl = TRUE), na.rm = TRUE)
total_non_na <- sum(!is.na(unlist(data_converted)))

cat("\nFinal dataframe status:\n")
cat("Total non-NA values:", total_non_na, "\n")
cat("Still Ensembl IDs:", remaining_ensembl, "\n")
cat("Successfully converted:", total_non_na - remaining_ensembl, "\n")
cat("Overall conversion rate:", round((total_non_na - remaining_ensembl)/total_non_na*100, 1), "%\n")

print("\nFirst few rows of converted dataframe:")
print(head(data_converted))

# Save unmapped genes to files for manual inspection
if(length(unmapped_rat) > 0) {
  write.csv(data.frame(unmapped_rat_genes = unmapped_rat), 
            "unmapped_rat_genes.csv", row.names = FALSE)
  cat("\nUnmapped rat genes saved to 'unmapped_rat_genes.csv'\n")
}

if(length(unmapped_human) > 0) {
  write.csv(data.frame(unmapped_human_genes = unmapped_human), 
            "unmapped_human_genes.csv", row.names = FALSE)
  cat("Unmapped human genes saved to 'unmapped_human_genes.csv'\n")
}


# Load manual mappings and complete the conversion
manual_mappings <- read.csv("unmapped_gene_symbols.txt", stringsAsFactors = FALSE)

# Create lookup table from manual mappings
manual_lookup <- setNames(manual_mappings$symbol, manual_mappings$unmapped_genes)

# Function to apply final conversion including manual mappings
final_convert <- function(x) {
  result <- character(length(x))
  
  for(i in seq_along(x)) {
    if(is.na(x[i])) {
      result[i] <- NA
    } else if(x[i] %in% names(manual_lookup)) {
      # Use manual mapping
      result[i] <- manual_lookup[x[i]]
    } else {
      # Keep the value from previous conversion
      result[i] <- x[i]
    }
  }
  
  return(result)
}

# Apply final conversion to all columns
data_final <- data.frame(lapply(data_converted, final_convert), stringsAsFactors = FALSE)

# Check for remaining unmapped genes in each column
cat("\n=== Remaining Unmapped Genes by Column ===\n")
ensembl_pattern <- "^ENS[RG]"

for(col_name in names(data_final)) {
  unmapped_in_col <- data_final[[col_name]][grepl(ensembl_pattern, data_final[[col_name]], perl = TRUE)]
  unmapped_in_col <- unmapped_in_col[!is.na(unmapped_in_col)]
  
  if(length(unmapped_in_col) > 0) {
    cat("\nColumn '", col_name, "' - ", length(unmapped_in_col), " unmapped genes:\n", sep = "")
    unique_unmapped <- unique(unmapped_in_col)
    cat(paste(unique_unmapped, collapse = ", "))
    cat("\n")
  } else {
    cat("Column '", col_name, "': All genes successfully mapped ✓\n", sep = "")
  }
}

# Final summary statistics
total_cells <- sum(!is.na(unlist(data_final)))
remaining_ensembl_final <- sum(grepl(ensembl_pattern, unlist(data_final), perl = TRUE), na.rm = TRUE)

cat("\n=== Final Summary ===\n")
cat("Total non-NA values:", total_cells, "\n")
cat("Successfully converted to symbols:", total_cells - remaining_ensembl_final, "\n")
cat("Still unmapped (Ensembl IDs):", remaining_ensembl_final, "\n")
cat("Final conversion rate:", round((total_cells - remaining_ensembl_final)/total_cells*100, 1), "%\n")

print("\nFirst few rows of final dataframe:")
print(head(data_final))

# Save the final converted dataframe
write.csv(data_final, "converted_gene_lists.csv", row.names = FALSE)
cat("\nFinal converted dataframe saved to 'converted_gene_lists.csv'\n")
