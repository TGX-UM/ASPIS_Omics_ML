library(httr)
library(dplyr)
library(purrr)
library(tidyr)
library(minervar)

if(packageVersion("minervar") < "0.8.3") {
  stop("Please install package:minervar 0.8.3")
}

#remotes::install_gitlab(repo = "uniluxembourg/lcsb/BioCore/disease_maps/tmcuration", host = "gitlab.com", dependencies = TRUE)

#########
### Load the bioentities from the map
components <- minervar::get_map_components("https://ontox.elixir-luxembourg.org/minerva/api/", project_id = "Liver_Lipid_Metabolism_Physiological_Map_August_2024")
bioentities <- components$map_elements[[1]]

# Access data
map_ids <- as.character(552:559)
full_elements_list <- map_df(map_ids, ~components$map_elements[[.x]])
full_reactions_list <- map_df(map_ids, ~components$map_reactions[[.x]])


# Using base R
filtered_df <- full_elements_list[full_elements_list$type %in% c("Protein", "Gene", "RNA"), ]

unique_symbols <- unique(filtered_df$symbol[!is.na(filtered_df$symbol)])

write.table(unique_symbols, file = "PM_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



