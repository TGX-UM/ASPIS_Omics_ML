#========================================================#
# Machine learning classification of steatotic compounds #
#========================================================#

# Below are useful tutorials on tidymodels

# 1. https://lsinks.github.io/posts/2023-04-10-tidymodels/tidymodels_tutorial.html
# 2. https://www.tmwr.org/workflow-sets.html
# 3. https://rebeccabarter.com/blog/2020-03-25_machine_learning
# 4. https://juliasilge.com/blog/spam-email/

# More info on recipes can be found here: https://recipes.tidymodels.org/reference/recipe.html

# If working on linux terminal, launch R using: R --max-ppsize=500000
# For windows users, use options(expressions = 200000)

library(tidyverse) 
library(tidymodels) 
library(glmnet) 
library(themis) 
library(ggfortify)
library(doParallel)
library(vip)
library(tidyr)
library(viridis)
library(stringr)
library(patchwork)
library(gridExtra)

setwd("~/Research/Data/ASPIS_OMICS/Ensembl_ID/Human/")

output_dir <- "~/Research/Figures_and_Results/ASPIS_OMICs_ML"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load files and perform column manipulation
Data <- read.delim("NormData_TG_Gates_Human_In_Vitro_Ens.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
metadata <- read.delim("TG_Human_In_Vitro_metadata.tsv", header = TRUE, sep = "\t", row.names=1)

Data <- Data[,order(colnames(Data))]
metadata <- metadata[order(rownames(metadata)),]
colnames(Data) <- metadata$Suggested_name

# Count the number of genes before filtering
count_before <- nrow(Data)

# Apply the filter to remove genes with low intensity
Data_filtered <- Data[rowMeans(Data) > 6, ]

# Count the number of genes after filtering
count_after <- nrow(Data_filtered)

# Calculate the number of genes removed
count_removed <- count_before - count_after

# Transpose the dataframe & add Steatogenic column
T_Data <- as.data.frame(x = t(Data_filtered), stringsAsFactors = FALSE)
T_Data <- cbind(T_Data, metadata$Steatogenic)
colnames(T_Data)[ncol(T_Data)] <- "Steatogenic"
T_Data$Steatogenic <- as.factor(T_Data$Steatogenic)

#====================#
# Summary stats plot #
#====================#

# Define the output directory for Summary stats
summary_stats_dir <- file.path(output_dir, "Summary_stats")

if (!dir.exists(summary_stats_dir)) {
  dir.create(summary_stats_dir, recursive = TRUE)
}

# Bar plot for number of genes before and after mean filter
plot_data_max <- data.frame(
  Status = factor(c("Before filter", "Removed by filter", "After filter"), 
                  levels = c("Before filter", "Removed by filter", "After filter")),
  Count = c(count_before, count_removed, count_after)
)

barplot <- ggplot(plot_data_max, aes(x = Status, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Before filter" = "blue", "Removed by filter" = "red", "After filter" = "darkgreen")) +
  labs(
    #title = "Number of genes before and after mean filter",
    x = "Filter Status",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    #plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18), 
    axis.text.x = element_text(angle = 0, hjust = 0.8, size = 14), 
    axis.text.y = element_text(size = 14), 
    legend.position = "none"
  ) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5)

print(barplot)
ggsave(filename = file.path(summary_stats_dir, "Number_of_genes_Barplot.pdf"), plot = barplot, width = 8, height = 6, dpi = 1000)

# PCA of filtered data
filtered_pca <- prcomp(T_Data[,-ncol(T_Data)])

# Calculate the percentage of variance explained by each principal component
explained_variance <- filtered_pca$sdev^2 / sum(filtered_pca$sdev^2) * 100

# Create a data frame for PCA plot
pca_df <- as.data.frame(filtered_pca$x)
pca_df$Steatogenic <- T_Data$Steatogenic

# Define custom shapes & colors for Steatogenic
steatogenic_colors <- c("No" = "darkgreen", "Yes" = "red")
steatogenic_shapes <- c("No" = 16, "Yes" = 17)  # 16 for circle, 17 for triangle

# Include time on scatter plot
pca_df$Time_Hours <- factor(metadata$Factor.Value.time., levels = c(2, 8, 24))

# Define custom colors for each time point
time_colors <- c("2" = "darkgreen", "8" = "blue", "24" = "red")

# After mean filter PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Time_Hours, shape = Steatogenic)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = time_colors) +  # Use manual colors for each time point
  scale_shape_manual(values = steatogenic_shapes) +  
  labs(title = "PCA: After Mean Filter",
       x = paste0("PC 1:", signif(explained_variance[1], 2), "%"),
       y = paste0("PC 2:", signif(explained_variance[2], 2), "%"),
       color = "Time (Hours)") +  # Label the color legend as Time (Hours)
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14)  
  )

print(pca_plot)
ggsave(filename = file.path(summary_stats_dir, "PCA_Plot_after_mean_filter.pdf"), plot = pca_plot, width = 8, height = 6, dpi = 1000)

### ====== End of summary stats ====== ###

#==========================#
# Machine Learning Section #
#==========================#

# Define the output directory for ML_results
ml_results_dir <- file.path(output_dir, "ML_results")

if (!dir.exists(ml_results_dir)) {
  dir.create(ml_results_dir, recursive = TRUE)
}
# For classification problem in tidymodels, outcome should be a factor
# NB: Tidymodels treats the first level as the event
levels(T_Data$Steatogenic)

# In our case, the first level is No, so we need to reorder the events
T_Data$Steatogenic <- relevel(T_Data$Steatogenic, ref = "Yes")
levels(T_Data$Steatogenic)

# Step 1: Split data into training and testing set

# The rsample package is used to create splits and cross validation folds from the data
# The resulting object is called an rsplit object, & it contains the 
# original data & information about whether a record goes to testing or training
# NB: This object is a nested list

#set.seed(123)
#Data_split <- initial_split(T_Data, strata = Steatogenic, prop = 3/4)

#Data_train <- training(Data_split)
#Data_test <- testing(Data_split)

set.seed(456)
Data_CV <- bootstraps(T_Data, times = 1000, strata = Steatogenic)
Data_CV 

# Supplementary: Visualize bootstraps
# Tidy Data_CV and extract numeric values correctly
Data_CV_boots <- tidy(Data_CV)

# Extract numeric part of the Resample and handle conversion safely
Data_CV_boots$Resample_Number <- as.numeric(str_extract(Data_CV_boots$Resample, "\\d+"))

# Remove rows with missing values to ensure proper plotting
Data_CV_boots <- na.omit(Data_CV_boots)

# Plot the heatmap 
heatmap <- ggplot(Data_CV_boots, aes(x = Resample_Number, y = Row, fill = Data)) +
  geom_tile(color = "white") +  
  scale_fill_viridis(discrete = TRUE, option = "C") + 
  labs(
    title = "Heatmap of 1000 Bootstraps",
    x = "Resample Number",
    y = "Number of Samples",
    fill = "Data Type"
  ) +
  theme_minimal(base_size = 12) +
  scale_x_continuous(breaks = seq(0, 1000, by = 50)) +  
  scale_y_continuous(breaks = seq(4, 108, by = 4)) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),  
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

print(heatmap)
ggsave(filename = file.path(summary_stats_dir, "Heatmap_of_Bootstraps.pdf"), plot = heatmap, width = 14, height = 12, dpi = 1000)

#==================== End of Heatmap =============================#
#==================================
# Step 2: Create a recipe

# This step bundles the formula, data and feature engineering
# steps into a recipe object
# All the feature engineering steps have the form step_*() 
# Feature engineering steps for each package are available at: https://www.tmwr.org/pre-proc-table.html
# Recipe consists of the formula (outcome ~ predictors)

Data_Recipe <- 
  recipe(Steatogenic ~ ., data = T_Data) %>%
  step_nzv(all_predictors()) %>%
  step_smote(Steatogenic) 

# Confirm that SMOTE fixes class imbalance 
Data_Recipe %>% prep() %>% bake(new_data = NULL) %>% count(Steatogenic) 

# Check how many features remain for model training
Data_Recipe %>% prep() %>% juice() %>% ncol()

#===========================================#
# Optional Step: PCA to check SMOTE Samples #
#===========================================#
# NB: This step adds an additional column for labelling original & SMOTE Data. If applied
# Check if it might impact before training the models 
# PCA after SMOTE will provide slightly varying scatterplots because we simulate samples randomly 

PCA_Recipe <- 
  recipe(Steatogenic ~ ., data = T_Data) %>%
  step_nzv(all_predictors()) %>%
  step_smote(Steatogenic) %>%
  prep()

# Apply the recipe to the data
prepped_data <- bake(PCA_Recipe, new_data = NULL)

# Identify original and SMOTE-generated samples
original_sample_count <- nrow(T_Data)
prepped_data$Sample_Type <- c(rep("Original", original_sample_count), rep("SMOTE", nrow(prepped_data) - original_sample_count))

# Perform PCA on the prepped data (after SMOTE)
# Exclude the last two columns (Steatogenic and Sample_Type)
pca_smote <- prcomp(prepped_data[, !(names(prepped_data) %in% c("Steatogenic", "Sample_Type"))])

# Calculate the percentage of variance explained by each principal component
explained_variance_smote <- pca_smote$sdev^2 / sum(pca_smote$sdev^2) * 100

# Create a data frame for PCA plot
pca_df_smote <- as.data.frame(pca_smote$x)
pca_df_smote$Steatogenic <- prepped_data$Steatogenic
pca_df_smote$Sample_Type <- prepped_data$Sample_Type

# Define custom colors for Sample_Type
sample_type_colors <- c("Original" = "darkgreen", "SMOTE" = "orange")

# Plot PCA using ggplot2 with the specified aesthetics
pca_plot_smote <- ggplot(pca_df_smote, aes(x = PC1, y = PC2, shape = Steatogenic, color = Sample_Type)) + 
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = sample_type_colors) +
  scale_shape_manual(values = steatogenic_shapes) +
  labs(title = "PCA: After SMOTE",
       x = paste0("PC 1:", signif(explained_variance_smote[1], 2), "%"),
       y = paste0("PC 2:", signif(explained_variance_smote[2], 2), "%")) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12)  
  ) +
  guides(
    color = guide_legend(title = "Sample Type"),
    shape = guide_legend(title = "Steatogenic")
  )

print(pca_plot_smote)
ggsave(filename = file.path(summary_stats_dir, "PCA_Plot_SMOTE.pdf"), plot = pca_plot_smote, width = 8, height = 6, dpi = 1000)

#========= End of optional step: SMOTE PCA ====================# 

# The functions prep(), bake(), and juice() can be used to apply the recipe object to the dataset
# Here, we use workflow procedure which handles these steps automatically
# The workflows package allow us to create bundles of models and fits

# Step 3: Define ML model specifications with parsnip
# The list of parsnip engines can be found here: https://www.tidymodels.org/find/parsnip/ 

# i. Elastic net regularization of logistic regression (ENLR)
# ENLR differs from logistic regression in that it has 2 hyperparameters that are tuned: penalty & mixture
glmnet_spec <-
  logistic_reg(penalty = tune(),
               mixture = tune()) %>%
  set_engine("glmnet")
glmnet_spec

# ii. Random forest model
rf_spec <- rand_forest() %>%
  set_args(mtry = sqrt(ncol(T_Data))) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")
rf_spec

# iii. XgBoost model
xgb_spec <-
  boost_tree(trees = 500, min_n = 10,
             mtry = sqrt(ncol(T_Data)),
             learn_rate = 0.01) %>%
  set_engine("xgboost") %>%
  set_mode("classification")
xgb_spec

# iv. KNN model
knn_spec <- nearest_neighbor() %>%
  set_engine("kknn") %>%
  set_mode("classification")
knn_spec

# vi. SVM linear model
svm_spec <- svm_linear() %>% 
  set_engine("kernlab", prob.model = TRUE) %>% 
  set_mode("classification")
svm_spec

# Specify metrics to be assessed for each of the models (uses yardstick package)
# Define metric sets for models
metrics <- metric_set(accuracy, sensitivity, specificity, roc_auc, j_index, mcc)

# Create a workflow set for each model
# This step matches each model to the recipe
# tunes them, then evaluates their performance efficiently
# Workflow sets take named lists of preprocessors and model specifications and combine them into an object containing multiple workflows
ML_models <- 
  workflow_set(
    preproc = list(Data_Recipe),
    models = list(
      glmnet = glmnet_spec,
      rf = rf_spec,
      xgb = xgb_spec,
      knn = knn_spec,
      svm = svm_spec  
    ),
    cross = TRUE # Useful if multiple preprocessing recipes are defined 
  )

ML_models

# Run all the models
# registerDoParallel(cores = detectCores() - 2)

ML_results <- 
  ML_models %>%
  workflow_map(
    "tune_grid",
    resamples = Data_CV,
    metrics = metrics,
    control = control_resamples(save_pred = TRUE),
    verbose = TRUE
  )

ML_results %>% collect_predictions()
ML_results %>% collect_metrics()

#===============#
# RoC-AuC Curve #
#===============#
# Define the label order of the legend
legend_order <- toupper(c("svm", "rf", "ENLR", "xgb", "knn"))

# Collect predictions
roc_predictions <- ML_results %>%
  mutate(preds = map(result, collect_predictions)) %>%
  select(wflow_id, preds) %>%
  unnest(preds)

# Generate ROC data for each model
roc_data <- roc_predictions %>%
  group_by(wflow_id) %>%
  roc_curve(truth = Steatogenic, .pred_Yes) %>%
  ungroup() %>%
  mutate(wflow_id = str_replace_all(wflow_id, c("recipe_" = "", "glmnet" = "ENLR"))) %>%
  mutate(wflow_id = toupper(wflow_id))  # Ensure all workflow names are uppercase

# Calculate AUC values and sort in ascending order
auc_values <- roc_predictions %>%
  group_by(wflow_id) %>%
  roc_auc(truth = Steatogenic, .pred_Yes) %>%
  ungroup() %>%
  mutate(wflow_id = str_replace_all(wflow_id, c("recipe_" = "", "glmnet" = "ENLR"))) %>%
  mutate(wflow_id = toupper(wflow_id)) %>%  
  arrange(.estimate)

# Ensure the factor levels are the same and in the desired order in both data frames
roc_data <- roc_data %>%
  mutate(wflow_id = factor(wflow_id, levels = legend_order))

auc_values <- auc_values %>%
  mutate(wflow_id = factor(wflow_id, levels = legend_order))

# Define colors for the models
model_colors <- setNames(
  c("purple", "orange", "red", "brown", "blue", "green"),  
  legend_order
)

# Plot ROC curves
roc_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = wflow_id)) +
  geom_path(linewidth = 1.5) +  
  geom_abline(lty = 2, color = "gray80", linewidth = 1.5) +
  coord_equal() +
  labs(
    title = "ROC AuC Curves of the different Models",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),   
    plot.title = element_text(hjust = 0.5)   
  ) +
  scale_color_manual(values = model_colors) +  
  geom_text(
    data = auc_values, 
    aes(
      x = 0.90,  
      y = 0.05 + seq(0, 0.05 * (nrow(auc_values) - 1), by = 0.05),  
      label = sprintf("%.3f", .estimate),  # Plot AUC values to three decimal places
      color = wflow_id
    ), 
    size = 5,
    hjust = 0.5  
  ) +
  annotate("text", x = 0.90, y = 0.05 + 0.05 * nrow(auc_values), label = "AUC", color = "black", size = 5, fontface = "bold", hjust = 0.5)

print(roc_plot)
ggsave(filename = file.path(ml_results_dir, "RoC_AuC.pdf"), plot = roc_plot, width = 8, height = 6, dpi = 1000)

#========= End of RoC Curve ====================#

#===================================================#
# Optional step: Rank all models with their metrics #
#===================================================#

model_ranking <- roc_predictions %>%
  group_by(wflow_id) %>%
  summarise(
    accuracy = accuracy_vec(truth = Steatogenic, estimate = .pred_class),
    sensitivity = sens_vec(truth = Steatogenic, estimate = .pred_class),
    specificity = spec_vec(truth = Steatogenic, estimate = .pred_class),
    roc_auc = roc_auc_vec(truth = Steatogenic, estimate = .pred_Yes),
    j_index = sens_vec(truth = Steatogenic, estimate = .pred_class) + 
              spec_vec(truth = Steatogenic, estimate = .pred_class) - 1,
    mcc = mcc_vec(truth = Steatogenic, estimate = .pred_class)
  ) %>%
  arrange(desc(roc_auc))   

print(model_ranking)

# Save the ranking to a file
write.table(model_ranking, file = file.path(ml_results_dir, "Ranking_of_models_and_metrics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

#========== End of Ranking step ============================#

# Step. 4. Extract and finalize workflows for the top 2 models based on accuracy

# Identify the top 3 models based on AUC values
top_models <- auc_values %>%
  arrange(desc(.estimate)) %>%  
  slice_head(n = 3)  

print(top_models)

# wflow_id names in ML_results differs with top models, hence the mapping
top_models <- top_models %>%
  mutate(wflow_id = case_when(
    wflow_id == "rf" ~ "RF",        
    wflow_id == "svm" ~ "SVM",       
    TRUE ~ wflow_id  
  )) %>%
  mutate(mapped_wflow_id = case_when(
    wflow_id == "ENLR" ~ "recipe_glmnet",
    wflow_id == "RF" ~ "recipe_rf",    
    wflow_id == "SVM" ~ "recipe_svm"   
  ))

# Create empty list to store the finalized workflows
final_workflows <- vector("list", length(top_models$mapped_wflow_id))

# Iterate over each top model and extract the corresponding workflow
for (i in seq_along(top_models$mapped_wflow_id)) {
  wflow_id <- top_models$mapped_wflow_id[i]

  # Extract the best hyperparameters for the current workflow
  best_hyperparam <- ML_results %>% 
    extract_workflow_set_result(wflow_id) %>%
    select_best(metric = "roc_auc")

  # Extract the workflow object from ML_results
  workflow <- ML_results %>%
    extract_workflow(id = wflow_id)  
  
  # Finalize the workflow with the best hyperparameters
  finalized_workflow <- finalize_workflow(workflow, best_hyperparam)
  final_workflows[[i]] <- finalized_workflow
}

names(final_workflows) <- top_models$wflow_id
print(final_workflows)

#=============================================# 
# To avoid having to re-run the entire script #
#=============================================#

# Define the output directory for Final_workflows RDS Files
final_workflows_dir <- file.path(output_dir, "Final_workflows")

if (!dir.exists(final_workflows_dir)) {
  dir.create(final_workflows_dir, recursive = TRUE)
}

# Save each finalized workflow as an RDS file
for (model_name in names(final_workflows)) {
  file_path <- file.path(final_workflows_dir, paste0(model_name, ".rds"))
  saveRDS(final_workflows[[model_name]], file_path)
  cat("Saved workflow for", model_name, "to", file_path, "\n")
}

# Reload the saved workflows
rds_files <- list.files(final_workflows_dir, pattern = "*.rds", full.names = TRUE)
reloaded_workflows <- list()

for (file_path in rds_files) {
  # Extract the model name from the filename (remove the directory and file extension)
  model_name <- tools::file_path_sans_ext(basename(file_path))
  reloaded_workflows[[model_name]] <- readRDS(file_path)
  cat("Loaded workflow for", model_name, "from", file_path, "\n")
}

print(reloaded_workflows)

# Check if final_workflows or reloaded_workflows exists and then fit the models to data
if (exists("final_workflows")) {
  cat("Using `final_workflows` to create fitted workflows...\n")
  workflows_to_fit <- final_workflows
} else if (exists("reloaded_workflows")) {
  cat("Using `reloaded_workflows` to create fitted workflows...\n")
  workflows_to_fit <- reloaded_workflows
} else {
  stop("Neither `final_workflows` nor `reloaded_workflows` exists. Please ensure that one is available.")
}

# Create an empty list to store the fitted workflows
fitted_workflows <- list()

# Loop through the workflows to fit each on the data 
# Make sure the Steatogenic label is set as reference -- relevel function performed earlier at line 146
for (model_name in names(workflows_to_fit)) {
  cat("Fitting workflow for model:", model_name, "\n")
  fitted_workflow <- workflows_to_fit[[model_name]] %>%
    fit(data = T_Data)
  fitted_workflows[[model_name]] <- fitted_workflow
}

# Save the fitted workflows/models in a folder for reloading later
# Define the output directory for the fit models
fitted_workflows_dir <- file.path(output_dir, "Fitted_workflows")

if (!dir.exists(fitted_workflows_dir)) {
  dir.create(fitted_workflows_dir, recursive = TRUE)
  cat("Created directory for fitted workflows:", fitted_workflows_dir, "\n")
}

# Save each fitted workflow as an RDS file
for (model_name in names(fitted_workflows)) {
  file_path <- file.path(fitted_workflows_dir, paste0(model_name, "_fitted.rds"))
  saveRDS(fitted_workflows[[model_name]], file_path)
  cat("Saved fitted workflow for", model_name, "to", file_path, "\n")
}

# Check if fitted_workflows exists and is not empty
if (exists("fitted_workflows") && length(fitted_workflows) > 0) {
  cat("`fitted_workflows` already exists and contains", length(fitted_workflows), "workflows. Skipping reload.\n")
} else {
  cat("`fitted_workflows` does not exist or is empty. Reloading workflows from the fitted_workflows directory...\n")
  
  # Initialize the fitted_workflows list
  fitted_workflows <- list()
  
  # Get the list of RDS files from the fitted_workflows_dir
  rds_files <- list.files(fitted_workflows_dir, pattern = "\\_fitted\\.rds$", full.names = TRUE)
  
  # Reload each fitted workflow from the RDS files
  for (file_path in rds_files) {
    # Extract the model name from the filename (remove the directory and file extension)
    model_name <- tools::file_path_sans_ext(basename(file_path))
    fitted_workflows[[model_name]] <- readRDS(file_path)
    cat("Loaded fitted workflow for", model_name, "from", file_path, "\n")
  }
  
  # Rename the workflows to drop "_fitted" from the names
  names(fitted_workflows) <- gsub("_fitted$", "", names(fitted_workflows))
  cat("Cleaned workflow names after reload:", paste(names(fitted_workflows), collapse = ", "), "\n")
}

# Print the final list of workflow names
cat("Final workflow names:", paste(names(fitted_workflows), collapse = ", "), "\n")


#===============================================#
# Extract Important Features from Random Forest #
#===============================================#
# Extract the workflow for the random forest model
rf_workflow <- fitted_workflows[["RF"]]

# Extract fitted model
rf_fit <- extract_fit_parsnip(rf_workflow)

# Extract feature importance
rf_importance <- vip::vi(rf_fit)
print(rf_importance)

# Sort the importance scores in descending order and select the top 1000 important genes
RF_top_1000_genes <- rf_importance %>%
  arrange(desc(Importance)) %>%
  head(1000)

write.table(RF_top_1000_genes, file.path(ml_results_dir, "RF_top_1000_importance_score_genes.tsv"), sep = "\t", row.names = FALSE)
write.table(rf_importance, file.path(ml_results_dir, "RF_importance_scores.tsv"), sep = "\t", row.names = FALSE)

# Plot RF feature importance
RF_vip_plot <- vip::vip(rf_fit)
ggsave(filename = file.path(ml_results_dir, "RF_ViP_plot.pdf"), plot = RF_vip_plot, width = 8, height = 6, dpi = 1000)

#==================================================#
# Important features for SVM have -ve & +ve values #
# Hence only the absolute values are considered    #
#==================================================# 

# Extract SVM workflow only from the fitted_workflows (Workflows contains everything including recipe)
svm_workflow <- fitted_workflows[["SVM"]]

# Extract the SVM model itself (This contains only the SVM model)
svm_fit <- extract_fit_parsnip(svm_workflow)$fit

# Check & ensure it is a linear SVM model
if (inherits(svm_fit, "ksvm") && class(svm_fit@kernelf) == "vanillakernel") {
  # Extract the coefficients (only for linear SVM)
  coefficients <- t(svm_fit@coef[[1]]) %*% svm_fit@xmatrix[[1]]
  # Convert the coefficients matrix to a numeric vector for importance calculation
  importance <- abs(as.vector(coefficients))
  # Normalize importance scores
  importance <- importance / max(importance)
  # Create a data frame with variable importance
  importance_df <- data.frame(
    Variable = colnames(T_Data)[-ncol(T_Data)],  # Exclude the target column
    Importance = importance
  )
  # Order by importance
  importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
  
  # Print the importance scores
  print(importance_df)
  
  # Save the importance scores to a TSV file
  output_file <- file.path(ml_results_dir, "SVM_importance_scores.tsv")
  write.table(importance_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved importance scores to", output_file, "\n")
} else {
  message("This method of extracting feature importance is only applicable to linear SVM models.")
}

# Sort the importance scores in descending order and select the top 1000 important genes
SVM_top_1000_genes <- importance_df %>%
  arrange(desc(Importance)) %>%
  head(1000)

write.table(SVM_top_1000_genes, file.path(ml_results_dir, "SVM_top_1000_importance_score_genes.tsv"), sep = "\t", row.names = FALSE)

#===== End of feature extraction based on importance scores =====#

#=====================================#
# Recursive Feature Elimination (RFE) #
#=====================================#
# Recursive Feature Elimination for SVM
library(caret)

# Define control parameters for RFE
svm_control <- rfeControl(functions = caretFuncs, 
  method = "boot",         
  number = 1000)           

# SVM Model Specifications
svm_model <- trainControl(method = "boot", number = 1000)  

# Perform RFE
rfe_results <- rfe(
x = T_Data[, -ncol(T_Data)],  # Features
y = T_Data$Steatogenic,       # Target variable
sizes = 1001, # Number of features to test
rfeControl = svm_control,
method = "svmLinear",
trControl = svm_model
)

print(rfe_results)


# RFE for Random forest
#registerDoParallel(cores = detectCores() - 2)
# Step 1: Define the control function for RFE with bootstrapping
rf_control <- rfeControl(functions = rfFuncs, method = "boot", number = 1000)  

# Step 2: Run RFE to select important features
set.seed(1234) 
results <- rfe(T_Data, T_Data$Steatogenic, sizes = 1001, rfeControl = rf_control)
write.table(results$optVariables, file.path(ml_results_dir, "Random_Forest_RFE_genes.tsv"), sep = "\t", row.names = FALSE)

# Get the selected features
selected_features <- predictors(results)

# Step 3: Define the recipe with the selected features
rfe_recipe <- recipe(Steatogenic ~ ., data = T_Data) %>%
  step_nzv(all_predictors()) %>%
  step_smote(Steatogenic) %>%
  step_select(all_of(selected_features))

rfe_recipe_prep <- rfe_recipe %>% 
  prep() %>% 
  juice()

print(head(rfe_recipe_prep))

# Define the RFE model specification
rfe_model <- rand_forest() %>%
  set_args(mtry = floor(sqrt(length(selected_features)))) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

rfe_model

# Step 4: Create the workflow
rf_workflow <- workflow() %>%
  add_recipe(rfe_recipe) %>%
  add_model(rfe_model)

# Step 5: Fit the final model to the entire data
final_rfe_fit <- fit(rf_workflow, data = T_Data)

#=================== End of RFE =======================#

#============================================================================================#
# Perform predictions using the fitted models on entirely new dataset e.g Full TG-GATES Data #
#============================================================================================#

# Define the output directory for the prediction plots
predictions_dir <- file.path(output_dir, "Predictions")

if (!dir.exists(predictions_dir)) {
  dir.create(predictions_dir, recursive = TRUE)
}

# Import the full TG-GATES Liver Data 
Full_TG_GATES <- read.delim("Normalised_TG_GATES_Full_Human.tsv", row.names = "Gene", sep = "\t", header = TRUE, check.names = FALSE)
TG_GATES_metadata <- read.delim("cleaned_Full_TG_Human_metadata.tsv", sep = "\t", header = TRUE, row.names = "BARCODE")

# Order the column & row names in full TG GATES data & metadata
Full_TG_GATES <- Full_TG_GATES[,order(colnames(Full_TG_GATES))]
TG_GATES_metadata <- TG_GATES_metadata[order(rownames(TG_GATES_metadata)),]
colnames(Full_TG_GATES) <- TG_GATES_metadata$Suggested_name

# Transpose the data
t_Full_TG_GATES <- as.data.frame(x = t(Full_TG_GATES), stringsAsFactors = FALSE) 

Predictions_on_full_data <- lapply(fitted_workflows, function(workflow) {
  predict(workflow, new_data = t_Full_TG_GATES, type = "prob")
})

# Check the prediction probabilities made by each model on the full TG-GATES Dataset
# Function to filter data and extract prediction probabilities
filter_data_and_predictions <- function(label, t_Full_TG_GATES, T_Data, Predictions_on_full_data) {
  samples <- rownames(T_Data[T_Data$Steatogenic == label, ])
  filtered_data <- t_Full_TG_GATES[rownames(t_Full_TG_GATES) %in% samples, ]
  indices <- match(rownames(filtered_data), rownames(t_Full_TG_GATES))
  
  if (label == "Yes") {
    filtered_data$Pred_SVM <- Predictions_on_full_data[[1]]$.pred_Yes[indices]
    filtered_data$Pred_RF <- Predictions_on_full_data[[2]]$.pred_Yes[indices]
    filtered_data$Pred_ENLR <- Predictions_on_full_data[[3]]$.pred_Yes[indices]
  } else {
    filtered_data$Pred_SVM <- Predictions_on_full_data[[1]]$.pred_No[indices]
    filtered_data$Pred_RF <- Predictions_on_full_data[[2]]$.pred_No[indices]
    filtered_data$Pred_ENLR <- Predictions_on_full_data[[3]]$.pred_No[indices]
  }
  
  return(filtered_data)
}

# Function to generate and save bar plots for the models
generate_barplot <- function(filtered_data, pred_col, model_name, label, output_dir) {
  filtered_data$cleaned_names <- gsub("TG_hum_vit_", "", rownames(filtered_data))
  
  plot_title <- ifelse(label == "Yes", "Steatogenic", "Non-steatogenic")
  
  plot <- ggplot(filtered_data, aes(x = reorder(cleaned_names, !!sym(pred_col)), y = !!sym(pred_col))) +
    geom_bar(stat = "identity", fill = ifelse(model_name == "SVM", "red", ifelse(model_name == "RF", "blue", "darkgreen"))) +
    labs(title = paste(plot_title, "Prediction Probabilities -", model_name),
         x = "Samples",
         y = "Prediction Probability") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
    )
  
  # Save the plot to a file
  file_name <- paste0(plot_title, "_Prediction_Probabilities_", model_name, ".pdf")
  ggsave(filename = file.path(predictions_dir, file_name), plot = plot, width = 10, height = 6, dpi = 1000)
}

# Apply the workflow for both "Yes" and "No" samples
labels <- c("Yes", "No")
model_names <- c("SVM", "RF", "ENLR")
prediction_cols <- c("Pred_SVM", "Pred_RF", "Pred_ENLR")

for (label in labels) {
  # Filter data and predictions based on the label
  filtered_data <- filter_data_and_predictions(label, t_Full_TG_GATES, T_Data, Predictions_on_full_data)
  
  # Generate and save bar plots for each model
  for (i in seq_along(model_names)) {
    generate_barplot(filtered_data, prediction_cols[i], model_names[i], label, output_dir)
  }
}

# Save predictions on the Full TG-GATES Data for each model in a tsv file
sample_names <- rownames(t_Full_TG_GATES)

# Ensure that sample_names length matches the prediction data
if (length(sample_names) != nrow(Predictions_on_full_data[[1]])) {
  stop("Mismatch between number of sample names and predictions.")
}

export_predictions <- function(predictions, model_name, sample_names, output_dir) {
  # Add sample names to the predictions
  predictions_with_samples <- as.data.frame(predictions)
  predictions_with_samples$Sample <- sample_names
  
  # Reorder columns so that "Sample" is the first column
  predictions_with_samples <- predictions_with_samples[, c("Sample", colnames(predictions))]
  file_path <- file.path(predictions_dir, paste0(model_name, "_predictions.tsv"))
  write.table(predictions_with_samples, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Exported predictions for", model_name, "to", file_path, "\n")
}

# Export predictions for each model in Predictions_on_full_data
for (model_name in names(Predictions_on_full_data)) {
  export_predictions(Predictions_on_full_data[[model_name]], model_name, sample_names, predictions_dir)
}

#=========== End of Prediction saving =====================#

#============================================================================================#
# Set a prediction threshold to define what prediction prabability is considered steatogenic #
#============================================================================================#

# Function to apply threshold and classify
apply_threshold <- function(df, model_name, threshold) {
  df <- df %>%
    mutate(classification = case_when(
      .pred_Yes > threshold ~ "Steatogenic",
      .pred_No > threshold ~ "Non-steatogenic",
      TRUE ~ "Uncertain"
    )) %>%
    mutate(model = model_name)
  return(df)
}

# Set the custom prediction probability threshold
threshold <- 0.7
  
prediction_files <- list.files(path = predictions_dir, pattern = "*predictions.tsv", full.names = TRUE)

# Function to extract model name from file name 
extract_model_name <- function(filename) {
  model_name <- tools::file_path_sans_ext(basename(filename))
  model_name <- sub("_predictions$", "", model_name)  
  return(model_name)
} 

# Read all the files, apply threshold, and combine into one dataframe
combined_df <- do.call(rbind, lapply(prediction_files, function(file) {
  df <- read.csv(file, sep = "\t")
  model_name <- extract_model_name(file)
  apply_threshold(df, model_name, threshold)
}))

# Set the desired order of models (SVM, RF, ENLR)
combined_df <- combined_df %>%
  mutate(model = factor(model, levels = c("SVM", "RF", "ENLR")))

# Summarize the classifications by model
classification_summary <- combined_df %>%
  group_by(model, classification) %>%
  summarise(count = n()) %>%
  ungroup()

print(classification_summary)

# Create a stacked bar plot to visualize the results
prediction_summary <- ggplot(classification_summary, aes(x = model, y = count, fill = classification)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "white", size = 5) +  
  theme_minimal() +
  labs(title = "Comparison of prediction probabilities across models (Threshold 0.7)",
       x = "Model",
       y = "Number of Samples",
       fill = "Classification") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),      
        axis.title.x = element_text(size = 16),                 
        axis.title.y = element_text(size = 16),                 
        axis.text.x = element_text(size = 14),                  
        axis.text.y = element_text(size = 14),                  
        legend.text = element_text(size = 14),                  
        legend.title = element_text(size = 16),                 
        legend.key.size = unit(1.5, 'cm')) +                     
  scale_fill_manual(values = c("Steatogenic" = "orange", "Non-steatogenic" = "darkblue", "Uncertain" = "gray"))

print(prediction_summary)
ggsave(filename = file.path(predictions_dir, "Prediction_summary_Barplot.pdf"), plot = prediction_summary, width = 10, height = 6, dpi = 1000)

#==========================================================================#
# How does each model predict control samples in the Full TG-GATES dataset #
#==========================================================================#

# Filter combined_df to only include control (Ctl) samples
ctl_df <- combined_df %>%
  filter(str_detect(Sample, "Ctl"))

# Boxplot for Steatogenic prediction probabilities for control (Ctl) samples
ctl_boxplot_yes <- ggplot(ctl_df, aes(x = model, y = .pred_Yes, fill = model)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of control samples predicted as steatogenic",
       x = NULL,  
       y = "Prediction probability (Steatogenic)",
       fill = "Model") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 14),
    axis.ticks.x = element_blank(),  
    legend.position = "none"  
  ) +
  scale_fill_manual(values = c("SVM" = "maroon", "RF" = "blue", "ENLR" = "darkgreen"))

# Boxplot for Non-steatogenic prediction probabilities for Ctl samples
ctl_boxplot_no <- ggplot(ctl_df, aes(x = model, y = .pred_No, fill = model)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of control samples predicted as non-steatogenic",
       x = "Model",  
       y = "Prediction probability (Non-steatogenic)",
       fill = "Model") +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14),
    legend.position = "none"  
  ) +
  scale_fill_manual(values = c("SVM" = "maroon", "RF" = "blue", "ENLR" = "darkgreen"))

# Use patchwork library to combine the plots vertically with shared x-axis label 
stacked_boxplot <- ctl_boxplot_yes / ctl_boxplot_no +
  plot_layout(heights = c(1, 1.2)) 

print(stacked_boxplot)
ggsave(filename = file.path(predictions_dir, "Boxplot_for_ctrl_prediction_Probability.pdf"), plot = stacked_boxplot, width = 10, height = 10, dpi = 1000)

#===== Prediction on RNA-seq Data =====#

#=====================================================#
# Are our top performing models platform independent? #
# Perform Predictions on RNA-seq data (TPM values)    #
#=====================================================#

hecatos_data <- read.delim("Hecatos_TPM_counts.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Transpose the Hecatos data
t_hecatos_data <- as.data.frame(t(hecatos_data), stringsAsFactors = FALSE)

# Make predictions on the filtered Hecatos data for all models in 'fitted_workflows'
Predictions_on_Hecatos <- lapply(fitted_workflows, function(workflow) {
  predict(workflow, new_data = t_hecatos_data, type = "prob")
})

# Define Hecatos_predictions directory for the plots
hecatos_predictions_dir <- file.path(predictions_dir, "Hecatos_predictions")

if (!dir.exists(hecatos_predictions_dir)) {
  dir.create(hecatos_predictions_dir, recursive = TRUE)
}

save_and_plot_predictions_separately <- function(predictions, model_name, output_dir) {
  # Add row names (sample names) to the predictions
  predictions_with_samples <- as.data.frame(predictions)
  predictions_with_samples$Sample <- rownames(t_hecatos_data)
  
  # Reorder columns so that "Sample" is the first column
  predictions_with_samples <- predictions_with_samples %>%
    select(Sample, everything())
  
  # Save predictions to a .tsv file
  output_file <- file.path(hecatos_predictions_dir, paste0(model_name, "_predictions.tsv"))
  write.table(predictions_with_samples, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved predictions for", model_name, "to", output_file, "\n")
  
  # Reshape predictions for plotting (gather Yes/No predictions)
  predictions_long <- predictions_with_samples %>%
    pivot_longer(cols = starts_with(".pred_"), names_to = "Prediction_Class", values_to = "Prediction_Probability")
  
  # Separate barplots for 'Yes' and 'No' predictions
  for (prediction_class in unique(predictions_long$Prediction_Class)) {
    # Filter the data for the current prediction class
    filtered_data <- predictions_long %>%
      filter(Prediction_Class == prediction_class)
    
    # Create the barplot for the current class
    plot <- ggplot(filtered_data, aes(x = reorder(Sample, Prediction_Probability), y = Prediction_Probability)) +
      geom_bar(stat = "identity", fill = ifelse(prediction_class == ".pred_Yes", "red", "blue")) +
      labs(title = paste("Prediction Probabilities for", ifelse(prediction_class == ".pred_Yes", "Yes", "No"), "-", model_name),
           x = "Samples",
           y = "Prediction Probability") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)
      )
    
    # Save the barplot as a PDF file
    plot_file <- file.path(hecatos_predictions_dir, paste0(model_name, "_", prediction_class, "_Prediction_Barplot.pdf"))
    ggsave(filename = plot_file, plot = plot, width = 10, height = 6, dpi = 1000)
    cat("Saved", prediction_class, "barplot for", model_name, "to", plot_file, "\n")
  }
}


# Loop over each model's predictions and generate separate barplots for 'Yes' and 'No'
model_names <- c("SVM", "RF", "ENLR")  
for (i in seq_along(Predictions_on_Hecatos)) {
  save_and_plot_predictions_separately(Predictions_on_Hecatos[[i]], model_names[i], output_dir)
}






































