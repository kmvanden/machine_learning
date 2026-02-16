# Random Forest Workflow for Metabolomics Data

# load libraries
library(ggplot2)
library(cowplot)
library(patchwork)
library(tidyverse)
library(readxl)
library(impute)
library(randomForest)
library(caret)
library(pROC)
library(Boruta)
library(ParBayesianOptimization)
library(doParallel)
library(fastshap)
library(viridis)
library(ggrepel)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
### metadata
meta <- read.table("metadata.txt", header = TRUE)

# subset meta to samples present in metabolomics data
sample_key <- read_excel("metabo_batch.xlsx", sheet = "Sample Meta Data") # sample name key
sample_key <- sample_key %>% mutate(CLIENT_SAMPLE_ID = sub("KMV_", "SMS_", CLIENT_SAMPLE_ID)) # rename sample ids
all(sample_key$CLIENT_SAMPLE_ID %in% meta$sample_id)
meta <- meta %>%
  left_join(sample_key %>% dplyr::select(CLIENT_SAMPLE_ID, PARENT_SAMPLE_NAME),
            by = c("sample_id" = "CLIENT_SAMPLE_ID")) %>% # add PARENT_SAMPLE_NAME to meta
  filter(!is.na(PARENT_SAMPLE_NAME)) # subset

# convert condition column into a factor
rownames(meta) <- meta$sample_name
meta$condition <- as.factor(meta$condition) 
table(meta$condition)


### feature table
feat <- read_excel("metabo_batch.xlsx", sheet = "Batch-normalized Data") # batch-normalized feature table
all(feat$PARENT_SAMPLE_NAME %in% meta$PARENT_SAMPLE_NAME)

# replace names in PARENT_SAMPLE_NAME column in feat with names from sample_name (from meta)
feat <- feat %>%
  left_join(meta %>% dplyr::select(PARENT_SAMPLE_NAME, sample_name), by = "PARENT_SAMPLE_NAME") %>%
  mutate(PARENT_SAMPLE_NAME = sample_name) %>%
  dplyr::select(-sample_name) %>%
  slice(match(meta$sample_name, PARENT_SAMPLE_NAME)) %>% # match order of samples in meta
  column_to_rownames(var = "PARENT_SAMPLE_NAME") # move sample id column to rownames
all(rownames(feat) == rownames(meta))

# replace CHEM_ID (colnames of feat) with names from CHEMICAL_NAME (from feat_key)
feat_key <- read_excel("metabo_batch.xlsx", sheet = "Chemical Annotation")
all(colnames(feat) == feat_key$CHEM_ID)
name_map <- feat_key %>% dplyr::select(CHEM_ID, CHEMICAL_NAME) %>% deframe() 
colnames(feat) <- name_map[colnames(feat)] %>% ifelse(is.na(.), colnames(feat), .)
all(colnames(feat) == feat_key$CHEMICAL_NAME)
colnames(feat) <- make.names(colnames(feat)) # make sure names are syntactically valid 


### feature filtering, zero imputation and log transformation 
# assess missingness (zeros listed as NAs in data.frame)
total_missing <- sum(is.na(feat)) / prod(dim(feat)) * 100 # overall missingness
missing_per_feature <- colMeans(is.na(feat)) * 100 # missingness per feature
summary(missing_per_feature)
missing_per_sample <- rowMeans(is.na(feat)) * 100 # missingness per sample
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
feat_filtered <- feat[, missing_per_feature <= 30]

# impute remaining zeros with kNN
feat_imputed <- impute.knn(as.matrix(t(feat_filtered)))$data # impute.knn expects features as rows
feat_imputed <- t(feat_imputed)

# log2 transformation (with pseudocount)
feat_log <- log2(feat_imputed + 1e-6)
feat_log <- as.data.frame(feat_log)

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat_log))
feat_log$sample_name <- rownames(feat_log) # add column sample_name
meta <- meta[, -c(3,4)] # remove meta columns that were used to format data.frames


### merge metadata and feature table
metabo <- merge(meta, feat_log, by = "sample_name", all.x = TRUE)
metabo <- metabo[,-1] # remove sample_name


###############################################################################
###   BASELINE RANDOM FOREST MODEL - 5-FOLD CROSS-VALIDATION + 50 REPEATS   ###
###############################################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metabo), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  set.seed(1234 + r*100)
  folds <- createFolds(metabo$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
    test_data <- metabo[test_idx, ] # testing data (fold f)
    
    # train random forest model
    rf_model <- randomForest(x = train_data[, all_feat_cols], 
                             y = as.factor(train_data$condition), 
                             ntree = 500, 
                             importance = TRUE) 
    
    # evaluate on test set
    test_predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
    test_probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
    # evaluate model on training set
    train_predictions <- predict(rf_model, newdata = train_data[, all_feat_cols], type = "response")
    train_probabilities <- predict(rf_model, newdata = train_data[, all_feat_cols], type = "prob")
    
    # calculate AUC on test set
    test_roc_obj <- roc(response = test_data$condition,
                        predictor = test_probabilities[, "disease"],
                        levels = c("healthy", "disease"),
                        direction = "<")
    test_auc <- auc(test_roc_obj)
    
    # store test ROC coordinates
    test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                              sensitivity = test_roc_obj$sensitivities,
                              Repeat = r, Fold = f, Set = "Test")
    
    # calculate AUC on train set
    train_roc_obj <- roc(response = train_data$condition,
                         predictor = train_probabilities[, "disease"],
                         levels = c("healthy", "disease"),
                         direction = "<")
    train_auc <- auc(train_roc_obj)
    
    # store train ROC coordinates
    train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                               sensitivity = train_roc_obj$sensitivities,
                               Repeat = r, Fold = f, Set = "Train")
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrices
    test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
    train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
    
    
    ### store with repeat (r) and fold (f) index
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                       train_cm = train_cm, train_auc = train_auc,
                                       test_roc_df = test_roc_df, train_roc_df = train_roc_df) # store performance metrics (test and train)
    feature_importances[[key]] <- importance(rf_model)  # store feature importances
  }
}


### calculate feature frequencies and relative frequency of feature selection
feature_split_summary <- bind_rows(feature_frequencies, .id = "Repeat_Fold") %>%
  rename(Feature = tree_split_vars) %>%
  group_by(Feature) %>% # group stability metrics by feature
  summarise(total_count = sum(Freq, na.rm = TRUE),
            mean_count = mean(Freq, na.rm = TRUE),
            n_models = n()) %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree)) %>%
  arrange(desc(total_count))
head(as.data.frame(feature_split_summary), 20)

# plot total number of models where feature was used at least once
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - models",
       x = "Feature", y = "Number of models")

# plot average number of times feature was used in a split per tree (across all models) 
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - splits",
       x = "Feature", y = "Average splits per tree")


### calculate feature importances
# combine all feature_importances data.frames into one data.frame
all_features_importances <- do.call(rbind, lapply(names(feature_importances), function(name) {
  df <- as.data.frame(feature_importances[[name]])
  df$Feature <- rownames(df)
  df$Repeat_Fold <- name
  return(df)
}))

# group importance metrics by feature and sort by overall importance
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(mean_importance, 20)


### plot metabolite with highest MeanDecreaseAccuracy
ggplot(metabo, aes(x = beta.alanine)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative metabolite",
       subtitle = "Beta alanine",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()


### plot split frequency versus feature importance
# top 20 features by split frequency
top_features <- feature_split_summary %>% 
  arrange(desc(avg_per_tree)) %>% 
  slice(1:20) # top 20 features by split frequency

# merge with mean_importance data.frame
top_features_importance <- top_features %>%
  left_join(mean_importance, by = c("Feature"))

# plot meanMDA versus split frequency
ggplot(top_features_importance, aes(x = avg_per_tree, y = mean_MeanDecreaseAccuracy, label = Feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA",
       title = "Split frequency versus permutation importance")

# long data.frame for plotting importance by condition
top_features_long <- top_features_importance %>%
  select(Feature, avg_per_tree, mean_healthy, mean_disease) %>%
  pivot_longer(cols = c(mean_healthy, mean_disease),
               names_to = "Condition",
               values_to = "Importance") %>%
  mutate(Condition = recode(Condition,
                            mean_healthy = "Healthy",
                            mean_disease = "Disease"))

# plot meanMDA by condition versus split frequency
ggplot(top_features_long, aes(x = avg_per_tree, y = Importance, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA by condition",
       title = "Split frequency versus permutation importance by condition", color = "Condition") +
  scale_color_manual(values = c("Healthy" = "steelblue", "Disease" = "indianred3"))


### calculate performance statistics
perf_stats <- function(performance_metrics, type = c("test", "train", "gap")) {
  
  # match type argument
  type <- match.arg(type)
  
  # create vectors to store metrics - test 
  test_balanced_accuracy <- numeric()
  test_f1_score <- numeric()
  test_sensitivity <- numeric()
  test_specificity <- numeric()
  test_auc_vals <- numeric()
  
  # create vectors to store metrics - train 
  train_balanced_accuracy <- numeric()
  train_f1_score <- numeric()
  train_sensitivity <- numeric()
  train_specificity <- numeric()
  train_auc_vals <- numeric()
  
  # extract metrics from the stored confusion matrices
  for (perf in performance_metrics) {
    test_cm <- perf$test_cm
    test_auc_val <- as.numeric(perf$test_auc[])
    train_cm <- perf$train_cm
    train_auc_val <- as.numeric(perf$train_auc[])
    
    
    # confusion matrix metrics and auc (test)
    test_balanced_accuracy <- c(test_balanced_accuracy, test_cm$byClass["Balanced Accuracy"])
    test_f1_score <- c(test_f1_score, test_cm$byClass["F1"])
    test_sensitivity <- c(test_sensitivity, test_cm$byClass["Sensitivity"])
    test_specificity <- c(test_specificity, test_cm$byClass["Specificity"])
    test_auc_vals <- c(test_auc_vals, test_auc_val)
    
    # confusion matrix metrics and auc (train)
    train_balanced_accuracy <- c(train_balanced_accuracy, train_cm$byClass["Balanced Accuracy"])
    train_f1_score <- c(train_f1_score, train_cm$byClass["F1"])
    train_sensitivity <- c(train_sensitivity, train_cm$byClass["Sensitivity"])
    train_specificity <- c(train_specificity, train_cm$byClass["Specificity"])
    train_auc_vals <- c(train_auc_vals, train_auc_val)
  }
  
  # test metric summary
  test_metric_summary <- data.frame(mean_bal_acc = mean(test_balanced_accuracy, na.rm = TRUE),
                                    sd_bal_acc = sd(test_balanced_accuracy, na.rm = TRUE),
                                    mean_f1 = mean(test_f1_score, na.rm = TRUE),
                                    sd_f1 = sd(test_f1_score, na.rm = TRUE),
                                    mean_sens = mean(test_sensitivity, na.rm = TRUE),
                                    sd_sens = sd(test_sensitivity, na.rm = TRUE),
                                    mean_spec = mean(test_specificity, na.rm = TRUE),
                                    sd_spec = sd(test_specificity, na.rm = TRUE),
                                    mean_auc = mean(test_auc_vals, na.rm = TRUE),
                                    sd_auc = sd(test_auc_vals, na.rm = TRUE))
  
  # train metric summary
  train_metric_summary <- data.frame(mean_bal_acc = mean(train_balanced_accuracy, na.rm = TRUE),
                                     sd_bal_acc = sd(train_balanced_accuracy, na.rm = TRUE),
                                     mean_f1 = mean(train_f1_score, na.rm = TRUE),
                                     sd_f1 = sd(train_f1_score, na.rm = TRUE),
                                     mean_sens = mean(train_sensitivity, na.rm = TRUE),
                                     sd_sens = sd(train_sensitivity, na.rm = TRUE),
                                     mean_spec = mean(train_specificity, na.rm = TRUE),
                                     sd_spec = sd(train_specificity, na.rm = TRUE),
                                     mean_auc = mean(train_auc_vals, na.rm = TRUE),
                                     sd_auc = sd(train_auc_vals, na.rm = TRUE))
  
  # gap metric summary
  gap_metric_summary <- data.frame(mean_bal_acc = mean(train_balanced_accuracy - test_balanced_accuracy, na.rm = TRUE),
                                   sd_bal_acc = sd(train_balanced_accuracy - test_balanced_accuracy, na.rm = TRUE),
                                   mean_f1 = mean(train_f1_score - test_f1_score, na.rm = TRUE),
                                   sd_f1 = sd(train_f1_score - test_f1_score, na.rm = TRUE),
                                   mean_sens = mean(train_sensitivity - test_sensitivity, na.rm = TRUE),
                                   sd_sens = sd(train_sensitivity - test_sensitivity, na.rm = TRUE),
                                   mean_spec = mean(train_specificity - test_specificity, na.rm = TRUE),
                                   sd_spec = sd(train_specificity - test_specificity, na.rm = TRUE),
                                   mean_auc = mean(train_auc_vals - test_auc_vals, na.rm = TRUE),
                                   sd_auc = sd(train_auc_vals - test_auc_vals, na.rm = TRUE))
  
  # summary to return
  result <- switch(type, "test" = test_metric_summary,
                   "train" = train_metric_summary,
                   "gap" = gap_metric_summary)
  
  return(result)
}

perf_stats(performance_metrics, type = "test")
perf_stats(performance_metrics, type = "train")
perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
plot_roc <- function(performance_metrics) {
  
  # combine train and test ROC data frames
  all_roc_curves <- bind_rows(lapply(performance_metrics, function(x) bind_rows(x$train_roc_df, x$test_roc_df)))
  
  # compute average ROC curves for train and test set
  fpr_grid <- seq(0, 1, length.out = 100)
  
  interp_roc <- all_roc_curves %>%
    group_by(Set, Repeat, Fold) %>%
    reframe(tpr_interp = approx(1 - specificity, sensitivity, xout = fpr_grid, ties = mean)$y,
            .groups = "drop") %>%
    mutate(fpr = rep(fpr_grid, times = n() / length(fpr_grid)))
  
  mean_roc <- interp_roc %>%
    group_by(Set, fpr) %>%
    summarise(mean_tpr = mean(tpr_interp, na.rm = TRUE),
              lower_tpr = quantile(tpr_interp, 0.025, na.rm = TRUE),
              upper_tpr = quantile(tpr_interp, 0.975, na.rm = TRUE),
              .groups = "drop")
  
  # plot train and test ROC curves
  p <- ggplot(mean_roc, aes(x = fpr, y = mean_tpr, color = Set, fill = Set)) +
    geom_line(linewidth = 1.2) + coord_equal() + theme_minimal() +
    geom_ribbon(aes(ymin = lower_tpr, ymax = upper_tpr), alpha = 0.2, color = NA) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Train" = "indianred3", "Test" = "steelblue")) +
    scale_fill_manual(values = c("Train" = "indianred3", "Test" = "steelblue")) +
    labs(title = "Average ROC curves across CV folds",
         x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)", 
         color = "Dataset", fill = "Dataset")
  
  return(p)
}

plot_roc(performance_metrics)


##################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - CLASSWT  ###
##################################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")

# create list of class weight settings
weight_grid <- list(healthy = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    disease = c(healthy = 1, disease = 2))

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop through class weight
for (i in names(weight_grid)) {
  classwt <- weight_grid[[i]]
  cat("Class weight:", i, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    set.seed(1234 + r*100)
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metabo$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
      test_data <- metabo[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, 
                               importance = TRUE, 
                               classwt = classwt)
      
      # evaluate on test set
      test_predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
      test_probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # evaluate model on training set
      train_predictions <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "response")
      train_probabilities <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "prob")
      
      # calculate AUC on test set
      test_roc_obj <- roc(response = test_data$condition,
                          predictor = test_probabilities[, "disease"],
                          levels = c("healthy", "disease"),
                          direction = "<")
      test_auc <- auc(test_roc_obj)
      
      # store test ROC coordinates
      test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                                sensitivity = test_roc_obj$sensitivities,
                                Repeat = r, Fold = f, Set = "Test")
      
      # calculate AUC on train set
      train_roc_obj <- roc(response = train_data$condition,
                           predictor = train_probabilities[, "disease"],
                           levels = c("healthy", "disease"),
                           direction = "<")
      train_auc <- auc(train_roc_obj)
      
      # store train ROC coordinates
      train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                                 sensitivity = train_roc_obj$sensitivities,
                                 Repeat = r, Fold = f, Set = "Train")
      
      # generate confusion matrices
      test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
      train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
      
      
      ### store with repeat (r) and fold (f) index
      key <- paste0("classwt", i, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                         train_cm = train_cm, train_auc = train_auc,
                                         test_roc_df = test_roc_df, train_roc_df = train_roc_df,
                                         param_value = i) # store performance metrics (test and train)
    }
  }
}

### calculate performance statistics
grid_perf_stats <- function(performance_metrics, type = c("test", "train", "gap")) {
  
  # match type argument
  type <- match.arg(type)
  
  # get unique hyperparameter values
  param_values <- unique(sapply(performance_metrics, function(x) x$param_value))
  
  # initialize list to store summaries
  summary_list <- list()
  
  for (i in param_values) {
    # filter metrics for this class weight
    perf_subset <- performance_metrics[sapply(performance_metrics, function(x) x$param_value == i)]
    
    # create vectors to store metrics - test 
    test_balanced_accuracy <- numeric()
    test_f1_score <- numeric()
    test_sensitivity <- numeric()
    test_specificity <- numeric()
    test_auc_vals <- numeric()
    
    # create vectors to store metrics - train 
    train_balanced_accuracy <- numeric()
    train_f1_score <- numeric()
    train_sensitivity <- numeric()
    train_specificity <- numeric()
    train_auc_vals <- numeric()
    
    # extract metrics from the stored confusion matrices
    for (perf in perf_subset) {
      test_cm <- perf$test_cm
      test_auc_val <- as.numeric(perf$test_auc[])
      train_cm <- perf$train_cm
      train_auc_val <- as.numeric(perf$train_auc[])
      
      # confusion matrix metrics and auc (test)
      test_balanced_accuracy <- c(test_balanced_accuracy, test_cm$byClass["Balanced Accuracy"])
      test_f1_score <- c(test_f1_score, test_cm$byClass["F1"])
      test_sensitivity <- c(test_sensitivity, test_cm$byClass["Sensitivity"])
      test_specificity <- c(test_specificity, test_cm$byClass["Specificity"])
      test_auc_vals <- c(test_auc_vals, test_auc_val)
      
      # confusion matrix metrics and auc (train)
      train_balanced_accuracy <- c(train_balanced_accuracy, train_cm$byClass["Balanced Accuracy"])
      train_f1_score <- c(train_f1_score, train_cm$byClass["F1"])
      train_sensitivity <- c(train_sensitivity, train_cm$byClass["Sensitivity"])
      train_specificity <- c(train_specificity, train_cm$byClass["Specificity"])
      train_auc_vals <- c(train_auc_vals, train_auc_val)
    }
    
    # compute summary per type
    if (type == "test") {
      df <- data.frame(param_value = i,
                       mean_bal_acc = mean(test_balanced_accuracy, na.rm = TRUE),
                       sd_bal_acc = sd(test_balanced_accuracy, na.rm = TRUE),
                       mean_f1 = mean(test_f1_score, na.rm = TRUE),
                       sd_f1 = sd(test_f1_score, na.rm = TRUE),
                       mean_sens = mean(test_sensitivity, na.rm = TRUE),
                       sd_sens = sd(test_sensitivity, na.rm = TRUE),
                       mean_spec = mean(test_specificity, na.rm = TRUE),
                       sd_spec = sd(test_specificity, na.rm = TRUE),
                       mean_auc = mean(test_auc_vals, na.rm = TRUE),
                       sd_auc = sd(test_auc_vals, na.rm = TRUE))
      
    } else if (type == "train") {
      df <- data.frame(param_value = i,
                       mean_bal_acc = mean(train_balanced_accuracy, na.rm = TRUE),
                       sd_bal_acc = sd(train_balanced_accuracy, na.rm = TRUE),
                       mean_f1 = mean(train_f1_score, na.rm = TRUE),
                       sd_f1 = sd(train_f1_score, na.rm = TRUE),
                       mean_sens = mean(train_sensitivity, na.rm = TRUE),
                       sd_sens = sd(train_sensitivity, na.rm = TRUE),
                       mean_spec = mean(train_specificity, na.rm = TRUE),
                       sd_spec = sd(train_specificity, na.rm = TRUE),
                       mean_auc = mean(train_auc_vals, na.rm = TRUE),
                       sd_auc = sd(train_auc_vals, na.rm = TRUE))
      
    } else if (type == "gap") {
      df <- data.frame(param_value = i,
                       mean_bal_acc = mean(train_balanced_accuracy - test_balanced_accuracy, na.rm = TRUE),
                       sd_bal_acc = sd(train_balanced_accuracy - test_balanced_accuracy, na.rm = TRUE),
                       mean_f1 = mean(train_f1_score - test_f1_score, na.rm = TRUE),
                       sd_f1 = sd(train_f1_score - test_f1_score, na.rm = TRUE),
                       mean_sens = mean(train_sensitivity - test_sensitivity, na.rm = TRUE),
                       sd_sens = sd(train_sensitivity - test_sensitivity, na.rm = TRUE),
                       mean_spec = mean(train_specificity - test_specificity, na.rm = TRUE),
                       sd_spec = sd(train_specificity - test_specificity, na.rm = TRUE),
                       mean_auc = mean(train_auc_vals - test_auc_vals, na.rm = TRUE),
                       sd_auc = sd(train_auc_vals - test_auc_vals, na.rm = TRUE))
      
    }
    
    summary_list[[i]] <- df
  }
  
  # combine all param value summaries
  summary_df <- bind_rows(summary_list)
  return(summary_df)
}

grid_perf_stats(performance_metrics, type = "test")
grid_perf_stats(performance_metrics, type = "train")
grid_perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
grid_plot_roc <- function(performance_metrics) {
  
  # get unique hyperparameter values
  param_values <- unique(sapply(performance_metrics, function(x) x$param_value))
  
  # function to compute mean ROC per hyperparameter value
  mean_roc_per_param <- function(param_val) {
    
    # filter performance metrics by this hyperparameter value
    perf_subset <- performance_metrics[sapply(performance_metrics, function(x) x$param_value == param_val)]
    
    # combine train and test ROC data frames
    all_roc_curves <- bind_rows(lapply(perf_subset, function(x) bind_rows(x$train_roc_df, x$test_roc_df)))
    
    fpr_grid <- seq(0, 1, length.out = 100)
    
    interp_roc <- all_roc_curves %>%
      group_by(Set, Repeat, Fold) %>%
      reframe(tpr_interp = approx(1 - specificity, sensitivity, xout = fpr_grid, ties = mean)$y,
              .groups = "drop") %>%
      mutate(fpr = rep(fpr_grid, times = n() / length(fpr_grid)))
    
    # compute mean and 95% CI
    mean_roc <- interp_roc %>%
      group_by(Set, fpr) %>%
      summarise(mean_tpr = mean(tpr_interp, na.rm = TRUE),
                lower_tpr = quantile(tpr_interp, 0.025, na.rm = TRUE),
                upper_tpr = quantile(tpr_interp, 0.975, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(param_value = param_val)
    
    return(mean_roc)
  }
  
  # compute mean ROC for all hyperparameter values
  roc_list <- lapply(param_values, mean_roc_per_param)
  roc_df <- bind_rows(roc_list)
  
  # plot train and test ROC curves
  p <- ggplot(roc_df, aes(x = fpr, y = mean_tpr, color = Set, fill = Set)) +
    geom_line(linewidth = 0.5) + coord_equal() + theme_minimal() +
    facet_wrap(~ param_value, ncol = 2, nrow = 3) + 
    geom_ribbon(aes(ymin = lower_tpr, ymax = upper_tpr), alpha = 0.2, color = NA) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Train" = "indianred3", "Test" = "steelblue")) +
    scale_fill_manual(values = c("Train" = "indianred3", "Test" = "steelblue")) +
    labs(title = "Average ROC curves across CV folds",
         x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)", 
         color = "Dataset", fill = "Dataset")
  
  return(p)
}

grid_plot_roc(performance_metrics)


################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - NTREE  ###
################################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")

# ntree values to test
ntree_values <- c(125, 250, 500, 1000, 2000)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop through ntree values
for (i in ntree_values) {
  cat("ntree:", i, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    set.seed(1234 + r*100)
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metabo$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
      test_data <- metabo[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = i, 
                               importance = TRUE)
      
      # evaluate on test set
      test_predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
      test_probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # evaluate model on training set
      train_predictions <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "response")
      train_probabilities <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "prob")
      
      # calculate AUC on test set
      test_roc_obj <- roc(response = test_data$condition,
                          predictor = test_probabilities[, "disease"],
                          levels = c("healthy", "disease"),
                          direction = "<")
      test_auc <- auc(test_roc_obj)
      
      # store test ROC coordinates
      test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                                sensitivity = test_roc_obj$sensitivities,
                                Repeat = r, Fold = f, Set = "Test")
      
      # calculate AUC on train set
      train_roc_obj <- roc(response = train_data$condition,
                           predictor = train_probabilities[, "disease"],
                           levels = c("healthy", "disease"),
                           direction = "<")
      train_auc <- auc(train_roc_obj)
      
      # store train ROC coordinates
      train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                                 sensitivity = train_roc_obj$sensitivities,
                                 Repeat = r, Fold = f, Set = "Train")
      
      # generate confusion matrices
      test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
      train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
      
      
      ### store with repeat (r) and fold (f) index
      key <- paste0("ntree_", i, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                         train_cm = train_cm, train_auc = train_auc,
                                         test_roc_df = test_roc_df, train_roc_df = train_roc_df,
                                         param_value = i) # store performance metrics (test and train)
    }
  }
}


### calculate performance statistics
grid_perf_stats(performance_metrics, type = "test")
grid_perf_stats(performance_metrics, type = "train")
grid_perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
grid_plot_roc(performance_metrics)


###############################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - MTRY  ###
###############################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")

# mtry values to test
mtry_values <- c(1, 2, 3, 4, 5, 6)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop through mtry values
for (i in mtry_values) {
  cat("mtry:", i, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    set.seed(1234 + r*100)
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metabo$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
      test_data <- metabo[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, 
                               importance = TRUE, 
                               mtry = i)
      
      # evaluate on test set
      test_predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
      test_probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # evaluate model on training set
      train_predictions <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "response")
      train_probabilities <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "prob")
      
      # calculate AUC on test set
      test_roc_obj <- roc(response = test_data$condition,
                          predictor = test_probabilities[, "disease"],
                          levels = c("healthy", "disease"),
                          direction = "<")
      test_auc <- auc(test_roc_obj)
      
      # store test ROC coordinates
      test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                                sensitivity = test_roc_obj$sensitivities,
                                Repeat = r, Fold = f, Set = "Test")
      
      # calculate AUC on train set
      train_roc_obj <- roc(response = train_data$condition,
                           predictor = train_probabilities[, "disease"],
                           levels = c("healthy", "disease"),
                           direction = "<")
      train_auc <- auc(train_roc_obj)
      
      # store train ROC coordinates
      train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                                 sensitivity = train_roc_obj$sensitivities,
                                 Repeat = r, Fold = f, Set = "Train")
      
      # generate confusion matrices
      test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
      train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
      
      
      ### store with repeat (r) and fold (f) index
      key <- paste0("mtry_", i, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                         train_cm = train_cm, train_auc = train_auc,
                                         test_roc_df = test_roc_df, train_roc_df = train_roc_df,
                                         param_value = i) # store performance metrics (test and train)
    }
  }
}

### calculate performance statistics
grid_perf_stats(performance_metrics, type = "test")
grid_perf_stats(performance_metrics, type = "train")
grid_perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
grid_plot_roc(performance_metrics)


###################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - NODESIZE  ###
###################################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")

# nodesize values to test
nodesize_values <- c(1, 2, 3, 4, 5, 6)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop through nodesize values
for (i in nodesize_values) {
  cat("nodesize:", i, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    set.seed(1234 + r*100)
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metabo$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
      test_data <- metabo[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, 
                               importance = TRUE, 
                               nodesize = i)
      
      # evaluate on test set
      test_predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
      test_probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # evaluate model on training set
      train_predictions <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "response")
      train_probabilities <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "prob")
      
      # calculate AUC on test set
      test_roc_obj <- roc(response = test_data$condition,
                          predictor = test_probabilities[, "disease"],
                          levels = c("healthy", "disease"),
                          direction = "<")
      test_auc <- auc(test_roc_obj)
      
      # store test ROC coordinates
      test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                                sensitivity = test_roc_obj$sensitivities,
                                Repeat = r, Fold = f, Set = "Test")
      
      # calculate AUC on train set
      train_roc_obj <- roc(response = train_data$condition,
                           predictor = train_probabilities[, "disease"],
                           levels = c("healthy", "disease"),
                           direction = "<")
      train_auc <- auc(train_roc_obj)
      
      # store train ROC coordinates
      train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                                 sensitivity = train_roc_obj$sensitivities,
                                 Repeat = r, Fold = f, Set = "Train")
      
      # generate confusion matrices
      test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
      train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
      
      
      ### store with repeat (r) and fold (f) index
      key <- paste0("nodesize_", i, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                         train_cm = train_cm, train_auc = train_auc,
                                         test_roc_df = test_roc_df, train_roc_df = train_roc_df,
                                         param_value = i) # store performance metrics (test and train)
    }
  }
}

### calculate performance statistics
grid_perf_stats(performance_metrics, type = "test")
grid_perf_stats(performance_metrics, type = "train")
grid_perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
grid_plot_roc(performance_metrics)


###########################################################################################################
########   RANDOM FOREST - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATION OF BAYES OPT  #######
###########################################################################################################

# data to be used in the model
str(metabo)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")

# create list of class weight settings
weight_grid <- list(healthy = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    disease = c(healthy = 1, disease = 2))

# define the set of categorical labels with numeric indices
label_keys <- c("healthy", "equal", "disease")

# scoring function
scoring_function <- function(mtry, ntree, nodesize, classwt_label) {
  
  # set seed
  set.seed(1234)
  
  # parameters
  mtry <- as.integer(mtry)
  ntree <- as.integer(ntree)
  nodesize <- as.integer(nodesize)
  
  # convert numeric index of classwt_label back to character label for scoring function
  classwt_label <- as.integer(classwt_label)
  label_str <- label_keys[classwt_label]
  classwt <- weight_grid[[label_str]]
  
  repeats <- 10
  repeat_auc <- numeric(repeats)
  
  # loop over each repeat
  for (r in 1:repeats) {
    
    #  create 5-folds for cross-validation (stratified on condition)
    folds <- caret::createFolds(metabo$condition, k = 5, list = TRUE, returnTrain = FALSE)
    fold_aucs <- numeric(length(folds))
    
    # loop through the folds
    for (f in seq_along(folds)) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
      test_data <- metabo[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               mtry = mtry,
                               ntree = ntree, 
                               nodesize = nodesize,
                               classwt = classwt,
                               importance = TRUE)
      
      # evaluate on test set
      probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # calculate AUC
      roc_obj <- tryCatch({
        roc(response = test_data$condition,
            predictor = probabilities[, "disease"],
            levels = c("healthy", "disease"),
            direction = "<")
      }, error = function(e) return(NULL))
      
      fold_aucs[f] <- if (!is.null(roc_obj)) auc(roc_obj) else NA
    }
    
    repeat_auc[r] <- mean(fold_aucs, na.rm = TRUE)
  }
  
  mean_auc <- mean(repeat_auc, na.rm = TRUE) # average across repeats
  return(list(Score = mean_auc)) # return mean auc values
}


# define parameter bounds
bounds <- list(mtry = c(1L, 30L),
               ntree = c(100L, 1000L),
               nodesize = c(1L, 6L),
               classwt_label = c(1L, length(label_keys))) # numeric range for classwt_label

# register back end
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 10,
                   iters.n = 10,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ()

# view results
optObj$stopStatus
score_summary = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
bestparams = getBestPars(optObj)


########################################################################################
########   RANDOM FOREST - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - NESTED CV  #######
########################################################################################

# data to be used in the model
str(metabo)
is.factor(metabo$condition) # condition needs to be a factor

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(metabo), "condition")


# create list of class weight settings
weight_grid <- list(healthy = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    disease = c(healthy = 1, disease = 2))

# define the set of categorical labels with numeric indices
label_keys <- c("healthy", "equal", "disease")


# set outer and inner loop parameters
outer_repeats <- 10
outer_folds <- 3
inner_repeats <- 10
inner_folds <- 3

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies
best_hyperparameters <-list() # list to store optimal hyperparameters

### outer loop
for (r in 1:outer_repeats) {
  
  # create outer_folds-folds for cross-validation (stratified on condition)
  set.seed(1234 + r*100)
  outer_folds_list <- createFolds(metabo$condition, k = outer_folds, list = TRUE)
  
  # loop through the folds
  for (f in 1:outer_folds) {
    cat("Fold ", f, " of Repeat ", r, "\n")
    
    # split dataset into training and testing sets for the current fold
    test_idx <- outer_folds_list[[f]] # test indices for the f-th fold
    train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
    test_data <- metabo[test_idx, ] # testing data (fold f)
    
    
    ### inner loop: Bayesian hyperparameter optimization
    
    # scoring function
    scoring_function <- function(mtry, ntree, nodesize, classwt_label) {
      
      # parameters
      mtry <- as.integer(mtry)
      ntree <- as.integer(ntree)
      nodesize <- as.integer(nodesize)
      
      # convert numeric index of classwt_label back to character label for scoring function
      classwt_label <- as.integer(classwt_label)
      label_str <- label_keys[classwt_label]
      classwt <- weight_grid[[label_str]]
      
      repeat_auc <- numeric(inner_repeats)
      
      # loop over each repeat
      for (inner_r in 1:inner_repeats) {
        
        # create folds for cross-validation (stratified on condition)
        inner_folds_list <- caret::createFolds(train_data$condition, k = inner_folds, list = TRUE, returnTrain = FALSE)
        fold_aucs <- numeric(length(inner_folds_list))
        
        # loop through the folds
        for (inner_f in seq_along(inner_folds_list)) {
          
          # splits the dataset into training and testing sets for the current fold
          inner_test_idx <- inner_folds_list[[inner_f]] # test indices for the f-th fold
          inner_train_data <- train_data[-inner_test_idx, ] # training data (all rows not in fold f)
          inner_test_data <- train_data[inner_test_idx, ] # testing data (fold f)
          
          # train random forest model using full features to rank features
          rf_model <- randomForest(x = inner_train_data[, subset_feat_cols], 
                                   y = as.factor(inner_train_data$condition),
                                   mtry = mtry,
                                   ntree = ntree, 
                                   nodesize = nodesize,
                                   classwt = classwt,
                                   importance = TRUE)
          
          # evaluate on test set
          probabilities <- predict(rf_model, newdata = inner_test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
          
          # calculate AUC
          roc_obj <- tryCatch({
            roc(response = inner_test_data$condition,
                predictor = probabilities[, "disease"],
                levels = c("healthy", "disease"),
                direction = "<")
          }, error = function(e) return(NULL))
          
          fold_aucs[inner_f] <- if (!is.null(roc_obj)) auc(roc_obj) else NA
        }
        
        repeat_auc[inner_r] <- mean(fold_aucs, na.rm = TRUE)
      }
      
      return(list(Score = mean(repeat_auc, na.rm = TRUE)))
    }
    
    # define parameter bounds
    bounds <- list(mtry = c(1L, 30L),
                   ntree = c(100L, 1000L),
                   nodesize = c(1L, 6L),
                   classwt_label = c(1L, length(label_keys))) # numeric range for classwt_label
    
    # register back end
    doParallel::registerDoParallel(parallel::detectCores() - 1)
    
    set.seed(1234 + r*1000 + f*10)
    optObj <- bayesOpt(FUN = scoring_function,
                       bounds = bounds,
                       initPoints = 10,
                       iters.n = 10,
                       acq = "ei",
                       parallel = TRUE,
                       verbose = 1)
    
    # unregister the backend
    registerDoSEQ()
    
    # get best hyperparameter values
    bestparams = getBestPars(optObj)
    best_hyperparameters[[paste0("Repeat_", r, "_Fold_", f)]] <- bestparams
    classwt_best <- weight_grid[[label_keys[bestparams$classwt_label]]]
    
    
    # train random forest model on the outer training data using best hyperparameters values
    rf_model <- randomForest(x = train_data[, subset_feat_cols],
                             y = as.factor(train_data$condition),
                             mtry = as.integer(bestparams$mtry),
                             ntree = as.integer(bestparams$ntree),
                             nodesize = as.integer(bestparams$nodesize),
                             classwt = classwt_best,
                             importance = TRUE) 
    
    # evaluate on outer test set
    test_predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
    test_probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
    # evaluate model on training set
    train_predictions <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "response")
    train_probabilities <- predict(rf_model, newdata = train_data[, subset_feat_cols], type = "prob")
    
    # calculate AUC on test set
    test_roc_obj <- roc(response = test_data$condition,
                        predictor = test_probabilities[, "disease"],
                        levels = c("healthy", "disease"),
                        direction = "<")
    test_auc <- auc(test_roc_obj)
    
    # store test ROC coordinates
    test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                              sensitivity = test_roc_obj$sensitivities,
                              Repeat = r, Fold = f, Set = "Test")
    
    # calculate AUC on train set
    train_roc_obj <- roc(response = train_data$condition,
                         predictor = train_probabilities[, "disease"],
                         levels = c("healthy", "disease"),
                         direction = "<")
    train_auc <- auc(train_roc_obj)
    
    # store train ROC coordinates
    train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                               sensitivity = train_roc_obj$sensitivities,
                               Repeat = r, Fold = f, Set = "Train")
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrices
    test_cm <- confusionMatrix(test_predictions, as.factor(test_data$condition), positive = "disease")
    train_cm <- confusionMatrix(train_predictions, as.factor(train_data$condition), positive = "disease")
    
    
    ### store with repeat (r) and fold (f) index
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                       train_cm = train_cm, train_auc = train_auc,
                                       test_roc_df = test_roc_df, train_roc_df = train_roc_df) # store performance metrics (test and train)
    feature_importances[[key]] <- importance(rf_model)  # store feature importances
  }
}


### calculate feature frequencies and relative frequency of feature selection
feature_split_summary <- bind_rows(feature_frequencies, .id = "Repeat_Fold") %>%
  rename(Feature = tree_split_vars) %>%
  group_by(Feature) %>% # group stability metrics by feature
  summarise(total_count = sum(Freq, na.rm = TRUE),
            mean_count = mean(Freq, na.rm = TRUE),
            n_models = n()) %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree)) %>%
  arrange(desc(total_count))
head(as.data.frame(feature_split_summary), 20)

# total number of models where feature was used at least once
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - split",
       x = "Feature", y = "Average splits per tree")


### calculate feature importances
# combine all feature_importances data.frames into one data.frame
all_features_importances <- do.call(rbind, lapply(names(feature_importances), function(name) {
  df <- as.data.frame(feature_importances[[name]])
  df$Feature <- rownames(df)
  df$Repeat_Fold <- name
  return(df)
}))

# group importance metrics by feature
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(as.data.frame(mean_importance), 20)

# plot features with highest MeanDecreaseAccuracy
ggplot(metabo, aes(x = beta.alanine)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Beta alanine",
       x = "Abundance", y = "Density of samples", fill = "Condition") +
  theme_minimal()


### combine feature stabilities and importances
feature_split_importance <- left_join(feature_split_summary, mean_importance, by = "Feature") %>%
  arrange(desc(avg_per_tree))
head(as.data.frame(feature_split_importance), 20)

# plot split frequency versus feature importance
# top 20 features by split frequency
top_features <- feature_split_importance %>% 
  slice(1:20) # top 20 features by split frequency

# plot meanMDA versus split frequency
ggplot(top_features, aes(x = avg_per_tree, y = mean_MeanDecreaseAccuracy, label = Feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA",
       title = "Split frequency versus permutation importance")

# long data.frame for plotting importance by condition
top_features_long <- top_features %>%
  select(Feature, avg_per_tree, mean_healthy, mean_disease) %>%
  pivot_longer(cols = c(mean_healthy, mean_disease),
               names_to = "Condition",
               values_to = "Importance") %>%
  mutate(Condition = recode(Condition,
                            mean_healthy = "Healthy",
                            mean_disease = "Disease"))

# plot meanMDA by condition versus split frequency
ggplot(top_features_long, aes(x = avg_per_tree, y = Importance, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA by condition",
       title = "Split frequency versus permutation importance by condition", color = "Condition") +
  scale_color_manual(values = c("Healthy" = "steelblue", "Disease" = "indianred3"))


### calculate overall performance statistics
perf_stats(performance_metrics, type = "test")
perf_stats(performance_metrics, type = "train")
perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
plot_roc(performance_metrics)


### best hyperparameters
hyperparameters <- do.call(rbind, lapply(names(best_hyperparameters), function(name) {
  df <- as.data.frame(t(best_hyperparameters[[name]]))
  df$Repeat_Fold <- name
  return(df)
}))
hyperparameters <- hyperparameters %>% select(Repeat_Fold, everything())
head(hyperparameters, 20)

# data frame of performance and hyperparameter values
hyper_perf <- do.call(rbind, lapply(names(performance_metrics), function(key) {
  perf <- performance_metrics[[key]]
  auc_test <- perf$test_auc
  auc_train <- perf$train_auc
  
  # get hyperparameters
  hyper <- best_hyperparameters[[key]]
  
  df <- data.frame(Repeat_Fold = key,
                   auc_test = auc_test,
                   auc_train = auc_train,
                   mtry = as.integer(hyper["mtry"]),
                   ntree = as.integer(hyper["ntree"]),
                   nodesize = as.integer(hyper["nodesize"]),
                   classwt_label = as.integer(hyper["classwt_label"]))
  return(df)
}))

hyper_perf

# long format for faceting
hyper_perf_long <- hyper_perf %>%
  pivot_longer(cols = c(mtry, ntree, nodesize, classwt_label),
               names_to = "hyperparameter",
               values_to = "param_value")

# facetted plot of hyperparameter versus test auc
ggplot(hyper_perf_long, aes(x = param_value, y = auc_test)) +
  geom_point(alpha = 0.7) + geom_smooth(method = "loess") +
  facet_wrap(~ hyperparameter, scales = "free_x") +  theme_minimal() +
  labs(x = "Hyperparameter value", y = "Test AUC",
       title = "Test AUC vs Hyperparameters")

hyper_cor <- sapply(c("mtry", "ntree", "nodesize", "classwt_label"), 
                    function(col) cor(hyper_perf$auc_test, hyper_perf[[col]], use = "complete.obs"))
hyper_cor


############################################################
########   RANDOM FOREST - BORUTA FEATURE SELECTION  #######
############################################################

# data to be used in the model
str(metabo)
is.factor(metabo$condition) # condition needs to be a factor

# set outer and inner loop parameters
outer_repeats <- 10
outer_folds <- 5
boruta_repeats <- 10

# create lists to store metrics
stable_feature_frequencies <- list() # list to store stable feature selection frequencies from Boruta
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# outer loop
for (r in 1:outer_repeats) {
  
  # create outer_folds-folds for cross-validation (stratified on condition)
  set.seed(1234 + r*100)
  outer_folds_list <- createFolds(metabo$condition, k = outer_folds, list = TRUE)
  
  # loop through the folds
  for (f in 1:outer_folds) {
    cat("Fold ", f, " of Repeat ", r, "\n")
    
    # split dataset into training and testing sets for the current fold
    test_idx <- outer_folds_list[[f]] # test indices for the f-th fold
    train_data <- metabo[-test_idx, ] # training data (all rows not in fold f)
    test_data <- metabo[test_idx, ] # testing data (fold f)
    
    ### Boruta feature selection
    boruta_metrics_list <- list() # list to store Boruta metrics
    
    for (b in 1:boruta_repeats) {
      cat("Running Boruta:", b, "\n")
      
      set.seed(1234 + r*1000 + f*100 + b)
      bor <- Boruta(condition ~ ., data = train_data, maxRuns = 250, ntree = 500, pValue = 0.01, doTrace = 0)
      bor <- TentativeRoughFix(bor)
      
      # boruta metrics
      imp_stats <- attStats(bor)
      imp_stats$Feature <- rownames(imp_stats)
      boruta_metrics_list[[b]] <- imp_stats
    }
    
    # aggregate median importance
    boruta_df <- do.call(rbind, boruta_metrics_list)
    
    # keep only confirmed features
    confirmed_boruta_df <- boruta_df %>%
      filter(decision == "Confirmed") %>%
      group_by(Feature) %>%
      summarise(mean_medianImp = mean(medianImp, na.rm = TRUE),
                sd_medianImp = sd(medianImp, na.rm = TRUE),
                count = n())
    
    # add proportion
    confirmed_boruta_df$proportion <- confirmed_boruta_df$count/boruta_repeats 
    
    # keep only stable features 
    stable_feats <- confirmed_boruta_df %>% 
      filter(proportion >= 0.4) 
    selected_features <- stable_feats$Feature
    length_feat <- length(selected_features)
    
    if (length(selected_features) < 2) { 
      cat(sprintf("Too few stable features: skipping Repeat %d Fold %d.\n", r, f)) 
      next 
    }
    
    # save stability info for features
    stable_feature_frequencies[[paste0("Repeat_", r, "_Fold_", f)]] <- stable_feats
    
    # subset train_data and test_data to just selected_features
    train_subset <- train_data[, c("condition", selected_features)]
    test_subset <- test_data[, c("condition", selected_features)]
    
    # column names for features to be included in model
    subset_feat_cols <- setdiff(colnames(train_subset), "condition")
    
    # train random forest model on the train_subset data
    set.seed(1234 + r*100000 + f*1000)
    rf_model <- randomForest(x = train_subset[, subset_feat_cols], 
                             y = as.factor(train_subset$condition), 
                             ntree = 500,
                             importance = TRUE) 
    
    # evaluate on test set
    test_predictions <- predict(rf_model, newdata = test_subset[, subset_feat_cols], type = "response") # predicted class labels for cm
    test_probabilities <- predict(rf_model, newdata = test_subset[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
    # evaluate model on training set
    train_predictions <- predict(rf_model, newdata = train_subset[, subset_feat_cols], type = "response")
    train_probabilities <- predict(rf_model, newdata = train_subset[, subset_feat_cols], type = "prob")
    
    # calculate AUC on test set
    test_roc_obj <- roc(response = test_subset$condition,
                        predictor = test_probabilities[, "disease"],
                        levels = c("healthy", "disease"),
                        direction = "<")
    test_auc <- auc(test_roc_obj)
    
    # store test ROC coordinates
    test_roc_df <- data.frame(specificity = test_roc_obj$specificities,
                              sensitivity = test_roc_obj$sensitivities,
                              Repeat = r, Fold = f, Set = "Test")
    
    # calculate AUC on train set
    train_roc_obj <- roc(response = train_subset$condition,
                         predictor = train_probabilities[, "disease"],
                         levels = c("healthy", "disease"),
                         direction = "<")
    train_auc <- auc(train_roc_obj)
    
    # store train ROC coordinates
    train_roc_df <- data.frame(specificity = train_roc_obj$specificities,
                               sensitivity = train_roc_obj$sensitivities,
                               Repeat = r, Fold = f, Set = "Train")
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrices
    test_cm <- confusionMatrix(test_predictions, as.factor(test_subset$condition), positive = "disease")
    train_cm <- confusionMatrix(train_predictions, as.factor(train_subset$condition), positive = "disease")
    
    
    ### store with repeat (r) and fold (f) index
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(test_cm = test_cm, test_auc = test_auc, 
                                       train_cm = train_cm, train_auc = train_auc,
                                       test_roc_df = test_roc_df, train_roc_df = train_roc_df) # store performance metrics (test and train)
    feature_importances[[key]] <- importance(rf_model)  # store feature importances
  }
}

# calculate frequency of Boruta feature selection
boruta_selected_summary <- bind_rows(stable_feature_frequencies, .id = "Repeat_Fold") %>%
  group_by(Feature) %>%
  summarise(n_selected_models = n(), # how many folds the feature appeared in
            overall_selection_rate = n_selected_models/length(stable_feature_frequencies), # stability of the feature
            mean_medianImp = mean(mean_medianImp)) %>% # average importance of feature
  arrange(desc(overall_selection_rate))
head(as.data.frame(boruta_selected_summary), 20)

# plot stable features ranked by average median importance
ggplot(boruta_selected_summary[1:30, ], aes(x = reorder(Feature, mean_medianImp), y = mean_medianImp)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Average median importance of confirmed features",
       x = "Feature", y = "Mean median importance")

# plot features selection frequency
ggplot(boruta_selected_summary[1:30, ], aes(x = reorder(Feature, overall_selection_rate), y = overall_selection_rate)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Proportion of Boruta repeats selecting each feature",
       x = "Feature", y = "Proportion of Runs")  


### calculate feature frequencies and relative frequency of feature selection
feature_split_summary <- bind_rows(feature_frequencies, .id = "Repeat_Fold") %>%
  rename(Feature = tree_split_vars) %>%
  group_by(Feature) %>% # group stability metrics by feature
  summarise(total_count = sum(Freq, na.rm = TRUE),
            mean_count = mean(Freq, na.rm = TRUE),
            n_models = n()) %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree)) %>%
  arrange(desc(total_count))
head(as.data.frame(feature_split_summary), 20)

# total number of models where feature was used at least once
ggplot(feature_split_summary[1:20, ], aes(x = reorder(Feature, n_models), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
ggplot(feature_split_summary[1:20, ], aes(x = reorder(Feature, avg_per_tree), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features - split",
       x = "Feature", y = "Average splits per tree")


### calculate feature importances
# combine all feature_importances data.frames into one data.frame
all_features_importances <- do.call(rbind, lapply(names(feature_importances), function(name) {
  df <- as.data.frame(feature_importances[[name]])
  df$Feature <- rownames(df)
  df$Repeat_Fold <- name
  return(df)
}))

# group importance metrics by feature
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(as.data.frame(mean_importance), 20)

# plot features with highest MeanDecreaseAccuracy
ggplot(metabo, aes(x = beta.alanine)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Beta alanine",
       x = "Abundance", y = "Density of samples", fill = "Condition") +
  theme_minimal()


### combine feature stabilities and importances
feature_split_importance <- left_join(feature_split_summary, mean_importance, by = "Feature") %>%
  arrange(desc(avg_per_tree))
head(as.data.frame(feature_split_importance), 20)

# plot split frequency versus feature importance
# top 20 features by split frequency
top_features <- feature_split_importance %>% 
  slice(1:20) # top 20 features by split frequency

# plot meanMDA versus split frequency
ggplot(top_features, aes(x = avg_per_tree, y = mean_MeanDecreaseAccuracy, label = Feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA",
       title = "Split frequency versus permutation importance")

# long data.frame for plotting importance by condition
top_features_long <- top_features %>%
  select(Feature, avg_per_tree, mean_healthy, mean_disease) %>%
  pivot_longer(cols = c(mean_healthy, mean_disease),
               names_to = "Condition",
               values_to = "Importance") %>%
  mutate(Condition = recode(Condition,
                            mean_healthy = "Healthy",
                            mean_disease = "Disease"))

# plot meanMDA by condition versus split frequency
ggplot(top_features_long, aes(x = avg_per_tree, y = Importance, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) + theme_minimal() +
  labs(x = "Split frequency", y = "Mean MDA by condition",
       title = "Split frequency versus permutation importance by condition", color = "Condition") +
  scale_color_manual(values = c("Healthy" = "steelblue", "Disease" = "indianred3"))


### calculate performance statistics
perf_stats(performance_metrics, type = "test")
perf_stats(performance_metrics, type = "train")
perf_stats(performance_metrics, type = "gap")


### plot average ROC curve across folds
plot_roc(performance_metrics)


################################################################################
########   RANDOM FOREST MODEL - TRAIN FINAL MODEL WITH BORUTA FEATURES  #######
################################################################################

# features to be used in the model
boruta_features_df <- boruta_selected_summary %>%
  arrange(desc(overall_selection_rate)) %>%
  slice_head(n = 20)
boruta_features <- boruta_features_df$Feature


# train the model on the full dataset
set.seed(1234)
final_model <- randomForest(x = metabo[, boruta_features], 
                            y = as.factor(metabo$condition), 
                            ntree = 500,
                            importance = TRUE)

# # save final model
# saveRDS(final_model, file = "final_rf_model_metabo.rds")


#########################################################################
###   RANDOM FOREST - SHAP VALUES - DEPENDENCE AND INTERACTION PLOTS  ###
#########################################################################

# function to get predicted disease probability
pred_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "disease"]
}

# compute fast SHAP values
set.seed(1234) # set seed

shap_values <- fastshap::explain(object = final_model,
                                 X = metabo[, boruta_features],
                                 pred_wrapper = pred_fun,
                                 nsim = 100, # number of Monte Carlo permutations
                                 adjust = TRUE) # centered SHAP values (baseline + sum of SHAPs = prediction)
shap_values <- as.data.frame(shap_values)     


### SHAP dependence plots (how SHAP values for a given feature vary as the input values for the feature vary)
dependence_plot <- function(shap_values, feature_matrix, feature_name, interaction_feature) {
  
  # rfvalues (abundance) wide
  rfvalue_wide <- feature_matrix %>%
    mutate(ID = row_number())
  
  # rfvalues (abundance) long
  rfvalue_long <- rfvalue_wide %>%
    pivot_longer(cols = -ID, names_to = "feature", values_to = "rfvalue")
  
  # summarize SHAP values
  shap_long <- shap_values %>%
    mutate(ID = row_number()) %>%
    pivot_longer(-ID, names_to = "feature", values_to = "value") %>%
    left_join(rfvalue_long, by = c("ID", "feature")) # add rfvalues (abundance) to shap_long
  
  # filter to target (feature_name)
  shap_dep <- shap_long %>% filter(feature == feature_name)
  
  # add interaction feature values
  plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
  
  
  p <- ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + 
    theme_minimal() + geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
    scale_color_distiller(palette = "RdBu", direction = 1) +
    labs(title = paste("SHAP dependence plot for", feature_name),
         x = paste(feature_name, "- log abun"), 
         y = paste(feature_name, "- SHAP value"),
         color = paste(interaction_feature, "- log abun"))
  
  return(p)
}

dependence_plot(shap_values = shap_values, feature_matrix = metabo[, boruta_features], 
                feature_name = "beta.alanine", 
                interaction_feature = "phosphate")

dependence_plot(shap_values = shap_values, feature_matrix = metabo[, boruta_features], 
                feature_name = "myristate..14.0.", 
                interaction_feature = "X.24683")


### bee swarm plot
beeswarm_plot <- function(shap_values, feature_matrix) {
  
  # rfvalues (abundance) wide
  rfvalue_wide <- feature_matrix %>%
    mutate(ID = row_number())
  
  # rfvalues (abundance) long
  rfvalue_long <- rfvalue_wide %>%
    pivot_longer(cols = -ID, names_to = "feature", values_to = "rfvalue")
  
  # summarize SHAP values
  shap_long <- shap_values %>%
    mutate(ID = row_number()) %>%
    pivot_longer(-ID, names_to = "feature", values_to = "value") %>%
    left_join(rfvalue_long, by = c("ID", "feature")) # add rfvalues (abundance) to shap_long
  
  # compute mean abs shap for ordering
  feature_order <- shap_long %>%
    group_by(feature) %>%
    summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_shap)) %>%
    pull(feature)
  
  # data for plotting
  shap_plot <- shap_long %>%
    group_by(feature) %>%
    mutate(scaled_rfvalue = scale(rfvalue)[,1]) %>%
    ungroup() %>%
    mutate(feature = factor(feature, levels = rev(feature_order))) 
  
  p <- ggplot(shap_plot, aes(x = value, y = feature, color = scaled_rfvalue)) +
    geom_jitter(height = 0.2, alpha = 0.7, size = 1.2) +
    scale_color_viridis_c(option = "plasma", direction = -1) + theme_minimal() +
    labs(title = "SHAP summary plot", x = "SHAP value (impact on model output)", 
         color = "Feature value")
  
  return(p)
}

beeswarm_plot(shap_values = shap_values, 
              feature_matrix = metabo[, boruta_features])


### mean and absolute mean SHAP values
shap_values_plot <- function(shap_values, type = "mean_abs_shap") {
  
  # summarize SHAP values
  shap_long <- shap_values %>%
    mutate(ID = row_number()) %>%
    pivot_longer(-ID, names_to = "feature", values_to = "value") 
  
  # mean and mean abs shap values
  shap_summary <- shap_long %>%
    group_by(feature) %>%
    summarise(mean_shap = mean(value, na.rm = TRUE),
              mean_abs_shap = mean(abs(value), na.rm = TRUE),
              .groups = "drop")
  
  # plot depending on type
  if (type == "mean_abs_shap") {
    p <- ggplot(shap_summary, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
      geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
      labs(title = "Mean absolute SHAP value", 
           x = "Feature", y = "Mean absolute SHAP value")
    
  } else if (type == "mean_shap") {
    p <- ggplot(shap_summary, aes(x = reorder(feature, mean_shap), y = mean_shap)) +
      geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
      labs(title = "Mean SHAP value", 
           x = "Feature", y = "Mean SHAP value")
    
  } else {
    stop("type must be either mean_shap or mean_abs_shap")
  }
  
  return(p)
}

shap_values_plot(shap_values = shap_values, type = "mean_shap")
shap_values_plot(shap_values = shap_values, type = "mean_abs_shap")


### mean or mean abs SHAP values versus MeanDecreaseAccuracy or MeanDecreaseGini
shap_importance <- function(final_model, shap_values, 
                            measure = c("MeanDecreaseAccuracy", "MeanDecreaseGini"),
                            type = c("mean_abs_shap", "mean_shap")) {
  
  # validate arguments
  measure <- match.arg(measure)
  type <- match.arg(type)
  
  # summarize SHAP values
  shap_long <- shap_values %>%
    mutate(ID = row_number()) %>%
    pivot_longer(-ID, names_to = "feature", values_to = "value") 
  
  # mean shap values
  shap_summary <- shap_long %>%
    group_by(feature) %>%
    summarise(mean_shap = mean(value, na.rm = TRUE),
              mean_abs_shap = mean(abs(value), na.rm = TRUE),
              .groups = "drop")
  
  # get importance and join to shap values
  importance <- as.data.frame(importance(final_model))
  importance$feature <- rownames(importance)
  shap_import_df <- shap_summary %>%
    left_join(importance, by = "feature")
  
  # dynamic axis labels 
  x_label <- if (measure == "MeanDecreaseAccuracy") "Mean decrease accuracy" else "Mean decrease in Gini"
  y_label <- if (type == "mean_abs_shap") "Mean absolute SHAP value" else "Mean SHAP value"
  title_label <- paste(y_label, "versus", tolower(x_label))
  
  # plot dynamically
  p <- ggplot(shap_import_df, aes(x = .data[[measure]], y = .data[[type]], label = feature)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.8) +
    geom_text_repel(size = 3) + theme_minimal() +
    labs(title = title_label, x = x_label, y = y_label)
  
  return(p)
}

shap_importance(final_model = final_model, shap_values = shap_values, 
                measure = "MeanDecreaseAccuracy", type = "mean_abs_shap")
shap_importance(final_model = final_model, shap_values = shap_values, 
                measure = "MeanDecreaseGini", type = "mean_abs_shap")


### SHAP values by label/class
shap_by_class <- function(shap_values, dataset, label, type = "mean_abs_shap") {
  
  # mean and mean abs shap value per class 
  shap_values[[label]] <- dataset[[label]]
  shap_by_class_df <- shap_values %>%
    pivot_longer(cols = -all_of(label)) %>%
    group_by(.data[[label]], name) %>%
    summarise(mean_shap = mean(value),
              mean_abs_shap = mean(abs(value)),
              .groups = "drop")
  
  if (type == "mean_abs_shap") {
    p <- ggplot(shap_by_class_df, aes(x = reorder(name, mean_abs_shap), y = mean_abs_shap, fill = .data[[label]])) +
      geom_col(position = "dodge") + coord_flip() + theme_minimal() +
      labs(title = "Class-specific mean absolute SHAP values",
           x = "Feature", y = "Mean abs SHAP", fill = label)
    
  } else if (type == "mean_shap") {
    p <- ggplot(shap_by_class_df, aes(x = reorder(name, mean_shap), y = mean_shap, fill = .data[[label]])) +
      geom_col(position = "dodge") + coord_flip() + theme_minimal() +
      labs(title = "Class-specific mean SHAP values",
           x = "Feature", y = "Mean SHAP", fill = label)
    
  } else {
    stop("type must be either mean_abs_shap or mean_shap")
  }
  
  return(p)
}

shap_by_class(shap_values = shap_values, dataset = metabo, 
              label = "condition", type = "mean_abs_shap")
shap_by_class(shap_values = shap_values, dataset = metabo, 
              label = "condition", type = "mean_shap")


### feature metrics
feature_metrics <- function(shap_values, dataset, label, feature_name) {
  
  # add label to shap_values
  shap_values[[label]] <- dataset[[label]]
  
  # ensure feature exists in both data frames
  if (!feature_name %in% colnames(shap_values)) {
    stop(paste("Feature", feature_name, "not found in shap_values"))
  }
  if (!feature_name %in% colnames(dataset)) {
    stop(paste("Feature", feature_name, "not found in dataset"))
  }
  
  # density plot of SHAP values
  p1 <- ggplot(shap_values, aes(x = .data[[feature_name]], fill = .data[[label]])) +
    geom_density(alpha = 0.6) + theme_minimal() +
    labs(title = "SHAP value distribution",
         x = "SHAP value", y = "Density", fill = label)
  
  # boxplot of SHAP values by label
  p2 <- ggplot(shap_values, aes(x = .data[[label]], y = .data[[feature_name]], fill = .data[[label]])) +
    geom_boxplot() + theme_minimal() +
    labs(title = paste("SHAP values by", label),
         y = "SHAP value", x = label)
  
  # box plot of log-abundance by label
  p3 <- ggplot(dataset, aes(x = .data[[label]], y = .data[[feature_name]], fill = .data[[label]])) +
    geom_boxplot() + theme_minimal() +
    labs(title = paste("Log abundance by", label),
         y = "Log Abundance", x = label)
  
  # combine plots vertically
  combined_plot <- (p1 / p2 / p3) +
    plot_annotation(title = paste("Feature metrics for", feature_name))
  
  return(combined_plot)
}

feature_metrics(shap_values = shap_values, dataset = metabo, label = "condition",
                feature_name = "beta.alanine")
feature_metrics(shap_values = shap_values, dataset = metabo, label = "condition",
                feature_name = "myristate..14.0.")


### SHAP values versus predicted probabilities
shap_pred_probs <- function(final_model, shap_values, dataset, feature_vector, label, positive_class, feature_name) {
  
  # compute predicted probabilities from the model
  prob <- as.data.frame(predict(final_model, newdata = dataset[, feature_vector], type = "prob"))
  
  # add predicted probability for positive class
  shap_values$pred_prob <- prob[[positive_class]]
  
  # add label to shap_values
  shap_values[[label]] <- dataset[[label]]
  
  p <- ggplot(shap_values, aes(x = .data[[feature_name]], y = pred_prob, color = .data[[label]])) +
    geom_point(alpha = 0.6) + theme_minimal() +
    labs(title = paste(feature_name,"SHAP value versus prediction probability"),
         x = paste("SHAP value for", feature_name), y = paste("Predicted probability (", positive_class,")"), 
         color = label)
  
  return(p)
}

shap_pred_probs(final_model = final_model, shap_values = shap_values, 
                dataset = metabo, feature_vector = boruta_features, 
                label = "condition", positive_class = "disease", 
                feature_name = "beta.alanine")

shap_pred_probs(final_model = final_model, shap_values = shap_values, 
                dataset = metabo, feature_vector = boruta_features, 
                label = "condition", positive_class = "disease", 
                feature_name = "myristate..14.0.")


### fastshap does not calculate second-order SHAP effects
# partial dependence (PDP) interaction plots (how much each feature interacts with others)
library(iml)
predictor <- Predictor$new(final_model, data = metabo[, boruta_features], y = metabo$condition, type = "prob")
interaction_strength <- Interaction$new(predictor)
plot(interaction_strength)


sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Tahoe 26.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Edmonton
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      tools     stats4    splines   compiler  parallel  stats     graphics  grDevices utils     datasets 
# [12] methods   base     
# 
# other attached packages:
#   [1] MASS_7.3-65                   hardhat_1.4.2                 lifecycle_1.0.4              
# [4] farver_2.1.2                  digest_0.6.39                 rstatix_0.7.3                
# [7] gtable_0.3.6                  lava_1.8.2                    cli_3.6.5                    
# [10] Formula_1.2-5                 ipred_0.9-15                  ggsignif_0.6.4               
# [13] gower_1.0.2                   ModelMetrics_1.2.2.2          data.table_1.17.8            
# [16] glue_1.8.0                    class_7.3-23                  globals_0.18.0               
# [19] scales_1.4.0                  hms_1.1.4                     dbscan_1.2.4                 
# [22] generics_0.1.4                checkmate_2.3.3               ggpubr_0.6.2                 
# [25] pillar_1.11.1                 proxy_0.4-28                  survival_3.8-3               
# [28] S7_0.2.1                      withr_3.0.2                   plyr_1.8.9                   
# [31] listenv_0.10.0                codetools_0.2-20              timeDate_4051.111            
# [34] abind_1.4-8                   dichromat_2.0-0.1             rstudioapi_0.17.1            
# [37] tidyselect_1.2.1              timechange_0.3.0              nnet_7.3-20                  
# [40] Matrix_1.7-4                  Metrics_0.1.4                 future.apply_1.20.1          
# [43] Rcpp_1.1.0                    cellranger_1.1.0              rpart_4.1.24                 
# [46] car_3.1-3                     carData_3.0-5                 parallelly_1.46.0            
# [49] ranger_0.17.0                 RColorBrewer_1.1-3            stringi_1.8.7                
# [52] R6_2.6.1                      broom_1.0.11                  recipes_1.3.1                
# [55] tzdb_0.5.0                    prodlim_2025.04.28            utf8_1.2.6                   
# [58] labeling_0.4.3                backports_1.5.0               crayon_1.5.3                 
# [61] pkgconfig_2.0.3               lhs_1.2.0                     reshape2_1.4.5               
# [64] vctrs_0.6.5                   mgcv_1.9-4                    nlme_3.1-168                 
# [67] e1071_1.7-17                  magrittr_2.0.4                rlang_1.1.6                  
# [70] gridExtra_2.3                 future_1.68.0                 iml_0.11.4                   
# [73] DiceKriging_1.6.1             ggrepel_0.9.6                 viridis_0.6.5                
# [76] viridisLite_0.4.2             fastshap_0.1.1                doParallel_1.0.17            
# [79] iterators_1.0.14              foreach_1.5.2                 ParBayesianOptimization_1.2.6
# [82] Boruta_9.0.0                  pROC_1.19.0.1                 caret_7.0-1                  
# [85] lattice_0.22-7                randomForest_4.7-1.2          impute_1.84.0                
# [88] readxl_1.4.5                  lubridate_1.9.4               forcats_1.0.1                
# [91] stringr_1.6.0                 dplyr_1.1.4                   purrr_1.2.0                  
# [94] readr_2.1.6                   tidyr_1.3.1                   tibble_3.3.0                 
# [97] tidyverse_2.0.0               patchwork_1.3.2               cowplot_1.2.0                
# [100] ggplot2_4.0.1

