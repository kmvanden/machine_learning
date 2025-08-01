# Random Forest Workflow - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(cowplot)
library(tidyverse)
library(compositions)
library(randomForest)
library(caret)
library(pROC)
library(e1071)
library(stringr)
library(Boruta)
library(ggpubr)
library(doParallel)
library(foreach)
library(ParBayesianOptimization)

# set.seed
set.seed(1234)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
str(meta)

# convert condition column into a factor
meta$condition <- as.factor(meta$condition)
table(meta$condition)


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

### filter species present in less than 10% of samples
feat_t <- as.data.frame(t(feat)) # transpose feature table
min_prevalence <- 0.10
feat_filtered <- feat_t[, colMeans(feat_t > 0) >= min_prevalence] # subset the feature table to only include features present in at least 10% of samples
dim(feat_filtered) # 70 935

# convert feature table to relative abundances
feat_rel_abund <- feat_filtered/rowSums(feat_filtered)

# add pseudocount and perform CLR transformation
feat_rel_abund[feat_rel_abund == 0] <- 1e-6
feat_clr <- clr(feat_rel_abund)
feat <- as.data.frame(feat_clr)


# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat))
feat$sample_name <- rownames(feat) # add column sample_name

### merge metadata and feature table
metagen <- merge(meta, feat, by = "sample_name", all.x = TRUE)
metagen <- metagen[,-1] # remove sample_name

# make sure names are syntactically valid 
colnames(metagen) <- make.names(colnames(metagen))


########################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS   ###
########################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(metagen$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- metagen[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by all_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, all_feat_cols], 
                             y = as.factor(train_data$condition), 
                             ntree = 500, importance = TRUE) 
    
    # evaluate on test set
    predictions <- predict(rf_model, newdata = test_data[, all_feat_cols])
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrix
    cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
    
    # store with repeat (r) and fold (f) index
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- cm # store performance metrics
    feature_importances[[key]] <- importance(rf_model)  # store feature importances
  }
}

### calculate feature frequencies
all_splits <- bind_rows(feature_frequencies, .id = "Repeat_Fold") # combine frequencies into a single data.frame
colnames(all_splits) <- c("Repeat_Fold", "Feature", "Count") # rename columns

# summarize total and average counts
feature_split_summary <- all_splits %>%
  group_by(Feature) %>%
  summarise(total_count = sum(Count, na.rm = TRUE),
            mean_count = mean(Count, na.rm = TRUE),
            n_models = n()) %>%
  arrange(desc(total_count))
head(feature_split_summary, 20)

# calculate relative frequency of feature selection
feature_split_summary <- feature_split_summary %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree))

# total number of models where feature was used at least once
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features by RF-based RFE",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features by RF-based RFE",
       x = "Feature", y = "Number of models")


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()

# extract metrics from the stored confusion matrices (50 repeats x 5 folds = 250 values)
for (cm in performance_metrics) {
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE))
metric_summary


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

### plot species with highest MeanDecreaseAccuracy
ggplot(metagen, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(metagen, aes(x = Anaerobutyricum_hallii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Anaerobutyricum hallii",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(metagen, aes(x = Clostridium_sp._M62.1)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Clostridium sp. M62.1",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()


#######################################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - RECURSIVE FEATURE ELIMINATION  ###
#######################################################################################################


# microbiome datasets have many weakly informative features --> RFE helps to find combinations that are informative together
# recursive feature selection using MeanDecreaseAccuracy
# full model trained on all features + extract top-ranked features based on MeanDecreaseAccuracy
# iterate over the defined feature sizes, training a model using only the top N features


# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# feature sizes to use in recursive feature selection
feature_sizes <- c(10, 25, 50, 75, 100, 200, 400, 600, 800, 935)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies
feature_importances <- list() # list to store feature importances

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(metagen$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- metagen[test_idx, ] # testing data (fold f)
    
    # train random forest model using full features to rank features
    full_rf <- randomForest(x = train_data[, all_feat_cols], 
                            y = as.factor(train_data$condition),
                            ntree = 500, importance = TRUE)
    
    # extract feature importances and arrange by MeanDecreaseAccuracy for feature selection
    full_importance <- importance(full_rf)
    importance_df <- as.data.frame(full_importance)
    importance_df$Feature <- rownames(importance_df)
    
    top_features_by_acc <- importance_df %>%
      arrange(desc(MeanDecreaseAccuracy)) %>%
      pull(Feature)
    
    # loop over feature size sets
    for (n_feat in feature_sizes) {
      
      selected_feats <- top_features_by_acc[1:n_feat]
      
      # train random forest model on selected features
      rf_model <- randomForest(x = train_data[, selected_feats], 
                               y = as.factor(train_data$condition), 
                               ntree = 500, importance = TRUE) 
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, selected_feats], type = "response") # predicted class labels for cm
      probabilities <- predict(rf_model, newdata = test_data[, selected_feats], type = "prob") # class probabilities (ROC/AUC)
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # count how often each feature is used in the trees
      tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
        tree <- getTree(rf_model, k = t, labelVar = TRUE)
        as.character(tree$`split var`[tree$`split var` != "<leaf>"])
      }))
      # count the occurrences of each feature 
      split_counts <- table(tree_split_vars)
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # store with repeat (r) and fold (f) index
      key <- paste0("Repeat_", r, "_Fold_", f, "_N_", n_feat)
      feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
      feature_importances[[key]] <- importance(rf_model) # store feature importances
    }
  }
}

### calculate feature frequencies
all_splits <- bind_rows(feature_frequencies, .id = "Repeat_Fold") # combine frequencies into a single data.frame
colnames(all_splits) <- c("Repeat_Fold", "Feature", "Count") # rename columns

# summarize total and average counts
feature_split_summary <- all_splits %>%
  group_by(Feature) %>%
  summarise(total_count = sum(Count, na.rm = TRUE),
            mean_count = mean(Count, na.rm = TRUE),
            n_models = n()) %>%
  arrange(desc(total_count))
head(feature_split_summary, 20)

# calculate relative frequency of feature selection
feature_split_summary <- feature_split_summary %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree))

# total number of models where feature was used at least once
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features by RF-based RFE",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features by RF-based RFE",
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

### plot species with highest MeanDecreaseAccuracy
ggplot(metagen, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(metagen, aes(x = Anaerobutyricum_hallii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Anaerobutyricum hallii",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(metagen, aes(x = Clostridium_sp._M62.1)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Clostridium sp. M62.1",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

### calculate performance statistics
# create vectors to store metrics
n_feat <- numeric()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract feature size from the key (e.g., "Repeat_1_Fold_3_N_100" → 100)
  n <- as.numeric(str_extract(key, "(?<=N_)\\d+"))
  
  n_feat <- c(n_feat, n)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc)
}

# combine metrics in a summary table
results_df <- data.frame(feature_count = n_feat,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         auc = auc_vals)

# group by feature_count and summarise
metric_summary <- results_df %>%
  group_by(feature_count) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            sd_bal_acc = sd(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            sd_f1 = sd(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            sd_sens = sd(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            sd_spec = sd(spec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(mean_bal_acc))
metric_summary

# balanced accuracy
ggplot(metric_summary, aes(x = feature_count, y = mean_bal_acc)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Performance versus number of selected features",
       x = "Number of features", y = "Balanced accuracy") +
  geom_errorbar(aes(ymin = mean_bal_acc - sd_bal_acc, ymax = mean_bal_acc + sd_bal_acc), width = 10, alpha = 0.4)

# f1 score
ggplot(metric_summary, aes(x = feature_count, y = mean_f1)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Performance versus number of selected features",
       x = "Number of features", y = "F1 score") +
  geom_errorbar(aes(ymin = mean_f1 - sd_f1, ymax = mean_f1 + sd_f1), width = 10, alpha = 0.4)

# sensitivity
ggplot(metric_summary, aes(x = feature_count, y = mean_sens)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Performance versus number of selected features",
       x = "Number of features", y = "Sensitivity") +
  geom_errorbar(aes(ymin = mean_sens - sd_sens, ymax = mean_sens + sd_sens), width = 10, alpha = 0.4)

# specificity
ggplot(metric_summary, aes(x = feature_count, y = mean_spec)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Performance versus number of selected features",
       x = "Number of features", y = "Specificity") +
  geom_errorbar(aes(ymin = mean_spec - sd_spec, ymax = mean_spec + sd_spec), width = 10, alpha = 0.4)

# auc
ggplot(metric_summary, aes(x = feature_count, y = mean_auc)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Performance versus number of selected features",
       x = "Number of features", y = "AUC") +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), width = 10, alpha = 0.4)


#######################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - CLASS WEIGHTS  ###
#######################################################################################


# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# create list of class weight settings
weight_grid <- list(high.h = c(healthy = 3, disease = 1),
                    med.h = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    med.d = c(healthy = 1, disease = 2))

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop for class weight
for (w in names(weight_grid)) {
  classwt <- weight_grid[[w]]
  cat("Class Weight:", w, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metagen$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- metagen[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, all_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, importance = TRUE, classwt = classwt)
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
      probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # store with repeat (r) and fold (f) index
      key <- paste0(w, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
    }
  }
}

### calculate performance statistics
# create vectors to store metrics
weight_setting <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract weight setting from key name
  cw <- str_extract(key, "^[^_]+_[^_]+_[^_]+") # extract everything before the first underscore
  
  weight_setting <- c(weight_setting, cw)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc)
}

# combine metrics in a summary table
results_df <- data.frame(weight_setting = weight_setting,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         auc = auc_vals)

# group by weight_setting and summarise
metric_summary <- results_df %>%
  group_by(weight_setting) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            sd_bal_acc = sd(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            sd_f1 = sd(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            sd_sens = sd(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            sd_spec = sd(spec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(weight_setting))
metric_summary

# add weight column for plotting
metric_summary$weight <- str_extract(metric_summary$weight_setting, "^[^_]+")
metric_summary$weight <- factor(metric_summary$weight, levels = c("high.h", "med.h", "equal", "med.d"))


# balanced accuracy
ggplot(metric_summary, aes(x = weight, y = mean_bal_acc, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("high.h" = "skyblue", "med.h" = "pink", "equal" = "plum", "med.d" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "Balanced accuracy", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# f1 score
ggplot(metric_summary, aes(x = weight, y = mean_f1, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("high.h" = "skyblue", "med.h" = "pink", "equal" = "plum", "med.d" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "F1 score", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# sensitivity
ggplot(metric_summary, aes(x = weight, y = mean_sens, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("high.h" = "skyblue", "med.h" = "pink", "equal" = "plum", "med.d" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "Sensitivity", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# specificity
ggplot(metric_summary, aes(x = weight, y = mean_spec, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("high.h" = "skyblue", "med.h" = "pink", "equal" = "plum", "med.d" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "Specificity", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# auc
ggplot(metric_summary, aes(x = weight, y = mean_auc, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("high.h" = "skyblue", "med.h" = "pink", "equal" = "plum", "med.d" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "AUC", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")


##########################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - NUMBER OF TREES   ###
##########################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# ntree values to test
ntree_values <- c(500, 1000, 2000)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop for ntree values
for (n in ntree_values) {
  cat("ntree =", n, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metagen$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- metagen[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, all_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = n, importance = TRUE)
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
      probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # store with repeat (r) and fold (f) index
      key <- paste0("ntree_", n, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
    }
  }
}

### calculate performance statistics
# create vectors to store metrics
ntree_number <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract ntree name
  nt <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+") 
  
  ntree_number <- c(ntree_number, nt)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc)
}

# combine metrics in a summary table
results_df <- data.frame(ntree_number = ntree_number,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         auc = auc_vals)

# summary of performance metrics
metric_summary <- results_df %>%
  group_by(ntree_number) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            sd_bal_acc = sd(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            sd_f1 = sd(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            sd_sens = sd(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            sd_spec = sd(spec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(ntree_number)

# add ntree column for plotting
metric_summary$ntree <- str_extract(metric_summary$ntree_number, "(?<=ntree_)\\d+")
metric_summary$ntree <- factor(metric_summary$ntree, levels = c("500", "1000", "2000"))

# balanced accuracy
ggplot(metric_summary, aes(x = ntree, y = mean_bal_acc, fill = ntree)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("500" = "skyblue", "1000" = "pink", "2000" = "plum")) +
  labs(title = "Performance versus number of trees", 
       y = "Balanced accuracy", x = "Number of trees") + 
  theme_minimal() + theme(legend.position = "none")

# f1 score
ggplot(metric_summary, aes(x = ntree, y = mean_f1, fill = ntree)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("500" = "skyblue", "1000" = "pink", "2000" = "plum")) +
  labs(title = "Performance versus number of trees", 
       y = "F1 score", x = "Number of trees") + 
  theme_minimal() + theme(legend.position = "none")

# sensitivity
ggplot(metric_summary, aes(x = ntree, y = mean_sens, fill = ntree)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("500" = "skyblue", "1000" = "pink", "2000" = "plum")) +
  labs(title = "Performance versus number of trees", 
       y = "Sensitivity", x = "Number of trees") + 
  theme_minimal() + theme(legend.position = "none")

# specificity
ggplot(metric_summary, aes(x = ntree, y = mean_spec, fill = ntree)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("500" = "skyblue", "1000" = "pink", "2000" = "plum")) +
  labs(title = "Performance versus number of trees", 
       y = "Specificity", x = "Number of trees") + 
  theme_minimal() + theme(legend.position = "none")

# auc
ggplot(metric_summary, aes(x = ntree, y = mean_auc, fill = ntree)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("500" = "skyblue", "1000" = "pink", "2000" = "plum")) +
  labs(title = "Performance versus number of trees", 
       y = "AUC", x = "Number of trees") + 
  theme_minimal() + theme(legend.position = "none")


######################################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - NUMBER OF FEATURES AT SPLIT   ###
######################################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# mtry values to test
mtry_values <- c(15, 30, 60, 120, 240)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop for mtry values
for (m in mtry_values) {
  cat("mtry =", m, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metagen$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- metagen[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, all_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, importance = TRUE, mtry = m)
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
      probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # store with repeat (r) and fold (f) index
      key <- paste0("mtry_", m, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
    }
  }
}

### calculate performance statistics
# create vectors to store metrics
mtry_number <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract mtry name
  mt <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+") 
  
  mtry_number <- c(mtry_number, mt)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc)
}

# combine metrics in a summary table
results_df <- data.frame(mtry_number = mtry_number,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         auc = auc_vals)

# summary of performance metrics
metric_summary <- results_df %>%
  group_by(mtry_number) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            sd_bal_acc = sd(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            sd_f1 = sd(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            sd_sens = sd(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            sd_spec = sd(spec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(mtry_number)

# add ntree column for plotting
metric_summary$mtry <- str_extract(metric_summary$mtry_number, "(?<=mtry_)\\d+")
metric_summary$mtry <- factor(metric_summary$mtry, levels = c("15", "30", "60", "120", "240"))

# balanced accuracy
ggplot(metric_summary, aes(x = mtry, y = mean_bal_acc, fill = mtry)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("15" = "skyblue", "30" = "deepskyblue3", "60" = "skyblue3", "120" = "deepskyblue4", "240" = "skyblue4")) +
  labs(title = "Performance versus number of features at splits", 
       y = "Balanced accuracy", x = "Number of features at splits") + 
  theme_minimal() + theme(legend.position = "none")


# f1 score
ggplot(metric_summary, aes(x = mtry, y = mean_f1, fill = mtry)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("15" = "skyblue", "30" = "deepskyblue3", "60" = "skyblue3", "120" = "deepskyblue4", "240" = "skyblue4")) +
  labs(title = "Performance versus number of features at splits", 
       y = "F1 score", x = "Number of features at splits") + 
  theme_minimal() + theme(legend.position = "none")

# sensitivity
ggplot(metric_summary, aes(x = mtry, y = mean_sens, fill = mtry)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("15" = "skyblue", "30" = "deepskyblue3", "60" = "skyblue3", "120" = "deepskyblue4", "240" = "skyblue4")) +
  labs(title = "Performance versus number of features at splits", 
       y = "Sensitivity", x = "Number of features at splits") + 
  theme_minimal() + theme(legend.position = "none")

# specificity
ggplot(metric_summary, aes(x = mtry, y = mean_spec, fill = mtry)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("15" = "skyblue", "30" = "deepskyblue3", "60" = "skyblue3", "120" = "deepskyblue4", "240" = "skyblue4")) +
  labs(title = "Performance versus number of features at splits", 
       y = "Specificity", x = "Number of features at splits") + 
  theme_minimal() + theme(legend.position = "none")

# auc
ggplot(metric_summary, aes(x = mtry, y = mean_auc, fill = mtry)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("15" = "skyblue", "30" = "deepskyblue3", "60" = "skyblue3", "120" = "deepskyblue4", "240" = "skyblue4")) +
  labs(title = "Performance versus number of features at splits", 
       y = "AUC", x = "Number of features at splits") + 
  theme_minimal() + theme(legend.position = "none")


####################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - NODE SIZE   ###
####################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# nodesize values to test
nodesize_values <- c(1, 5, 10, 20, 30)

# create list to store performance metrics
performance_metrics <- list() # list to store performance metrics

# loop for nodesize values
for (n in nodesize_values) {
  cat("nodesize =", n, "\n")
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    # create 5-folds for cross-validation (stratified on condition)
    folds <- createFolds(metagen$condition, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- metagen[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, all_feat_cols], 
                               y = as.factor(train_data$condition),
                               ntree = 500, importance = TRUE, nodesize = n)
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
      probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # store with repeat (r) and fold (f) index
      key <- paste0("nodesize_", n, "_Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
    }
  }
}

### calculate performance statistics
# create vectors to store metrics
node_number <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract mtry name
  ns <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+") 
  
  node_number <- c(node_number, ns)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc)
}

# combine metrics in a summary table
results_df <- data.frame(node_number = node_number,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         auc = auc_vals)

# summary of performance metrics
metric_summary <- results_df %>%
  group_by(node_number) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            sd_bal_acc = sd(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            sd_f1 = sd(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            sd_sens = sd(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            sd_spec = sd(spec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            sd_auc = sd(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(node_number)

# add ntree column for plotting
metric_summary$nodesize <- str_extract(metric_summary$node_number, "(?<=nodesize_)\\d+")
metric_summary$nodesize <- factor(metric_summary$nodesize, levels = c("1", "5", "10", "20", "30"))

# balanced accuracy
ggplot(metric_summary, aes(x = nodesize, y = mean_bal_acc, fill = nodesize)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("1" = "skyblue", "5" = "deepskyblue3", "10" = "skyblue3", "20" = "deepskyblue4", "30" = "skyblue4")) +
  labs(title = "Performance versus node size", 
       y = "Balanced accuracy", x = "Node size") + 
  theme_minimal() + theme(legend.position = "none")


# f1 score
ggplot(metric_summary, aes(x = nodesize, y = mean_f1, fill = nodesize)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("1" = "skyblue", "5" = "deepskyblue3", "10" = "skyblue3", "20" = "deepskyblue4", "30" = "skyblue4")) +
  labs(title = "Performance versus node size", 
       y = "F1 score", x = "Node size") + 
  theme_minimal() + theme(legend.position = "none")

# sensitivity
ggplot(metric_summary, aes(x = nodesize, y = mean_sens, fill = nodesize)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("1" = "skyblue", "5" = "deepskyblue3", "10" = "skyblue3", "20" = "deepskyblue4", "30" = "skyblue4")) +
  labs(title = "Performance versus node size", 
       y = "Sensitivity", x = "Node size") + 
  theme_minimal() + theme(legend.position = "none")

# specificity
ggplot(metric_summary, aes(x = nodesize, y = mean_spec, fill = nodesize)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("1" = "skyblue", "5" = "deepskyblue3", "10" = "skyblue3", "20" = "deepskyblue4", "30" = "skyblue4")) +
  labs(title = "Performance versus node size", 
       y = "Specificity", x = "Node size") + 
  theme_minimal() + theme(legend.position = "none")

# auc
ggplot(metric_summary, aes(x = nodesize, y = mean_auc, fill = nodesize)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("1" = "skyblue", "5" = "deepskyblue3", "10" = "skyblue3", "20" = "deepskyblue4", "30" = "skyblue4")) +
  labs(title = "Performance versus node size", 
       y = "AUC", x = "Node size") + 
  theme_minimal() + theme(legend.position = "none")


##########################################################################################################
###   RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - HYPERPARAMETER TUNING - PARALLELIZATION   ###
##########################################################################################################

# data to be used in the model
str(metagen)

set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# convert label to a factor 
metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease"))

# set hyperparameter grids
feature_sizes <- c(50, 75, 100, 200) # number of features to select
mtry_values <- c(30, 60, 120) # number of features per split
nodesize_values <- c(1, 5, 10) # number of samples for split to occur
weight_grid <- list(equal = c(healthy = 1, disease = 1),
                    med.d = c(healthy = 1, disease = 2)) # class weight settings

# outer cross-validation folds
outer_control <- createMultiFolds(metagen$condition, k = 5, times = 10)
outer_folds <- names(outer_control) # extract fold names 
outer_results <- list() # creat list to store results from outer folds

# create and register parallel backend
num_cores <- parallel::detectCores() - 1 # number of available cores (minus 1)
cl <- makeCluster(num_cores) # create cluster
registerDoParallel(cl) # resister parallel backend to enable parallel %dopar% operations

# outer loop with parallelization
outer_results <- foreach(outer = outer_folds, .packages = c("randomForest", "pROC", "caret")) %dopar% {
  
  # training and testing data for outer fold
  train_idx <- outer_control[[outer]] # training index
  test_idx  <- setdiff(seq_len(nrow(metagen)), train_idx) # testing index
  train_df <- metagen[train_idx, ] # training data
  test_df  <- metagen[test_idx, ] # testing data
  
  # re-factor levels
  train_df$condition <- factor(train_df$condition, levels = c("healthy", "disease"))
  test_df$condition <- factor(test_df$condition, levels = c("healthy", "disease"))
  
  # initialize variables to keep track of best model
  best_auc <- -Inf
  best_params <- NULL
  
  # nested hyperparameter tuning loop
  for (cw_lab in names(weight_grid)) {
    for (fs in feature_sizes) {
      for (m in mtry_values) {
        if (m > fs) next # skip invalid mtry values
        for (ns in nodesize_values) {
          
          # for each combination of hyperparametrs: inner 5-fold cross-validation
          auc_scores <- numeric(5) # initialize vector to store auc values
          inner_folds <- createFolds(train_df$condition, k = 5, list = TRUE) # create 5 stratified folds on training data
          
          # inner fold training and evaluation
          # split data into inner training and testing
          for (idx in seq_along(inner_folds)) {
            inner_val_idx <- inner_folds[[idx]]
            inner_tr <- train_df[-inner_val_idx, ]
            inner_val <- train_df[inner_val_idx, ]
            
            # re-factor labels
            inner_tr$condition <- factor(inner_tr$condition, levels = c("healthy", "disease"))
            inner_val$condition <- factor(inner_val$condition, levels = c("healthy", "disease"))
            
            # feature selection inside inner training data set
            # extract feature importance by meanDecreaseAccuracy
            rf_full <- randomForest(x = inner_tr[, all_feat_cols],
                                    y = inner_tr$condition,
                                    ntree = 500,
                                    classwt = weight_grid[[cw_lab]],
                                    importance = TRUE)
            
            imp <- importance(rf_full)[, "MeanDecreaseAccuracy"]
            feats <- names(sort(imp, decreasing = TRUE))[1:fs]
            
            # train inner random forest model on selected features and current hyperparameters
            rf_i <- randomForest(x = inner_tr[, feats],
                                 y = inner_tr$condition,
                                 ntree = 500,
                                 classwt = weight_grid[[cw_lab]],
                                 mtry = m,
                                 nodesize = ns)
            
            # evaluate on inner test fold
            # predict probabilties for class "disease"
            pr <- predict(rf_i, inner_val[, feats], type = "prob")[, "disease"]
            
            # calculate ROC curve and extract and store auc
            roc_i <- roc(inner_val$condition, pr,
                         levels = c("healthy", "disease"),
                         direction = "<")
            auc_scores[idx] <- auc(roc_i)
          }
          
          # average inner auc (over 5 inner folds)
          # update hyperparameters and auc if better than previous best
          mean_auc <- mean(auc_scores)
          if (mean_auc > best_auc) {
            best_auc <- mean_auc
            best_params <- list(feature_size = fs,
                                classwt_label = cw_lab,
                                mtry = m,
                                nodesize = ns)
          }
        }
      }
    }
  }
  
  # train final model (on full training set) with best hyperparameters from inner loop
  fs <- best_params$feature_size
  cw <- weight_grid[[best_params$classwt_label]]
  
  rf_full <- randomForest(x = train_df[, all_feat_cols],
                          y = train_df$condition,
                          ntree = 500,
                          classwt = cw,
                          importance = TRUE)
  
  # select top features again
  feats <- names(sort(importance(rf_full)[, "MeanDecreaseAccuracy"], decreasing = TRUE))[1:fs]
  
  # train final random forest model on selected features and best hyperparameters
  rf_final <- randomForest(x = train_df[, feats],
                           y = train_df$condition,
                           ntree = 500,
                           classwt = cw,
                           mtry = best_params$mtry,
                           nodesize = best_params$nodesize)
  
  # evaluate final model on outer test set
  # predict probabilities on outer testset
  prt <- predict(rf_final, test_df[, feats], type = "prob")[, "disease"]
  
  # calculate ROC and auc
  roc_final <- roc(test_df$condition, prt,
                   levels = c("healthy", "disease"),
                   direction = "<")
  
  # generate confusion matrix 
  cm <- confusionMatrix(predict(rf_final, test_df[, feats]),
                        test_df$condition,
                        positive = "disease")
  
  # return results for this outer fold
  # each parallel worker returns list of best hyperparameters, test AUC and confusion matrix for the fold
  list(best_params = best_params,
       test_auc = auc(roc_final),
       confusion = cm)
}

stopCluster(cl) # shut down the parallel cluster
registerDoSEQ() # unregister the parallel backend


# summarize and aggregate results
auc <- sapply(outer_results, `[[`, "test_auc") # extract all test AUCs
param_list <- lapply(outer_results, `[[`, "best_params") # extract best hyperparameters from each fol
params_df <- bind_rows(param_list, .id = "repeat_fold") # combine hyperparmeters into one data.frame with their fold ids

# extract balanced accuracy from each outer fold
bal_acc <- sapply(outer_results, function(res) {
  cm <- res$confusion
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  mean(c(sens, spec))
})

# # combine hyperparameters, auc and balanced accuracy into data.frame
# results_df <- data.frame(param_summary, auc = auc, bal_acc = bal_acc)

# combine hyperparameters, auc and balanced accuracy into data.frame
results_df <- data.frame(params_df, auc = auc, bal_acc = bal_acc)


# hyperparameter combination with the highest mean auc
best_combo_auc <- results_df %>%
  group_by(feature_size, classwt_label, mtry, nodesize) %>%
  summarise(mean_auc = mean(auc), .groups = "drop") %>%
  arrange(desc(mean_auc))

# hyperparameter combination with the highest mean balanced accuracy
best_combo_bal_acc <- results_df %>%
  group_by(feature_size, classwt_label, mtry, nodesize) %>%
  summarise(mean_bal_acc = mean(bal_acc), .groups = "drop") %>%
  arrange(desc(mean_bal_acc))


################################################################################################################
###   RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - HYPERPARAMETER TUNING - BAYESIAN OPTIMIZATION   ###
################################################################################################################

# data to be used in the model
str(metagen)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# convert label to a factor 
metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease"))

# define the scoring function for the Bayesian optimizer
scoring_function <- function(feature_size, mtry, nodesize, classwt_label) {
  fs <- as.integer(feature_size)
  m <- as.integer(mtry)
  ns <- as.integer(nodesize)
  cw <- c(healthy = 1, disease = 1 + classwt_label)
  
  # skip invalid combinations
  if (m > fs || fs > length(all_feat_cols)) {
    return(list(Score = 0))
  }
  
  # perform 5-fold cross-validation
  inner_folds <- caret::createMultiFolds(metagen$condition, k = 5, times = 10)
  bal_acc_values <- numeric(length(inner_folds))
  
  for (i in seq_along(inner_folds)) {
    
    # feature selection on inner training split
    idx <- inner_folds[[i]]
    inner_tr <- metagen[idx, ]
    inner_val <- metagen[-idx, ]
    
    rf_full <- tryCatch({
      randomForest(inner_tr[, all_feat_cols],
                   inner_tr$condition,
                   ntree = 500,
                   classwt = cw,
                   importance = TRUE)
      }, error = function(e) return(NULL))
    
    if (is.null(rf_full)) return(list(Score = 0)) # handle model failure
    
    imp <- importance(rf_full)[, "MeanDecreaseAccuracy"]
    feats <- names(sort(imp, decreasing = TRUE))[1:fs]
    
    rf_mod <- tryCatch({
      randomForest(inner_tr[, feats],
                   inner_tr$condition,
                   ntree = 500,
                   classwt = cw,
                   mtry = m,
                   nodesize = ns)
      }, error = function(e) return(NULL))
    
    if (is.null(rf_mod)) return(list(Score = 0)) # handle model failure
    
    preds <- predict(rf_mod, inner_val[, feats])
    cm <- caret::confusionMatrix(preds, inner_val$condition, positive = "disease")
    
    bal_acc_values[i] <- mean(c(cm$byClass["Sensitivity"], cm$byClass["Specificity"]), na.rm = TRUE)
  }
  mean_bal_acc <- mean(bal_acc_values, na.rm = TRUE)
  return(list(Score = ifelse(is.na(mean_bal_acc), 0, mean_bal_acc)))
}

# define parameter bounds for the Bayesian optimizer
bounds <- list(feature_size = c(50L, 200L),
               mtry = c(30L, 120L),
               nodesize = c(1L, 10L),
               classwt_label = c(0L, 1L)) # dummy numeric mapping: 0 for equal, 1 for med.d

# set seed
set.seed(1234)

# initialize and run Bayesian optimization
# uses one training set and performs multiple evaluations across hyperparameter space via Bayesian optimization (don't need to run inner cv for each outer fold)

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, varlist = c("metagen", "all_feat_cols"))
clusterEvalQ(cl, {
  library(randomForest)
  library(caret)
})

optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 20,
                   iters.n = 50,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

stopCluster(cl) # shut down the parallel cluster
registerDoSEQ() # unregister the parallel backend

### output summary
optObj$stopStatus # stop status of optimization
head(optObj$scoreSummary[order(-Score)]) # score = balanced accuracy
getBestPars(optObj) # optimal combination of hyperparameters


### loop over each outer fold to evaluate true generalization performance of final model using best hyperparameters
# set best values for hyperparameters
fs_best <- 150 # feature size
m_best  <- 100 # mtry
ns_best <- 10 # node size
cw_best <- c(healthy = 1, disease = 2) # class weight settings

# create lists to store metrics
final_feature_importances <- list() # list to store final feature importances
final_performance_metrics <- list() # list to store final model performances
final_feature_frequencies <- list() # list to store final feature frequencies
final_auc <- numeric() # to store auc values

# set seed
set.seed(1234)

# 10 repeats of stratified 5-fold cross-validation (50 models)
outer_folds <- createMultiFolds(metagen$condition, k = 5, times = 10)

# loop over the folds
for (f in seq_along(outer_folds)) {
  cat("Outer Fold:", f, "\n")
  
  # splits the dataset into training and testing sets
  train_idx <- outer_folds[[f]] # train indices
  train_data <- metagen[train_idx, ] # training data
  test_data  <- metagen[-train_idx, ] # testing data
  
  # train full random forest model only for feature selection
  rf_full <- randomForest(train_data[, all_feat_cols],
                          train_data$condition,
                          ntree = 500,
                          classwt = cw_best,
                          importance = TRUE)
  
  # select top fs_best features by MeanDecreaseAccuracy
  imp_vals <- importance(rf_full)[, "MeanDecreaseAccuracy"]
  top_feats <- names(sort(imp_vals, decreasing = TRUE))[1:fs_best]
  
  # train final random forest model with selected features
  rf_final <- randomForest(train_data[, top_feats],
                           train_data$condition,
                           ntree = 500,
                           classwt = cw_best,
                           mtry = m_best,
                           nodesize = ns_best,
                           importance = TRUE)
  
  # evaluate on test set
  preds <- predict(rf_final, test_data[, top_feats])
  probs <- predict(rf_final, test_data[, top_feats], type = "prob")
  
  # generate confusion matrix
  cm <- confusionMatrix(preds, test_data$condition, positive = "disease")
  
  # feature frequency from split variables
  split_vars <- unlist(lapply(1:rf_final$ntree, function(t) {
    tree <- getTree(rf_final, k = t, labelVar = TRUE)
    as.character(tree$`split var`[tree$`split var` != "<leaf>"])
  }))
  split_counts <- table(split_vars) # count the occurrences of each feature
  feature_freq_df <- as.data.frame(split_counts)
  colnames(feature_freq_df) <- c("Feature", "Count")
  
  # auc values
  roc_obj <- pROC::roc(response = test_data$condition,
                       predictor = probs[, "disease"],
                       levels = c("healthy", "disease"),
                       direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # store with repeat/fold (f) index
  key <- paste0("Model_", f)
  final_feature_importances[[key]] <- importance(rf_final) # store feature importances
  final_performance_metrics[[key]] <- cm # store confusion matrices
  final_feature_frequencies[[f]] <- feature_freq_df # store feature frequencies
  final_auc <- c(final_auc, auc_val) # store auc values
}

### calculate performance statistics
# create vectors to store metrics from confusion matrices
final_bal_acc <- final_f1 <- final_sens <- final_spec <- numeric()

# loop over confusion matrices
for (cm in final_performance_metrics) {
  final_bal_acc <- c(final_bal_acc, cm$byClass["Balanced Accuracy"])
  final_f1 <- c(final_f1, cm$byClass["F1"])
  final_sens <- c(final_sens, cm$byClass["Sensitivity"])
  final_spec <- c(final_spec, cm$byClass["Specificity"])
}

# summary dat.frame of performance metrics
final_metric_summary <- data.frame(mean_bal_acc = mean(final_bal_acc, na.rm = TRUE),
                                   sd_bal_acc = sd(final_bal_acc, na.rm = TRUE),
                                   mean_f1 = mean(final_f1, na.rm = TRUE),
                                   sd_f1 = sd(final_f1, na.rm = TRUE),
                                   mean_sens = mean(final_sens, na.rm = TRUE),
                                   sd_sens = sd(final_sens, na.rm = TRUE),
                                   mean_spec = mean(final_spec, na.rm = TRUE),
                                   sd_spec = sd(final_spec, na.rm = TRUE),
                                   mean_auc = mean(final_auc, na.rm = TRUE),
                                   sd_auc = sd(final_auc, na.rm = TRUE),
                                   min_auc = min(final_auc, na.rm = TRUE),
                                   max_auc = max(final_auc, na.rm = TRUE))
final_metric_summary

### combine all feature importances
final_all_importances <- do.call(rbind, lapply(names(final_feature_importances), function(name) {
  df <- as.data.frame(final_feature_importances[[name]])
  df$Feature <- rownames(df)
  df$Model <- name
  return(df)
}))

# calculate average MeanDecreaseAccuracy across all 50 models
final_mean_importance <- final_all_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(final_mean_importance, 20)

# plot species with highest MeanDecreaseAccuracy
ggplot(head(final_mean_importance, 20), 
       aes(x = reorder(Feature, mean_MeanDecreaseAccuracy), y = mean_MeanDecreaseAccuracy)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() + 
  labs(title = "Top 20 most important features by meanDecreaseAccuracy",
       x = "Feature", y = "MeanDecreaseAccuracy")

# plot abundance of top species
ggplot(metagen, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "Abundance", y = "Density of samples", fill = "Condition") +
  theme_minimal()

# plot abundance of top species
ggplot(metagen, aes(x = Acutalibacter_muris)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Acutalibacter muris",
       x = "Abundance", y = "Density of samples", fill = "Condition") +
  theme_minimal()

### combine all feature frequencies into one data.frame
final_all_splits <- bind_rows(final_feature_frequencies, .id = "Model")

# summarize how often each feature appeared and how frequently it was used
final_feature_split_summary <- final_all_splits %>%
  group_by(Feature) %>%
  summarise(total_count = sum(Count, na.rm = TRUE),
            mean_count = mean(Count, na.rm = TRUE),
            n_models = n()) %>%
  arrange(desc(total_count))
head(final_feature_split_summary, 20)

# calculate relative feature frequency selection
final_feature_split_summary <- final_feature_split_summary %>%
  mutate(prop_models = n_models/length(final_feature_frequencies),
         avg_per_tree = total_count/(length(final_feature_frequencies) * rf_final$ntree)) 

# plot total number of models where feature was used (top 30 features by frequency)
ggplot(head(final_feature_split_summary, 30), aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features in splits",
       x = "Feature", y = "Number of models")


# plot average number of models where feature was used (top 30 features by frequency)
ggplot(head(final_feature_split_summary, 30), aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features in splits",
       x = "Feature", y = "Frequency of selection")


### final model training - full dataset with best hyperparameters and best feature subset (fs_best)
# feature selection
rf_full_final <- randomForest(x = metagen[, all_feat_cols],
                              y = metagen$condition,
                              ntree = 500,
                              classwt = cw_best,
                              importance = TRUE)

# select features
imp_final <- importance(rf_full_final)[, "MeanDecreaseAccuracy"]
final_feats <- names(sort(imp_final, decreasing = TRUE))[1:fs_best]

# final model using selected features and optimal parameters
rf_final_model <- randomForest(x = metagen[, final_feats],
                               y = metagen$condition,
                               ntree = 500,
                               mtry = m_best,
                               nodesize = ns_best,
                               classwt = cw_best)

### save final model 
saveRDS(rf_final_model, file = "rf_final_model.rds")


################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - BORUTA  ###
################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies
boruta_selected_features <- list() # list to store Boruta selected features per fold

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(metagen$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # split dataset into training and testing sets for current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- metagen[test_idx, ] # testing data (fold f)
    
    # identify relevant features in training data with Boruta 
    boruta_result <- Boruta(x = train_data[, all_feat_cols], 
                            y = as.factor(train_data$condition),
                            doTrace = 0)
    
    # resolve tentative features
    boruta_final <- TentativeRoughFix(boruta_result)
    
    # select confirmed Boruta features only
    confirmed_feats <- getSelectedAttributes(boruta_final, withTentative = FALSE)
    
    # use all features if there are no Boruta confirmed features
    if (length(confirmed_feats) == 0) {
      warning(paste("No Boruta confirmed features in Repeat", r, "Fold", f))
      confirmed_feats <- all_feat_cols
    }
    
    # train random forest model using Boruta selected features
    rf_model <- randomForest(x = train_data[, confirmed_feats, drop = FALSE], 
                             y = as.factor(train_data$condition), 
                             ntree = 500, importance = TRUE) 
    
    # evaluate on test set (using on Boruta selected features)
    predictions <- predict(rf_model, newdata = test_data[, confirmed_feats, drop = FALSE])
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrix
    cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
    
    # store with repeat (r) and fold (f) index
    key <- paste0("Repeat_", r, "_Fold_", f) 
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- cm # store performance metrics
    feature_importances[[key]] <- importance(rf_model) # store feature importances
    boruta_selected_features[[key]] <- confirmed_feats # store Botuta selected features
  }
}

### feature frequency and importance analyses only use Boruta-selected features

### calculate feature frequencies
all_splits <- bind_rows(feature_frequencies, .id = "Repeat_Fold") # combine frequencies into a single data.frame
colnames(all_splits) <- c("Repeat_Fold", "Feature", "Count") # rename columns

# summarize total and average counts
feature_split_summary <- all_splits %>%
  group_by(Feature) %>%
  summarise(total_count = sum(Count, na.rm = TRUE),
            mean_count = mean(Count, na.rm = TRUE),
            n_models = n()) %>%
  arrange(desc(total_count))
head(feature_split_summary, 20)

# calculate relative frequency of feature selection
feature_split_summary <- feature_split_summary %>%
  mutate(prop_models = n_models / length(feature_frequencies),
         avg_per_tree = total_count / (length(feature_frequencies) * rf_model$ntree))

# total number of models where feature was used at least once
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = n_models)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features (Boruta) RF",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features (Boruta) RF",
       x = "Feature", y = "Number of models")


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()

# extract metrics from the stored confusion matrices (50 repeats x 5 folds = 250 values)
for (cm in performance_metrics) {
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE))
metric_summary


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

### plot species with highest MeanDecreaseAccuracy
ggplot(metagen, aes(x = Atlantibacter_hermannii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Atlantibacter hermannii",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggboxplot(metagen, x = "condition", y = "Atlantibacter_hermannii",
          color = "condition", add = "jitter") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "Abundance of discriminative species",
       subtitle = "Atlantibacter hermannii",
       y = "Relative Abundance") +
  theme_minimal()

ggplot(metagen, aes(x = Filifactor_alocis)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Filifactor alocis",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggboxplot(metagen, x = "condition", y = "Filifactor_alocis",
          color = "condition", add = "jitter") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "Abundance of discriminative species",
       subtitle = "Filifactor alocis",
       y = "Relative Abundance") +
  theme_minimal()

ggplot(metagen, aes(x = Actinomyces_wuliandei)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative species",
       subtitle = "Actinomyces wuliandei",
       x = "Abundance", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggboxplot(metagen, x = "condition", y = "Actinomyces_wuliandei",
          color = "condition", add = "jitter") +
  stat_compare_means(method = "wilcox.test") +
  labs(title = "Abundance of discriminative species",
       subtitle = "Actinomyces wuliandei",
       y = "Relative Abundance") +
  theme_minimal()


### Boruta feature selection frequency (across repeats and folds)
# how often feature was selected across models (50 repeats x 5 folds = 250)
boruta_feature_counts <- table(unlist(boruta_selected_features))
boruta_feature_summary <- as.data.frame(boruta_feature_counts)
colnames(boruta_feature_summary) <- c("feature", "times_selected")
boruta_feature_summary$prop_selected <- boruta_feature_summary$times_selected/length(boruta_selected_features) # proportion selected
boruta_feature_summary <- boruta_feature_summary[order(-boruta_feature_summary$times_selected), ]
head(boruta_feature_summary, 20)

### plot feature selection stability
ggplot(boruta_feature_summary[1:20, ], aes(x = reorder(feature, times_selected), y = times_selected)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 20 most frequently selected features by Boruta",
       x = "Feature", y = "Times selected across 250 models")

# features never selected (never considered relevant by Boruta)
never_selected <- setdiff(all_feat_cols, unique(unlist(boruta_selected_features)))


############################ SPLIT DEPTH OF FEATURES ###########################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# column names for features to be included in model (full predictor set)
all_feat_cols <- setdiff(colnames(metagen), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies
boruta_selected_features <- list() # list to store Boruta selected features per fold
feature_split_depths <- list() # list to store split depths

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(metagen$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # split dataset into training and testing sets for current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- metagen[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- metagen[test_idx, ] # testing data (fold f)
    
    # identify relevant features in training data with Boruta 
    boruta_result <- Boruta(x = train_data[, all_feat_cols], 
                            y = as.factor(train_data$condition),
                            doTrace = 0)
    
    # resolve tentative features
    boruta_final <- TentativeRoughFix(boruta_result)
    
    # select confirmed Boruta features only
    confirmed_feats <- getSelectedAttributes(boruta_final, withTentative = FALSE)
    
    # use all features if there are no Boruta confirmed features
    if (length(confirmed_feats) == 0) {
      warning(paste("No Boruta confirmed features in Repeat", r, "Fold", f))
      confirmed_feats <- all_feat_cols
    }
    
    # train random forest model using Boruta selected features
    rf_model <- randomForest(x = train_data[, confirmed_feats, drop = FALSE], 
                             y = as.factor(train_data$condition), 
                             ntree = 500, importance = TRUE) 
    
    # evaluate on test set (using on Boruta selected features)
    predictions <- predict(rf_model, newdata = test_data[, confirmed_feats, drop = FALSE])
    
    # count how often each feature is used in the trees
    tree_split_vars <- unlist(lapply(1:rf_model$ntree, function(t) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      as.character(tree$`split var`[tree$`split var` != "<leaf>"])
    }))
    
    # collect feature split depths
    split_depths <- data.frame()
    
    for (t in 1:rf_model$ntree) {
      tree <- getTree(rf_model, k = t, labelVar = TRUE)
      n_nodes <- nrow(tree)
      node_depths <- rep(NA, n_nodes)
      node_depths[1] <- 0  # root node
      
      # calculate depth of each node
      for (n in 2:n_nodes) {
        parent <- which(tree$`left daughter` == n | tree$`right daughter` == n)
        if (length(parent) == 1) {
          node_depths[n] <- node_depths[parent] + 1
        }
      }
      
      # record feature and depth
      for (n in 1:n_nodes) {
        feature <- as.character(tree$`split var`[n])
        if (!is.na(feature) && feature != "<leaf>") {
          split_depths <- rbind(split_depths, data.frame(
            Repeat = r,
            Fold = f,
            Tree = t,
            Feature = feature,
            Depth = node_depths[n]
          ))
        }
      }
    }
    
    # count the occurrences of each feature
    split_counts <- table(tree_split_vars)
    
    # generate confusion matrix
    cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
    
    # store with repeat (r) and fold (f) index
    key <- paste0("Repeat_", r, "_Fold_", f) 
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- cm # store performance metrics
    feature_importances[[key]] <- importance(rf_model) # store feature importances
    boruta_selected_features[[key]] <- confirmed_feats # store Botuta selected features
    feature_split_depths[[key]] <- split_depths # store split depths of features
  }
}

# combine all depth data
all_depths <- do.call(rbind, feature_split_depths)

# subset all depths for top 20 important features when using Boruta
all_depths_Boruta <- all_depths[all_depths$Feature %in% Boruta_important, ] 
all_depths_Boruta$importance <- "Boruta"

# subset all depths for top 20 important features when not using Boruta
all_depths_noBoruta <- all_depths[all_depths$Feature %in% noBoruta_important, ]
all_depths_noBoruta$importance <- "noBoruta"

# combine data.frames to plot
all_depths_importance <- rbind(all_depths_Boruta, all_depths_noBoruta)

# violin plot to see full depth distributions
ggplot(all_depths_importance, aes(x = importance, y = Depth, fill = importance)) +
  geom_violin(trim = FALSE) +     
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  labs(title = "Split depth distribution by feature category",
       y = "Split Depth",
       x = "Feature category") +
  theme_minimal()


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] future_1.58.0                 ParBayesianOptimization_1.2.6 doParallel_1.0.17            
# [4] iterators_1.0.14              foreach_1.5.2                 ggpubr_0.6.0                 
# [7] Boruta_9.0.0                  e1071_1.7-16                  pROC_1.18.5                  
# [10] caret_7.0-1                   lattice_0.22-7                randomForest_4.7-1.2         
# [13] lubridate_1.9.4               forcats_1.0.0                 stringr_1.5.1                
# [16] dplyr_1.1.4                   purrr_1.0.4                   readr_2.1.5                  
# [19] tidyr_1.3.1                   tibble_3.3.0                  tidyverse_2.0.0              
# [22] cowplot_1.1.3                 ggplot2_3.5.2                
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1     timeDate_4041.110    farver_2.1.2         digest_0.6.37       
# [5] rpart_4.1.24         timechange_0.3.0     lifecycle_1.0.4      survival_3.8-3      
# [9] dbscan_1.2.2         magrittr_2.0.3       compiler_4.5.0       rlang_1.1.6         
# [13] tools_4.5.0          data.table_1.17.4    ggsignif_0.6.4       plyr_1.8.9          
# [17] RColorBrewer_1.1-3   abind_1.4-8          withr_3.0.2          nnet_7.3-20         
# [21] grid_4.5.0           stats4_4.5.0         globals_0.18.0       scales_1.4.0        
# [25] MASS_7.3-65          dichromat_2.0-0.1    cli_3.6.5            crayon_1.5.3        
# [29] generics_0.1.4       rstudioapi_0.17.1    future.apply_1.20.0  reshape2_1.4.4      
# [33] tzdb_0.5.0           proxy_0.4-27         splines_4.5.0        vctrs_0.6.5         
# [37] hardhat_1.4.1        Matrix_1.7-3         carData_3.0-5        car_3.1-3           
# [41] hms_1.1.3            rstatix_0.7.2        Formula_1.2-5        listenv_0.9.1       
# [45] DiceKriging_1.6.0    gower_1.0.2          recipes_1.3.1        glue_1.8.0          
# [49] parallelly_1.45.0    codetools_0.2-20     stringi_1.8.7        gtable_0.3.6        
# [53] pillar_1.10.2        ipred_0.9-15         lava_1.8.1           R6_2.6.1            
# [57] lhs_1.2.0            backports_1.5.0      broom_1.0.8          class_7.3-23        
# [61] Rcpp_1.0.14          nlme_3.1-168         prodlim_2025.04.28   ModelMetrics_1.2.2.2
# [65] pkgconfig_2.0.3    

