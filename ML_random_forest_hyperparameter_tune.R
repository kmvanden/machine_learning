# Random Forest Workflow - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(cowplot)
library(tidyverse)
library(compositions)
library(randomForest)
library(caret)
library(pROC)
library(Boruta)
library(ParBayesianOptimization)
library(doParallel)
library(fastshap)
library(viridis)
library(ggrepel)


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
feat_clr <- clr(feat_rel_abund + 1e-6)
feat <- as.data.frame(feat_clr)


# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat))
feat$sample_name <- rownames(feat) # add column sample_name
meta <- meta[, -3] # remove unnecessary meta column


### merge metadata and feature table
metagen <- merge(meta, feat, by = "sample_name", all.x = TRUE)
metagen <- metagen[,-1] # remove sample_name

# make sure names are syntactically valid 
colnames(metagen) <- make.names(colnames(metagen))


###############################################################################
###   BASELINE RANDOM FOREST MODEL - 5-FOLD CROSS-VALIDATION + 50 REPEATS   ###
###############################################################################

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
    predictions <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "response") # predicted class labels for cm
    probabilities <- predict(rf_model, newdata = test_data[, all_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
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
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
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
  labs(title = "Top 30 most frequently selected features - models",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "TTop 30 most frequently selected features - splits",
       x = "Feature", y = "Number of models")


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
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# extract metrics from the stored confusion matrices
for (perf in performance_metrics) {
  cm <- perf$cm
  auc_val <- as.numeric(perf$auc[])
  
  # confusion matrix metrics
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc_val)
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE),
                             mean_auc = mean(auc_vals, na.rm = TRUE),
                             sd_auc = sd(auc_vals, na.rm = TRUE))
metric_summary


#####################################################################################
###   RANDOM FOREST - RECURSIVE FEATURE ELIMINATION USING MEAN DECREASE ACCURACY  ###
#####################################################################################

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
    test_data <- metagen[test_idx, ] # testing data (fold f)
    
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
  labs(title = "Top 30 most frequently selected features by RF-based RFE - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
ggplot(feature_split_summary[1:30, ], aes(x = reorder(Feature, total_count), y = avg_per_tree)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Top 30 most frequently selected features by RF-based RFE - split",
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


##################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - CLASSWT  ###
##################################################################

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


################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - NTREE  ###
################################################################

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


###############################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - MTRY  ###
###############################################################

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


###################################################################
###   RANDOM FOREST - OPTIMAL HYPERPARAMETER VALUES - NODESIZE  ###
###################################################################

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


#################################################################################################
###   BORUTA FEATURE SELECTION + RANDOM FOREST MODEL - 5-FOLD CROSS-VALIDATION + 50 REPEATS   ###
#################################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

n_repeats <- 25
boruta_metrics_list <- list() # list to store Boruta metrics

for (i in 1:n_repeats) {
  cat("Running Boruta:", i, "\n")
  bor <- Boruta(condition ~ ., data = metagen, maxRuns = 250, ntree = 500, pValue = 0.01, doTrace = 0)
  bor <- TentativeRoughFix(bor)
  
  # boruta metrics
  imp_stats <- attStats(bor)
  imp_stats$Feature <- rownames(imp_stats)
  boruta_metrics_list[[i]] <- imp_stats
}


# aggregate median importance
boruta_df <- do.call(rbind, boruta_metrics_list)

# keep only confirmed features
confirmed_boruta_df <- boruta_df %>%
  filter(decision == "Confirmed") %>%
  group_by(Feature) %>%
  summarise(mean_medianImp = mean(medianImp, na.rm = TRUE),
            sd_medianImp = sd(medianImp, na.rm = TRUE),
            count = n()) %>%
  arrange(desc(mean_medianImp))

# plot stable features ranked by average median importance
ggplot(confirmed_boruta_df, aes(x = reorder(Feature, mean_medianImp), y = mean_medianImp)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Average median importance of confirmed features",
       x = "Feature", y = "Mean median importance")

# add proportion and sort by count
confirmed_boruta_df$proportion <- confirmed_boruta_df$count/n_repeats 
confirmed_boruta_df <- confirmed_boruta_df %>% arrange(desc(count))

# plot features selection frequency
ggplot(confirmed_boruta_df, aes(x = reorder(Feature, count), y = count)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Feature selection frequency across Boruta repeats",
       x = "Feature", y = paste("Number of runs selected (out of", n_repeats, ")"))  

# plot by proportion
ggplot(confirmed_boruta_df, aes(x = reorder(Feature, proportion), y = proportion)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Proportion of Boruta repeats selecting each feature",
       x = "Feature", y = "Proportion of Runs")  


# # extract confirmed features
# confirmed_feats <- getSelectedAttributes(boruta_resolved, withTentative = FALSE)

# keep only stable features
stable_feats <- confirmed_boruta_df %>% 
  filter(proportion >= 0.4) # model performance plateaus at 0.4
selected_features <- stable_feats$Feature

# subset metagen to just confirmed features
boruta_df <- metagen[, c("condition", selected_features)]


# data to be used in the model
str(boruta_df)

# set seed
set.seed(1234)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(boruta_df$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- boruta_df[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- boruta_df[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by subset_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                             y = as.factor(train_data$condition), 
                             ntree = 500, importance = TRUE) 
    
    # evaluate on test set
    predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
    probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
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
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
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
  labs(title = "Top 30 most frequently selected features - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
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

# group importance metrics by feature and sort by overall importance
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(mean_importance, 20)


### plot features with highest MeanDecreaseAccuracy
ggplot(boruta_df, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Anaerobutyricum_hallii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Anaerobutyricum hallii",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Clostridium_sp._M62.1)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Clostridium sp. M62.1",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# extract metrics from the stored confusion matrices
for (perf in performance_metrics) {
  cm <- perf$cm
  auc_val <- as.numeric(perf$auc[])
  
  # confusion matrix metrics
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc_val)
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE),
                             mean_auc = mean(auc_vals, na.rm = TRUE),
                             sd_auc = sd(auc_vals, na.rm = TRUE))
metric_summary


################################################################################################################
########   RANDOM FOREST - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF BAYES OPT - AUC  #######
################################################################################################################

# data to be used in the model
str(boruta_df)

# set seed
set.seed(1234)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# create list of class weight settings
weight_grid <- list(high.h = c(healthy = 3, disease = 1),
                    med.h = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    med.d = c(healthy = 1, disease = 2))

# define the set of categorical labels with numeric indices (classwt_label is a numeric index from 1 to 4)
label_keys <- c("high.h", "med.h", "equal", "med.d")

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
    folds <- caret::createFolds(boruta_df$condition, k = 5, list = TRUE, returnTrain = FALSE)
    fold_aucs <- numeric(length(folds))
    
    # loop through the folds
    for (f in seq_along(folds)) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- boruta_df[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- boruta_df[test_idx, ] # testing data (fold f)
      
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
bounds <- list(mtry = c(1L, 9L),
               ntree = c(100L, 1000L),
               nodesize = c(1L, 10L),
               classwt_label = c(1L, length(label_keys))) # numeric range for classwt_label

# resister back end
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 20,
                   iters.n = 50,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ()

# view results
stopstatus_auc = optObj$stopStatus
scoresum_auc = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
bestparams_auc = getBestPars(optObj)


##################################################################################################################
########   RANDOM FOREST - EVALUATION OF MODEL WITH BEST HYPERPARAMETERS FROM BAYESIAN OPTIMISATION - AUC  #######
##################################################################################################################

# data to be used in the model
str(boruta_df)

# set seed
set.seed(1234)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# optimal hyperparameter values
mtry_opt <- bestparams_auc$mtry
ntree_opt <- bestparams_auc$ntree
nodesize_opt <- bestparams_auc$nodesize

# # create list of class weight settings
# weight_grid <- list(high.h = c(healthy = 3, disease = 1),
#                     med.h = c(healthy = 2, disease = 1),
#                     equal = c(healthy = 1, disease = 1),
#                     med.d = c(healthy = 1, disease = 2))
# 
# # define the set of categorical labels with numeric indices
# label_keys <- c("high.h", "med.h", "equal", "med.d")
classwt_opt <- weight_grid[[label_keys[bestparams_auc$classwt_label]]]

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(boruta_df$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- boruta_df[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- boruta_df[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by subset_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                             y = as.factor(train_data$condition), 
                             mtry = mtry_opt,
                             ntree = ntree_opt,
                             nodesize = nodesize_opt,
                             classwt = classwt_opt,
                             importance = TRUE) 
    
    # evaluate on test set
    predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
    probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
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
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
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
  labs(title = "Top 30 most frequently selected features - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
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

# group importance metrics by feature and sort by overall importance
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(mean_importance, 20)


### plot features with highest MeanDecreaseAccuracy
ggplot(boruta_df, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Anaerobutyricum_hallii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Anaerobutyricum hallii",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Clostridium_sp._M62.1)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Clostridium sp. M62.1",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# extract metrics from the stored confusion matrices
for (perf in performance_metrics) {
  cm <- perf$cm
  auc_val <- as.numeric(perf$auc[])
  
  # confusion matrix metrics
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc_val)
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE),
                             mean_auc = mean(auc_vals, na.rm = TRUE),
                             sd_auc = sd(auc_vals, na.rm = TRUE))
metric_summary


###########################################################################################
########   RANDOM FOREST MODEL - TRAIN FINAL MODEL WITH BEST HYPERPARAMETERS - AUC  #######
###########################################################################################

# data to be used in the model
str(boruta_df)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# best hyperparameter values determined by Bayesian optimization
# optimal hyperparameter values
mtry_opt <- bestparams_auc$mtry
ntree_opt <- bestparams_auc$ntree
nodesize_opt <- bestparams_auc$nodesize

# # create list of class weight settings
# weight_grid <- list(high.h = c(healthy = 3, disease = 1),
#                     med.h = c(healthy = 2, disease = 1),
#                     equal = c(healthy = 1, disease = 1),
#                     med.d = c(healthy = 1, disease = 2))
# 
# # define the set of categorical labels with numeric indices
# label_keys <- c("high.h", "med.h", "equal", "med.d")
classwt_opt <- weight_grid[[label_keys[bestparams_auc$classwt_label]]]

# train the model on the full dataset

set.seed(1234)
final_model_auc <- randomForest(x = boruta_df[, subset_feat_cols], 
                                y = as.factor(boruta_df$condition), 
                                mtry = mtry_opt,
                                ntree = ntree_opt,
                                nodesize = nodesize_opt,
                                classwt = classwt_opt,
                                importance = TRUE)

# # save final model
# saveRDS(final_model_auc, file = "final_rf_model_metagen_auc.rds")


#################################################################################
###   XGBOOST LOGLOSS MODEL - SHAP VALUES - DEPENDENCE AND INTERACTION PLOTS  ###
#################################################################################

# function to get predicted disease probability
pred_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "disease"]
}

# compute fast SHAP values
set.seed(1234) # set seed

shap_values <- fastshap::explain(object = final_model_auc,
                                 X = boruta_df[, subset_feat_cols],
                                 pred_wrapper = pred_fun,
                                 nsim = 100, # number of Monte Carlo permutations
                                 adjust = TRUE) # centered SHAP values (baseline + sum of SHAPs = prediction)
shap_values <- as.data.frame(shap_values)     


### summarize SHAP values
shap_long <- shap_values %>%
  mutate(ID = row_number()) %>%
  pivot_longer(-ID, names_to = "feature", values_to = "value")

# add rfvalues (abundance) to shap_long
rfvalue_long <- boruta_df[, subset_feat_cols] %>%
  mutate(ID = row_number()) %>%
  pivot_longer(cols = -ID, names_to = "feature", values_to = "rfvalue")
shap_long <- shap_long %>%
  left_join(rfvalue_long, by = c("ID", "feature"))

# mean shap values
mean_shap_df <- shap_long %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_shap = mean(value, na.rm = TRUE)) 

# join mean abs values and mean values
shap_summary <- shap_long %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE)) %>%
  dplyr::arrange(desc(mean_abs_shap)) %>%
  dplyr::left_join(mean_shap_df, by = "feature")

# plot mean absolute SHAP value per feature
ggplot(shap_summary, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Mean absolute SHAP value", 
       x = "Feature", y = "Mean absolute SHAP value")

# plot mean SHAP value per feature
ggplot(shap_summary, aes(x = reorder(feature, mean_shap), y = mean_shap)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Mean SHAP value", 
       x = "Feature", y = "Mean SHAP value")

# recreate shap.plot.summary from treeshap (beeswarm-style plot)
shap_plot <- shap_long %>%
  dplyr::group_by(feature) %>%
  dplyr::mutate(scaled_rfvalue = scale(rfvalue)[,1]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(feature = factor(feature, levels = rev(shap_summary$feature)))

ggplot(shap_plot, aes(x = value, y = feature, color = scaled_rfvalue)) +
  geom_jitter(height = 0.2, alpha = 0.7, size = 1.2) +
  scale_color_viridis_c(option = "plasma", direction = -1) + theme_minimal() +
  labs(title = "SHAP summary plot", x = "SHAP value (impact on model output)", 
       color = "Feature value")


### SHAP dependence plots (how SHAP values for a given feature vary as the input values for the feature vary)
# wide table of log-transformed relative abundance
rfvalue_wide <- boruta_df[, subset_feat_cols] %>%
  mutate(ID = row_number())

feature_name <- "Lachnoclostridium_sp._YL32"
shap_dep <- shap_long %>% filter(feature == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "Clostridium_sp._M62.1"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- log abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- log abun"))


feature_name <- "Petrimonas_mucosa"
shap_dep <- shap_long %>% filter(feature == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "Enterocloster_clostridioformis"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- log abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- log abun"))


### mean absolute SHAP value per class per feature
shap_values$condition <- boruta_df$condition
abs_shap_by_class <- shap_values %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

# plot mean absolute SHAP value per class per feature
ggplot(abs_shap_by_class, aes(x = reorder(name, mean_abs_shap), y = mean_abs_shap, fill = condition)) +
  geom_col(position = "dodge") + coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean absolute SHAP values",
       x = "Feature", y = "Mean abs SHAP", fill = "Condition")


### mean SHAP value per class per feature
shap_values$condition <- boruta_df$condition
abs_shap_by_class <- shap_values %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_shap = mean(value), .groups = "drop")

# plot mean SHAP value per class per feature
ggplot(abs_shap_by_class, aes(x = reorder(name, mean_shap), y = mean_shap, fill = condition)) +
  geom_col(position = "dodge") + coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean SHAP values",
       x = "Feature", y = "Mean SHAP", fill = "Condition")


### how the impact of a given feature on model prediction varies between healthy and disease samples
# histogram of SHAP values
ggplot(shap_values, aes(x = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - Lachnoclostridium_sp._YL32",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_values, aes(x = condition, y = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for Lachnoclostridium_sp._YL32 by condition",
       y = "SHAP value", x = "Condition")

# box plot of log-abundance
ggplot(boruta_df, aes(x = condition, y = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "Log abundance of Lachnoclostridium_sp._YL32 by condition",
       y = "Log Abundance", x = "Condition")


# histogram of SHAP values
ggplot(shap_values, aes(x = Petrimonas_mucosa, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - Petrimonas_mucosa",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_values, aes(x = condition, y = Petrimonas_mucosa, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for Petrimonas_mucosa by condition",
       y = "SHAP value", x = "Condition")

# box plot of log-abundance
ggplot(boruta_df, aes(x = condition, y = Petrimonas_mucosa, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "Log abundance of Petrimonas_mucosa by condition",
       y = "Log abundance", x = "Condition")


### SHAP values versus predicted probabilities
shap_values$ID <- rownames(shap_values)
prob <- as.data.frame(predict(final_model_auc, newdata = boruta_df[, subset_feat_cols], type = "prob"))
shap_values$pred_prob <- prob$disease

# Petrimonas_mucosa
pred_plot <- shap_values %>%
  select(ID, Petrimonas_mucosa, condition, pred_prob)

ggplot(pred_plot, aes(x = Petrimonas_mucosa, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Petrimonas_mucosa SHAP value versus prediction probability",
       x = "SHAP value for Petrimonas_mucosa", y = "Predicted probability (disease)")

# Lachnoclostridium_sp._YL32
pred_plot <- shap_values %>%
  select(ID, Lachnoclostridium_sp._YL32, condition, pred_prob)

ggplot(pred_plot, aes(x = Lachnoclostridium_sp._YL32, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Lachnoclostridium_sp._YL32 SHAP value versus prediction probability",
       x = "SHAP value for Lachnoclostridium_sp._YL32", y = "Predicted probability (disease)")


### plot SHAP values versus importance 
importance <- as.data.frame(importance(final_model_auc))
importance$feature <- rownames(importance)
shap_import_df <- shap_summary %>%
  left_join(importance, by = "feature")

# plot SHAP values verus MeanDecreaseAccuracy
ggplot(shap_import_df, aes(x = MeanDecreaseAccuracy, y = mean_abs_shap, label = feature)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_text_repel(size = 3) + theme_minimal() +
  labs(title = "Mean absolute SHAP versus mean decrease accuracy",
       x = "Mean decrease accuracy",
       y = "Mean absolute SHAP value")

# plot SHAP values verus MeanDecreaseGini
ggplot(shap_import_df, aes(x = MeanDecreaseGini, y = mean_abs_shap, label = feature)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_text_repel(size = 3) + theme_minimal() +
  labs(title = "Mean absolute SHAP versus mean decrease in Gini",
       x = "Mean decrease in Gini",
       y = "Mean absolute SHAP value")


### fastshap does not calculate second-order SHAP effects
# partial dependence (PDP) interaction plots (how much each feature interacts with others)
library(iml)
predictor <- Predictor$new(final_model_auc, data = boruta_df[, subset_feat_cols], y = boruta_df$condition, type = "prob")
interaction_strength <- Interaction$new(predictor)
plot(interaction_strength)


####################################################################################################################
########   RANDOM FOREST - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF BAYES OPT - BAL_ACC  #######
####################################################################################################################

# set seed
set.seed(1234)

# data to be used in the model
str(boruta_df)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# create list of class weight settings
weight_grid <- list(high.h = c(healthy = 3, disease = 1),
                    med.h = c(healthy = 2, disease = 1),
                    equal = c(healthy = 1, disease = 1),
                    med.d = c(healthy = 1, disease = 2))

# define the set of categorical labels with numeric indices (classwt_label is a numeric index from 1 to 4)
label_keys <- c("high.h", "med.h", "equal", "med.d")

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
  repeat_bal_acc <- numeric(repeats)
  
  # loop over each repeat
  for (r in 1:repeats) {
    
    #  create 5-folds for cross-validation (stratified on condition)
    folds <- caret::createFolds(boruta_df$condition, k = 5, list = TRUE, returnTrain = FALSE)
    fold_bal_acc <- numeric(length(folds))
    
    # loop through the folds
    for (f in seq_along(folds)) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- boruta_df[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- boruta_df[test_idx, ] # testing data (fold f)
      
      # train random forest model using full features to rank features
      rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                               y = as.factor(train_data$condition),
                               mtry = mtry,
                               ntree = ntree, 
                               nodesize = nodesize,
                               classwt = classwt,
                               importance = TRUE)
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      fold_bal_acc[f] <- cm$byClass["Balanced Accuracy"]
    }
    
    repeat_bal_acc[r] <- mean(fold_bal_acc, na.rm = TRUE)
  }
  
  mean_bal_acc <- mean(repeat_bal_acc, na.rm = TRUE) # average across repeats
  return(list(Score = mean_bal_acc)) # return mean auc values
}


# define parameter bounds
bounds <- list(mtry = c(1L, 9L),
               ntree = c(100L, 1000L),
               nodesize = c(1L, 10L),
               classwt_label = c(1L, length(label_keys))) # numeric range for classwt_label

# resister back end
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 20,
                   iters.n = 50,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ()

# view results
stopstatus_bal_acc = optObj$stopStatus
scoresum_bal_acc = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
bestparams_bal_acc = getBestPars(optObj)


######################################################################################################################
########   RANDOM FOREST - EVALUATION OF MODEL WITH BEST HYPERPARAMETERS FROM BAYESIAN OPTIMISATION - BAL_ACC  #######
######################################################################################################################

# data to be used in the model
str(boruta_df)

# set seed
set.seed(1234)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# optimal hyperparameter values
mtry_opt <- bestparams_bal_acc$mtry
ntree_opt <- bestparams_bal_acc$ntree
nodesize_opt <- bestparams_bal_acc$nodesize

# # create list of class weight settings
# weight_grid <- list(high.h = c(healthy = 3, disease = 1),
#                     med.h = c(healthy = 2, disease = 1),
#                     equal = c(healthy = 1, disease = 1),
#                     med.d = c(healthy = 1, disease = 2))
# 
# # define the set of categorical labels with numeric indices
# label_keys <- c("high.h", "med.h", "equal", "med.d")
classwt_opt <- weight_grid[[label_keys[bestparams_bal_acc$classwt_label]]]

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on condition)
  folds <- createFolds(boruta_df$condition, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- boruta_df[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- boruta_df[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by subset_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, subset_feat_cols], 
                             y = as.factor(train_data$condition), 
                             mtry = mtry_opt,
                             ntree = ntree_opt,
                             nodesize = nodesize_opt,
                             classwt = classwt_opt,
                             importance = TRUE) 
    
    # evaluate on test set
    predictions <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "response") # predicted class labels for cm
    probabilities <- predict(rf_model, newdata = test_data[, subset_feat_cols], type = "prob") # class probabilities (ROC/AUC)
    
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
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature frequencies
    performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
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
  labs(title = "Top 30 most frequently selected features - model",
       x = "Feature", y = "Number of models")

# average number of times feature was used in a split per tree (across all models) 
# 250 models (50 repeats x 5-fold CV) each with 500 trees (125,000 trees in total)
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

# group importance metrics by feature and sort by overall importance
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_disease = mean(disease, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseAccuracy))
head(mean_importance, 20)


### plot features with highest MeanDecreaseAccuracy
ggplot(boruta_df, aes(x = Lachnoclostridium_sp._YL32)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Lachnoclostridium sp. YL32",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Anaerobutyricum_hallii)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Anaerobutyricum hallii",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()

ggplot(boruta_df, aes(x = Clostridium_sp._M62.1)) +
  geom_density(aes(fill = condition), alpha = 0.5) +
  labs(title = "Abundance of discriminative features",
       subtitle = "Clostridium sp. M62.1",
       x = "log2(Abundance)", y = "Density of Samples", fill = "Condition") +
  theme_minimal()


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_vals <- numeric()

# extract metrics from the stored confusion matrices
for (perf in performance_metrics) {
  cm <- perf$cm
  auc_val <- as.numeric(perf$auc[])
  
  # confusion matrix metrics
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_vals <- c(auc_vals, auc_val)
}

# combine metrics in a summary table
metric_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE),
                             mean_auc = mean(auc_vals, na.rm = TRUE),
                             sd_auc = sd(auc_vals, na.rm = TRUE))
metric_summary


#####################################################################################
########   RANDOM FOREST MODEL - TRAIN FINAL MODEL WITH BEST HYPERPARAMETERS  #######
#####################################################################################

# data to be used in the model
str(boruta_df)

# set seed
set.seed(1234)

# column names for features to be included in model
subset_feat_cols <- setdiff(colnames(boruta_df), "condition")

# best hyperparameter values determined by Bayesian optimization
# optimal hyperparameter values
mtry_opt <- bestparams_bal_acc$mtry
ntree_opt <- bestparams_bal_acc$ntree
nodesize_opt <- bestparams_bal_acc$nodesize

# # create list of class weight settings
# weight_grid <- list(high.h = c(healthy = 3, disease = 1),
#                     med.h = c(healthy = 2, disease = 1),
#                     equal = c(healthy = 1, disease = 1),
#                     med.d = c(healthy = 1, disease = 2))
# 
# # define the set of categorical labels with numeric indices
# label_keys <- c("high.h", "med.h", "equal", "med.d")
classwt_opt <- weight_grid[[label_keys[bestparams_bal_acc$classwt_label]]]

# train the model on the full dataset

set.seed(1234)
final_model_bal_acc <- randomForest(x = boruta_df[, subset_feat_cols], 
                                    y = as.factor(boruta_df$condition), 
                                    mtry = mtry_opt,
                                    ntree = ntree_opt,
                                    nodesize = nodesize_opt,
                                    classwt = classwt_opt,
                                    importance = TRUE)

# # save final model
# saveRDS(final_model_bal_acc, file = "final_rf_model_metagen_bal_acc.rds")


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
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
#   [1] DiceKriging_1.6.0             doParallel_1.0.17            
# [3] iterators_1.0.14              foreach_1.5.2                
# [5] ParBayesianOptimization_1.2.6 Boruta_9.0.0                 
# [7] pROC_1.19.0.1                 caret_7.0-1                  
# [9] lattice_0.22-7                randomForest_4.7-1.2         
# [11] compositions_2.0-8            lubridate_1.9.4              
# [13] forcats_1.0.0                 stringr_1.5.1                
# [15] dplyr_1.1.4                   purrr_1.1.0                  
# [17] readr_2.1.5                   tidyr_1.3.1                  
# [19] tibble_3.3.0                  tidyverse_2.0.0              
# [21] cowplot_1.2.0                 ggplot2_3.5.2                
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1     timeDate_4041.110    farver_2.1.2         tensorA_0.36.2.1    
# [5] digest_0.6.37        rpart_4.1.24         timechange_0.3.0     lifecycle_1.0.4     
# [9] survival_3.8-3       magrittr_2.0.3       dbscan_1.2.2         compiler_4.5.0      
# [13] rlang_1.1.6          tools_4.5.0          utf8_1.2.6           data.table_1.17.8   
# [17] ggsignif_0.6.4       labeling_0.4.3       plyr_1.8.9           RColorBrewer_1.1-3  
# [21] abind_1.4-8          withr_3.0.2          nnet_7.3-20          grid_4.5.0          
# [25] stats4_4.5.0         ggpubr_0.6.1         e1071_1.7-16         future_1.67.0       
# [29] globals_0.18.0       scales_1.4.0         MASS_7.3-65          dichromat_2.0-0.1   
# [33] cli_3.6.5            crayon_1.5.3         generics_0.1.4       rstudioapi_0.17.1   
# [37] future.apply_1.20.0  robustbase_0.99-4-1  reshape2_1.4.4       tzdb_0.5.0          
# [41] proxy_0.4-27         bayesm_3.1-6         splines_4.5.0        vctrs_0.6.5         
# [45] hardhat_1.4.1        Matrix_1.7-3         carData_3.0-5        car_3.1-3           
# [49] hms_1.1.3            rstatix_0.7.2        Formula_1.2-5        listenv_0.9.1       
# [53] gower_1.0.2          recipes_1.3.1        glue_1.8.0           parallelly_1.45.1   
# [57] DEoptimR_1.1-4       codetools_0.2-20     stringi_1.8.7        gtable_0.3.6        
# [61] pillar_1.11.0        ipred_0.9-15         lava_1.8.1           R6_2.6.1            
# [65] lhs_1.2.0            backports_1.5.0      broom_1.0.9          class_7.3-23        
# [69] Rcpp_1.1.0           nlme_3.1-168         prodlim_2025.04.28   mgcv_1.9-3          
# [73] ranger_0.17.0        pkgconfig_2.0.3      ModelMetrics_1.2.2.2
