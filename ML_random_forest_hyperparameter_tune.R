# Random Forest Workflow - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(cowplot)
library(tidyverse)
library(randomForest)
library(caret)
library(pROC)
library(e1071)
library(stringr)

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
dim(feat) # 2302   70
presence_threshold <- 0.10 * ncol(feat) # needs to be in at least 7 samples
features_to_keep <- rowSums(feat != 0) >= presence_threshold # determine if the feature is present >=10% of samples
feat_filt <- feat[features_to_keep, ] # subset the feature table to only include features present in at least 10% of samples
dim(feat_filt) # 935  70

# convert feature table to relative abundances
feat_filt <- apply(feat_filt, 2, function(x) x / sum(x))
feat_filt <- as.data.frame(feat_filt)

### log.unit normalization
log_n0 <- 1e-6 # pseudocount
n_p <- 2 # L2 norm 
# add pseudocount and perform log-transform
feat_log <- log(feat + log_n0)
# perform row-wise L2 normalization
row_norms <- sqrt(rowSums(feat_log^n_p))
feat_norm <- sweep(feat_log, 1, row_norms, FUN = "/")

# rownames of metadata need to match the column names of the feature table
all(rownames(meta) == colnames(feat_norm))


### merge metadata and feature table
feat <- t(feat_norm) # transpose feature table
feat <- as.data.frame(feat) 
feat$sample_name <- rownames(feat) # add column sample_name

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
head(mean_importance, 10)

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
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$condition), positive = "disease")
      
      # calculate AUC
      roc_obj <- roc(response = test_data$condition,
                     predictor = probabilities[, "disease"],
                     levels = c("healthy", "disease"),
                     direction = "<")
      auc_value <- auc(roc_obj)
      
      # store with repeat (r) and fold (f) index
      key <- paste0("Repeat_", r, "_Fold_", f, "_N_", n_feat)
      performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
    }
  }
}

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
weight_grid <- list(equal = c(healthy = 1, disease = 1),
                    mild = c(healthy = 1, disease = 2),
                    med = c(healthy = 1, disease = 3),
                    high = c(healthy = 1, disease = 5))

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
metric_summary$weight <- factor(metric_summary$weight, levels = c("equal", "mild", "med", "high"))


# balanced accuracy
ggplot(metric_summary, aes(x = weight, y = mean_bal_acc, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("equal" = "skyblue", "mild" = "pink", "med" = "plum", "high" = "salmon1")) +
  labs(title = "Performance versus class weights", 
       y = "Balanced accuracy", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# f1 score
ggplot(metric_summary, aes(x = weight, y = mean_f1, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("equal" = "skyblue", "mild" = "pink", "med" = "plum", "high" = "salmon1")) +
  labs(title = "erformance versus class weights", 
       y = "F1 score", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# sensitivity
ggplot(metric_summary, aes(x = weight, y = mean_sens, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("equal" = "skyblue", "mild" = "pink", "med" = "plum", "high" = "salmon1")) +
  labs(title = "erformance versus class weights", 
       y = "Sensitivity", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# specificity
ggplot(metric_summary, aes(x = weight, y = mean_spec, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("equal" = "skyblue", "mild" = "pink", "med" = "plum", "high" = "salmon1")) +
  labs(title = "erformance versus class weights", 
       y = "Specificity", x = "Class weights") + 
  theme_minimal() + theme(legend.position = "none")

# auc
ggplot(metric_summary, aes(x = weight, y = mean_auc, fill = weight)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("equal" = "skyblue", "mild" = "pink", "med" = "plum", "high" = "salmon1")) +
  labs(title = "erformance versus class weights", 
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

# ntree values to test
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] e1071_1.7-16         pROC_1.18.5          caret_7.0-1          lattice_0.22-7      
# [5] randomForest_4.7-1.2 lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
# [9] dplyr_1.1.4          purrr_1.0.4          readr_2.1.5          tidyr_1.3.1         
# [13] tibble_3.3.0         tidyverse_2.0.0      cowplot_1.1.3        ggplot2_3.5.2       
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6         recipes_1.3.1        tzdb_0.5.0           vctrs_0.6.5         
# [5] tools_4.5.0          generics_0.1.4       stats4_4.5.0         parallel_4.5.0      
# [9] proxy_0.4-27         pkgconfig_2.0.3      ModelMetrics_1.2.2.2 Matrix_1.7-3        
# [13] data.table_1.17.4    RColorBrewer_1.1-3   lifecycle_1.0.4      compiler_4.5.0      
# [17] farver_2.1.2         codetools_0.2-20     class_7.3-23         prodlim_2025.04.28  
# [21] pillar_1.10.2        MASS_7.3-65          gower_1.0.2          iterators_1.0.14    
# [25] rpart_4.1.24         foreach_1.5.2        nlme_3.1-168         parallelly_1.45.0   
# [29] lava_1.8.1           tidyselect_1.2.1     digest_0.6.37        stringi_1.8.7       
# [33] future_1.58.0        reshape2_1.4.4       listenv_0.9.1        labeling_0.4.3      
# [37] splines_4.5.0        grid_4.5.0           cli_3.6.5            magrittr_2.0.3      
# [41] utf8_1.2.6           dichromat_2.0-0.1    survival_3.8-3       future.apply_1.20.0 
# [45] withr_3.0.2          scales_1.4.0         timechange_0.3.0     globals_0.18.0      
# [49] nnet_7.3-20          timeDate_4041.110    hms_1.1.3            hardhat_1.4.1       
# [53] rlang_1.1.6          Rcpp_1.0.14          glue_1.8.0           ipred_0.9-15        
# [57] rstudioapi_0.17.1    R6_2.6.1             plyr_1.8.9  

