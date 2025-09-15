# XGBoost Workflow - Metabolomis - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(tidyverse)
library(readxl)
library(impute)
library(xgboost)
library(caret)
library(pROC)
library(MLmetrics)
library(Matrix)
library(ParBayesianOptimization)
library(doParallel)
library(foreach)
library(SHAPforxgboost)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
### metadata
meta <- read.table("metadata.txt", header = TRUE)

# subset meta to samples present in metabolomics data
sample_key <- read_excel("metabo_batch.xlsx", sheet = "Sample Meta Data") # sample name key
sample_key <- sample_key %>% mutate(CLIENT_SAMPLE_ID = sub("KMV_", "SMS_", CLIENT_SAMPLE_ID)) # rename sample ids
all(sample_key$CLIENT_SAMPLE_ID %in% meta$sample_id)
meta <- merge(meta, sample_key[, c("CLIENT_SAMPLE_ID", "PARENT_SAMPLE_NAME")], 
              by.x = "sample_id", by.y = "CLIENT_SAMPLE_ID", all.x = TRUE) %>% # add PARENT_SAMPLE_NAME to meta
  filter(!is.na(PARENT_SAMPLE_NAME)) # subset

# convert condition column into a factor
rownames(meta) <- meta$sample_name
meta$condition <- as.factor(meta$condition) 
table(meta$condition)


### feature table
metabo <- read_excel("metabo_batch.xlsx", sheet = "Batch-normalized Data") # batch-normalized feature table
all(metabo$PARENT_SAMPLE_NAME %in% meta$PARENT_SAMPLE_NAME)

# replace names in PARENT_SAMPLE_NAME column in metabo with names from sample_name (from meta)
metabo <- metabo %>%
  left_join(meta %>% dplyr::select(PARENT_SAMPLE_NAME, sample_name), by = "PARENT_SAMPLE_NAME") %>%
  mutate(PARENT_SAMPLE_NAME = sample_name) %>%
  dplyr::select(-sample_name) %>%
  dplyr::slice(match(meta$sample_name, PARENT_SAMPLE_NAME)) %>% # match order of samples in meta
  column_to_rownames(var = "PARENT_SAMPLE_NAME") # move sample id column to rownames
all(rownames(metabo) == rownames(meta))

# replace CHEM_ID (colnames of metabo) with names from CHEMICAL_NAME (from metabo_key)
metabo_key <- read_excel("metabo_batch.xlsx", sheet = "Chemical Annotation")
all(colnames(metabo) == metabo_key$CHEM_ID)
name_map <- metabo_key %>% dplyr::select(CHEM_ID, CHEMICAL_NAME) %>% deframe() 
colnames(metabo) <- name_map[colnames(metabo)] %>% ifelse(is.na(.), colnames(metabo), .)
all(colnames(metabo) == metabo_key$CHEMICAL_NAME)
colnames(metabo) <- make.names(colnames(metabo)) # make sure names are syntactically valid 

### feature filtering, zero imputation and log transformation 
# assess missingness (zeros listed as NAs in data.frame)
total_missing <- sum(is.na(metabo)) / prod(dim(metabo)) * 100 # overall missingness
missing_per_feature <- colMeans(is.na(metabo)) * 100 # missingness per feature
summary(missing_per_feature)
missing_per_sample <- rowMeans(is.na(metabo)) * 100 # missingness per sample
hist(missing_per_feature, breaks = 50, main = "Percent missing per metabolite", xlab = "Percent missing")

# remove features with greater than 30% missing values
metabo_filtered <- metabo[, missing_per_feature <= 30]

# impute remaining zeros with kNN
metabo_imputed <- impute.knn(as.matrix(t(metabo_filtered)))$data # impute.knn expects features as rows
metabo_imputed <- t(metabo_imputed)

# log2 transformation (with pseudocount)
metabo <- log2(metabo_imputed + 1e-6)


# rownames of metadata need to match the rownames of the feature table
metabo <- as.data.frame(metabo)
all(rownames(meta) == rownames(metabo))
metabo$sample_name <- rownames(metabo) # add column sample_name
meta <- meta[, -c(1,4)] # remove meta columns that were used to format data.frames

# merge metadata and feature table
metabo_df <- merge(meta, metabo, by = "sample_name", all.x = TRUE)
metabo_df <- metabo_df[,-1] # remove sample_name


### for xgboost analysis
# convert condition to numeric (0 = healthy, 1 = disease)
metabo_df$condition_numeric <- ifelse(metabo_df$condition == "disease", 1, 0)
# set factor labels (healthy = negative, disease = positive)
metabo_df$condition <- factor(metabo_df$condition, levels = c("healthy", "disease"))

# set predictor columns
all_feat_cols <- setdiff(colnames(metabo_df), c("condition", "condition_numeric"))


############################################################
###   BASELINE XGBOOST MODEL - DEFAULT HYPERPARAMETERS   ###
############################################################

# data to be used in the model
str(metabo_df)

# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation (250 models)
folds <- createMultiFolds(metabo_df$condition, k = 5, times = 50)

# create list to store performance metrics
xgb_metrics <- list() # list to store performance metrics
xgb_importances <- list() # list to store feature importances

# loop through the folds
for (key in names(folds)) {
  
  # splits the dataset into training and testing sets
  train_idx <- folds[[key]] # train indices
  train_data <- metabo_df[train_idx, ] # training data
  test_data <- metabo_df[-train_idx, ] # testing data
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), label = train_data$condition_numeric)
  dtest  <- xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), label = test_data$condition_numeric)
  
  # train XGBoost model (binary classification with default parameters)
  xgb_model <- xgboost(data = dtrain,
                       objective = "binary:logistic",
                       eval_metric = "logloss",
                       nrounds = 5000,
                       verbose = 0)
  
  # evaluate on test set
  preds_prob <- predict(xgb_model, dtest)
  preds_label <- ifelse(preds_prob > 0.5, "disease", "healthy")
  preds_label <- factor(preds_label, levels = c("healthy", "disease"))
  
  # calculate AUC
  auc_val <- auc(response = test_data$condition,
                 predictor = preds_prob,
                 levels = c("healthy", "disease"),
                 direction = "<")
  
  # generate confusion matrix
  cm <- confusionMatrix(preds_label, test_data$condition, positive = "disease")
  
  # store metrics
  xgb_metrics[[key]] <- list(cm = cm, auc = auc_val) # store confusion matrices
  imp <- xgb.importance(feature_names = all_feat_cols, model = xgb_model) # feature importances
  imp$Repeat_Fold <- key # indices for feature importances
  xgb_importances[[key]] <- imp # store feature importances
}


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_values <- numeric()

# loop through stored metrics to extract values
for (key in names(xgb_metrics)) {
  cm <- xgb_metrics[[key]]$cm
  auc_val <- xgb_metrics[[key]]$auc
  
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_values <- c(auc_values, auc_val)
}

# data.frame of performance metrics
xgb_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                          sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                          mean_f1 = mean(f1_score, na.rm = TRUE),
                          sd_f1 = sd(f1_score, na.rm = TRUE),
                          mean_sens = mean(sensitivity, na.rm = TRUE),
                          sd_sens = sd(sensitivity, na.rm = TRUE),
                          mean_spec = mean(specificity, na.rm = TRUE),
                          sd_spec = sd(specificity, na.rm = TRUE),
                          mean_auc = mean(auc_values, na.rm = TRUE),
                          sd_auc = sd(auc_values, na.rm = TRUE))
xgb_summary


### combine all feature importances metrics
all_xgb_importances <- bind_rows(xgb_importances)

# data.frame of feature importances
# gain = total improvement in the model loss function (gives more importance to early splits in trees)
# cover = number of observations affected by that feature's splits (weighted average)
# frequency = how often the feature was used at split points
mean_importance <- all_xgb_importances %>%
  group_by(Feature = Feature) %>%
  summarise(mean_gain = mean(Gain, na.rm = TRUE),
            mean_cover = mean(Cover, na.rm = TRUE),
            mean_freq = mean(Frequency, na.rm = TRUE),
            freq_selected = n()) %>%
  arrange(desc(mean_gain))
head(mean_importance, 20)


#############################################################################
###   BASELINE XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING   ###
#############################################################################

# data to be used in the model
str(metabo_df)

# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation (250 models)
folds <- createMultiFolds(metabo_df$condition, k = 5, times = 50)

# create list to store performance metrics
xgb_metrics <- list() # list to store performance metrics
xgb_importances <- list() # list to store feature importances
best_nrounds_list <- list() # list to store best nrounds

# loop through the folds
for (key in names(folds)) {
  
  # splits the dataset into training and testing sets
  train_idx <- folds[[key]] # train indices
  train_data <- metabo_df[train_idx, ] # training data
  test_data <- metabo_df[-train_idx, ] # testing data
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), label = train_data$condition_numeric)
  dtest  <- xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), label = test_data$condition_numeric)
  
  # cross-validation with early stopping
  xgb_cv <- xgb.cv(data = dtrain,
                   nrounds = 5000,
                   nfold = 5,
                   early_stopping_rounds = 25,
                   maximize = FALSE, # lower logloss better
                   stratified = TRUE,
                   showsd = FALSE,
                   objective = "binary:logistic",
                   eval_metric = "logloss",
                   verbose = 0)
  
  # best number of boosting rounds
  best_nrounds <- xgb_cv$best_iteration
  
  # train XGBoost model with best number of boosting rounds
  xgb_model <- xgboost(data = dtrain,
                       objective = "binary:logistic",
                       eval_metric = "logloss",
                       nrounds = best_nrounds,
                       verbose = 0)
  
  # evaluate on test set
  preds_prob <- predict(xgb_model, dtest)
  preds_label <- ifelse(preds_prob > 0.5, "disease", "healthy")
  preds_label <- factor(preds_label, levels = c("healthy", "disease"))
  
  # calculate AUC
  auc_val <- auc(response = test_data$condition,
                 predictor = preds_prob,
                 levels = c("healthy", "disease"),
                 direction = "<")
  
  # generate confusion matrix
  cm <- confusionMatrix(preds_label, test_data$condition, positive = "disease")
  
  # store metrics
  xgb_metrics[[key]] <- list(cm = cm, auc = auc_val) # store confusion matrices
  imp <- xgb.importance(feature_names = all_feat_cols, model = xgb_model) # feature importances
  imp$Repeat_Fold <- key # indices for feature importances
  xgb_importances[[key]] <- imp # store feature importances
  best_nrounds_list[[key]] <- best_nrounds
}


### calculate performance statistics
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
auc_values <- numeric()

# loop through stored metrics to extract values
for (key in names(xgb_metrics)) {
  cm <- xgb_metrics[[key]]$cm
  auc_val <- xgb_metrics[[key]]$auc
  
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  auc_values <- c(auc_values, auc_val)
}

# data.frame of performance metrics
xgb_summary <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                          sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                          mean_f1 = mean(f1_score, na.rm = TRUE),
                          sd_f1 = sd(f1_score, na.rm = TRUE),
                          mean_sens = mean(sensitivity, na.rm = TRUE),
                          sd_sens = sd(sensitivity, na.rm = TRUE),
                          mean_spec = mean(specificity, na.rm = TRUE),
                          sd_spec = sd(specificity, na.rm = TRUE),
                          mean_auc = mean(auc_values, na.rm = TRUE),
                          sd_auc = sd(auc_values, na.rm = TRUE))
xgb_summary


### combine all feature importances metrics
all_xgb_importances <- bind_rows(xgb_importances)

# data.frame of feature importances
# gain = total improvement in the model loss function (gives more importance to early splits in trees)
# cover = number of observations affected by that feature's splits (weighted average)
# frequency = how often the feature was used at split points
mean_importance <- all_xgb_importances %>%
  group_by(Feature = Feature) %>%
  summarise(mean_gain = mean(Gain, na.rm = TRUE),
            mean_cover = mean(Cover, na.rm = TRUE),
            mean_freq = mean(Frequency, na.rm = TRUE),
            freq_selected = n()) %>%
  arrange(desc(mean_gain))
head(mean_importance, 20)


### early stopping metrics
# stability of logloss (early_stopping_rounds = 25)
head(xgb_cv$evaluation_log, 10)

# best nrounds from the model
summary(unlist(best_nrounds_list))


########################################################################
###   XGBOOST MODEL - HYPERPARAMETER TUNING + EARLY STOPPING - ETA   ###
########################################################################

### use tune_xgb_param function and associated analysis functions from ML_xgboost

# data to be used in the model (data.frame of all labels and features)
str(metabo_df)

# param name and grid of param values to test
param_grid_name <- "eta"
param_grid_values <- c(0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

# list of baseline hyperparameters to use (excluding param to test)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_eta <- tune_xgb_param(param_grid_name = "eta", 
                                  param_grid_values = param_grid_values,
                                  base_params = base_params,
                                  metagen = metabo_df,
                                  all_feat_cols = all_feat_cols,
                                  target_var = "condition",
                                  target_var_numeric = "condition_numeric",
                                  n_repeats = 50)


### functions to analyze the output of tune_xgb_param

### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
# uses output of tune_xgb_param
tune_summarize_performance(all_results_eta)


### feature importance (importance of features across parameter settings)
# uses output of tune_xgb_param, creates the all_importances_df
feature_freq_eta <- tune_feature_importance(all_results_eta, param_name = "eta")

# plot feature selection by frequency/selection across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_freq_eta, feature = "beta.alanine", param_name = "eta")

# plot feature frequency/selection by importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_gain", color_by = "eta")
tune_plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_cover", color_by = "eta")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_freq_eta, threshold_frac = 0.4, n_repeats = 50, param_name = "eta")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_freq_eta, param_name = "eta")


### logloss and overfitting analysis
# uses output of tune_xgb_param, creates logloss_df
logloss_eta <- tune_extract_logloss_df(all_results_eta, param_name = "eta")

# plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
# uses logloss_df 
tune_plot_logloss_curve(logloss_eta, type = "test", show_mean = TRUE, param_name = "eta")
tune_plot_logloss_curve(logloss_eta, type = "train", show_mean = FALSE, param_name = "eta")
tune_plot_logloss_curve(logloss_eta, type = "both", show_mean = FALSE, param_name = "eta")

# prepare logloss gap (calculates the generalization gap (test loss - train loss) over boosting rounds and parameter values)
# uses logloss_df, creates logloss_gap_df
logloss_gap_eta <- tune_logloss_gap(logloss_eta, param_name = "eta")

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_eta, param_name = "eta")


################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MAX DEPTH   ###
################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(2, 3, 4, 5, 6, 7)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_max_depth <- tune_xgb_param(param_grid_name = "max_depth", 
                                        param_grid_values = param_grid_values,
                                        base_params = base_params,
                                        metagen = metabo_df, 
                                        all_feat_cols = all_feat_cols,
                                        target_var = "condition",
                                        target_var_numeric = "condition_numeric",
                                        n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_max_depth)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_max_depth <- tune_feature_importance(all_results_max_depth, param_name = "max_depth")
feature_import_max_depth %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_max_depth, feature = "beta.alanine", param_name = "max_depth")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_max_depth, x = "freq_selected", y = "mean_gain", color_by = "max_depth")
tune_plot_feature_stability(feature_import_max_depth, x = "freq_selected", y = "mean_cover", color_by = "max_depth")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_max_depth, threshold_frac = 0.4, n_repeats = 10, param_name = "max_depth")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_max_depth, param_name = "max_depth")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_max_depth <- tune_extract_logloss_df(all_results_max_depth, param_name = "max_depth")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_max_depth, type = "test", show_mean = TRUE, param_name = "max_depth")
tune_plot_logloss_curve(logloss_max_depth, type = "train", show_mean = FALSE, param_name = "max_depth")
tune_plot_logloss_curve(logloss_max_depth, type = "both", show_mean = FALSE, param_name = "max_depth")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_max_depth <- tune_logloss_gap(logloss_max_depth, param_name = "max_depth")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_max_depth, param_name = "max_depth")


#######################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MIN CHILD WEIGHT   ###
#######################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(1, 2, 3, 4, 5, 6)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_min_child_weight <- tune_xgb_param(param_grid_name = "min_child_weight", 
                                               param_grid_values = param_grid_values,
                                               base_params = base_params,
                                               metagen = metabo_df, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_min_child_weight)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_min_child_weight <- tune_feature_importance(all_results_min_child_weight, param_name = "min_child_weight")
feature_import_min_child_weight %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_min_child_weight, feature = "beta.alanine", param_name = "min_child_weight")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_min_child_weight, x = "freq_selected", y = "mean_gain", color_by = "min_child_weight")
tune_plot_feature_stability(feature_import_min_child_weight, x = "freq_selected", y = "mean_cover", color_by = "min_child_weight")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_min_child_weight, threshold_frac = 0.4, n_repeats = 10, param_name = "min_child_weight")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_min_child_weight, param_name = "min_child_weight")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_min_child_weight <- tune_extract_logloss_df(all_results_min_child_weight, param_name = "min_child_weight")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_min_child_weight, type = "test", show_mean = FALSE, param_name = "min_child_weight")
tune_plot_logloss_curve(logloss_min_child_weight, type = "train", show_mean = FALSE, param_name = "min_child_weight")
tune_plot_logloss_curve(logloss_min_child_weight, type = "both", show_mean = FALSE, param_name = "min_child_weight")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_min_child_weight <- tune_logloss_gap(logloss_min_child_weight, param_name = "min_child_weight")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_min_child_weight, param_name = "min_child_weight")


################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - SUBSAMPLE   ###
################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_subsample <- tune_xgb_param(param_grid_name = "subsample", 
                                        param_grid_values = param_grid_values,
                                        base_params = base_params,
                                        metagen = metabo_df, 
                                        all_feat_cols = all_feat_cols,
                                        target_var = "condition",
                                        target_var_numeric = "condition_numeric",
                                        n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_subsample)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_subsample <- tune_feature_importance(all_results_subsample, param_name = "subsample")
feature_import_subsample %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_subsample, feature = "beta.alanine", param_name = "subsample")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_subsample, x = "freq_selected", y = "mean_gain", color_by = "subsample")
tune_plot_feature_stability(feature_import_subsample, x = "freq_selected", y = "mean_cover", color_by = "subsample")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_subsample, threshold_frac = 0.4, n_repeats = 10, param_name = "subsample")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_subsample, param_name = "subsample")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_subsample <- tune_extract_logloss_df(all_results_subsample, param_name = "subsample")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_subsample, type = "test", show_mean = TRUE, param_name = "subsample")
tune_plot_logloss_curve(logloss_subsample, type = "train", show_mean = FALSE, param_name = "subsample")
tune_plot_logloss_curve(logloss_subsample, type = "both", show_mean = FALSE, param_name = "subsample")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_subsample <- tune_logloss_gap(logloss_subsample, param_name = "subsample")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_subsample, param_name = "subsample")


#######################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - COLSAMPLE_BYTREE   ###
#######################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_colsample_bytree <- tune_xgb_param(param_grid_name = "colsample_bytree", 
                                               param_grid_values = param_grid_values,
                                               base_params = base_params,
                                               metagen = metabo_df, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_colsample_bytree)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_colsample_bytree <- tune_feature_importance(all_results_colsample_bytree, param_name = "colsample_bytree")
feature_import_colsample_bytree %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_colsample_bytree, feature = "beta.alanine", param_name = "colsample_bytree")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_colsample_bytree, x = "freq_selected", y = "mean_gain", color_by = "colsample_bytree")
tune_plot_feature_stability(feature_import_colsample_bytree, x = "freq_selected", y = "mean_cover", color_by = "colsample_bytree")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_colsample_bytree, threshold_frac = 0.4, n_repeats = 10, param_name = "colsample_bytree")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_colsample_bytree, param_name = "colsample_bytree")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_colsample_bytree <- tune_extract_logloss_df(all_results_colsample_bytree, param_name = "colsample_bytree")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_colsample_bytree, type = "test", show_mean = TRUE, param_name = "colsample_bytree")
tune_plot_logloss_curve(logloss_colsample_bytree, type = "train", show_mean = FALSE, param_name = "colsample_bytree")
tune_plot_logloss_curve(logloss_colsample_bytree, type = "both", show_mean = FALSE, param_name = "colsample_bytree")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_colsample_bytree <- tune_logloss_gap(logloss_colsample_bytree, param_name = "colsample_bytree")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_colsample_bytree, param_name = "colsample_bytree")


#######################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - COLSAMPLE_BYNODE   ###
#######################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_colsample_bynode <- tune_xgb_param(param_grid_name = "colsample_bynode", 
                                               param_grid_values = param_grid_values,
                                               base_params = base_params,
                                               metagen = metabo_df, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_colsample_bynode)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_colsample_bynode <- tune_feature_importance(all_results_colsample_bynode, param_name = "colsample_bynode")
feature_import_colsample_bynode %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_colsample_bynode, feature = "beta.alanine", param_name = "colsample_bynode")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_colsample_bynode, x = "freq_selected", y = "mean_gain", color_by = "colsample_bynode")
tune_plot_feature_stability(feature_import_colsample_bynode, x = "freq_selected", y = "mean_cover", color_by = "colsample_bynode")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_colsample_bynode, threshold_frac = 0.4, n_repeats = 10, param_name = "colsample_bynode")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_colsample_bynode, param_name = "colsample_bynode")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_colsample_bynode <- tune_extract_logloss_df(all_results_colsample_bynode, param_name = "colsample_bynode")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_colsample_bynode, type = "test", show_mean = TRUE, param_name = "colsample_bynode")
tune_plot_logloss_curve(logloss_colsample_bynode, type = "train", show_mean = FALSE, param_name = "colsample_bynode")
tune_plot_logloss_curve(logloss_colsample_bynode, type = "both", show_mean = FALSE, param_name = "colsample_bynode")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_colsample_bynode <- tune_logloss_gap(logloss_colsample_bynode, param_name = "colsample_bynode")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_colsample_bynode, param_name = "colsample_bynode")


############################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - GAMMA   ###
############################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 3)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    max_delta_step = 0)

all_results_gamma <- tune_xgb_param(param_grid_name = "gamma", 
                                    param_grid_values = param_grid_values,
                                    base_params = base_params,
                                    metagen = metabo_df, 
                                    all_feat_cols = all_feat_cols,
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric",
                                    n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_gamma)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_gamma <- tune_feature_importance(all_results_gamma, param_name = "gamma")
feature_import_gamma %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_gamma, feature = "beta.alanine", param_name = "gamma")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_gamma, x = "freq_selected", y = "mean_gain", color_by = "gamma")
tune_plot_feature_stability(feature_import_gamma, x = "freq_selected", y = "mean_cover", color_by = "gamma")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_gamma, threshold_frac = 0.4, n_repeats = 10, param_name = "gamma")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_gamma, param_name = "gamma")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_gamma <- tune_extract_logloss_df(all_results_gamma, param_name = "gamma")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_gamma, type = "test", show_mean = TRUE, param_name = "gamma")
tune_plot_logloss_curve(logloss_gamma, type = "train", show_mean = FALSE, param_name = "gamma")
tune_plot_logloss_curve(logloss_gamma, type = "both", show_mean = FALSE, param_name = "gamma")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_gamma <- tune_logloss_gap(logloss_gamma, param_name = "gamma")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_gamma, param_name = "gamma") 


#############################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - LAMBDA   ###
#############################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(0.1, 0.5, 1, 1.5, 3, 5)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_lambda <- tune_xgb_param(param_grid_name = "lambda", 
                                     param_grid_values = param_grid_values,
                                     base_params = base_params,
                                     metagen = metabo_df, 
                                     all_feat_cols = all_feat_cols,
                                     target_var = "condition",
                                     target_var_numeric = "condition_numeric",
                                     n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_lambda)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_lambda <- tune_feature_importance(all_results_lambda, param_name = "lambda")
feature_import_lambda %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_lambda, feature = "beta.alanine", param_name = "lambda")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_lambda, x = "freq_selected", y = "mean_gain", color_by = "lambda")
tune_plot_feature_stability(feature_import_lambda, x = "freq_selected", y = "mean_cover", color_by = "lambda")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_lambda, threshold_frac = 0.4, n_repeats = 10, param_name = "lambda")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_lambda, param_name = "lambda")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_lambda <- tune_extract_logloss_df(all_results_lambda, param_name = "lambda")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_lambda, type = "test", show_mean = TRUE, param_name = "lambda")
tune_plot_logloss_curve(logloss_lambda, type = "train", show_mean = FALSE, param_name = "lambda")
tune_plot_logloss_curve(logloss_lambda, type = "both", show_mean = FALSE, param_name = "lambda")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_lambda <- tune_logloss_gap(logloss_lambda, param_name = "lambda")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_lambda, param_name = "lambda") 


############################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - ALPHA   ###
############################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 3)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    gamma = 0,
                    max_delta_step = 0)

all_results_alpha <- tune_xgb_param(param_grid_name = "alpha", 
                                    param_grid_values = param_grid_values,
                                    base_params = base_params,
                                    metagen = metabo_df, 
                                    all_feat_cols = all_feat_cols,
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric",
                                    n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_alpha)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_alpha <- tune_feature_importance(all_results_alpha, param_name = "alpha")
feature_import_alpha %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_alpha, feature = "beta.alanine", param_name = "alpha")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_alpha, x = "freq_selected", y = "mean_gain", color_by = "alpha")
tune_plot_feature_stability(feature_import_alpha, x = "freq_selected", y = "mean_cover", color_by = "alpha")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_alpha, threshold_frac = 0.4, n_repeats = 10, param_name = "alpha")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_alpha, param_name = "alpha")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_alpha <- tune_extract_logloss_df(all_results_alpha, param_name = "alpha")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_alpha, type = "test", show_mean = TRUE, param_name = "alpha")
tune_plot_logloss_curve(logloss_alpha, type = "train", show_mean = FALSE, param_name = "alpha")
tune_plot_logloss_curve(logloss_alpha, type = "both", show_mean = FALSE, param_name = "alpha")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_alpha <- tune_logloss_gap(logloss_alpha, param_name = "alpha")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_alpha, param_name = "alpha") 


#####################################################################################
###   XGBOOST MODEL - HYPERPARAMETER TUNING + EARLY STOPPING - SCALE POS WEIGHT   ###
#####################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0,
                    max_delta_step = 0)

all_results_class_weight <- tune_xgb_param(param_grid_name = "scale_pos_weight", 
                                           param_grid_values = param_grid_values,
                                           base_params = base_params,
                                           metagen = metabo_df, 
                                           all_feat_cols = all_feat_cols,
                                           target_var = "condition",
                                           target_var_numeric = "condition_numeric",
                                           n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_class_weight)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_class_weight <- tune_feature_importance(all_results_class_weight, param_name = "scale_pos_weight")
feature_import_class_weight %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_class_weight, feature = "beta.alanine", param_name = "scale_pos_weight")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_class_weight, x = "freq_selected", y = "mean_gain", color_by = "scale_pos_weight")
tune_plot_feature_stability(feature_import_class_weight, x = "freq_selected", y = "mean_cover", color_by = "scale_pos_weight")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_class_weight, threshold_frac = 0.4, n_repeats = 10, param_name = "scale_pos_weight")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_class_weight, param_name = "scale_pos_weight")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_class_weight <- tune_extract_logloss_df(all_results_class_weight, param_name = "scale_pos_weight")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_class_weight, type = "test", show_mean = TRUE, param_name = "scale_pos_weight")
tune_plot_logloss_curve(logloss_class_weight, type = "train", show_mean = FALSE, param_name = "scale_pos_weight")
tune_plot_logloss_curve(logloss_class_weight, type = "both", show_mean = FALSE, param_name = "scale_pos_weight")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_class_weight <- tune_logloss_gap(logloss_class_weight, param_name = "scale_pos_weight")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_class_weight, param_name = "scale_pos_weight")


#####################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MAX_DELTA_STEP   ###
#####################################################################################

# data to be used in the model
str(metabo_df)

### run xgboost model
param_grid_values <- c(0, 1, 2, 3, 5, 10)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.005,
                    scale_pos_weight = 1,
                    max_depth = 6,
                    min_child_weight = 1,
                    subsample = 1,
                    colsample_bytree = 1,
                    colsample_bynode = 1,
                    lambda = 1,
                    alpha = 0,
                    gamma = 0)

all_results_max_delta_step <- tune_xgb_param(param_grid_name = "max_delta_step", 
                                             param_grid_values = param_grid_values,
                                             base_params = base_params,
                                             metagen = metabo_df, 
                                             all_feat_cols = all_feat_cols,
                                             target_var = "condition",
                                             target_var_numeric = "condition_numeric",
                                             n_repeats = 10)


### summarize performance metrics
# uses output of tune_xgb_param
tune_summarize_performance(all_results_max_delta_step)


### feature importance
# combine feature importances
# uses output of tune_xgb_param, creates the all_importances_df
feature_import_max_delta_step <- tune_feature_importance(all_results_max_delta_step, param_name = "max_delta_step")
feature_import_max_delta_step %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
tune_plot_feature_selection_frequency(feature_import_max_delta_step, feature = "beta.alanine", param_name = "max_delta_step")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
tune_plot_feature_stability(feature_import_max_delta_step, x = "freq_selected", y = "mean_gain", color_by = "max_delta_step")
tune_plot_feature_stability(feature_import_max_delta_step, x = "freq_selected", y = "mean_cover", color_by = "max_delta_step")

# number of features selected in more than threshold_frac models
# uses all_importances_df
tune_feature_stability_table(feature_import_max_delta_step, threshold_frac = 0.4, n_repeats = 10, param_name = "max_delta_step")

# mean frequency of feature selection per parameter value
# uses all_importances_df
tune_mean_feature_frequency(feature_import_max_delta_step, param_name = "max_delta_step")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_max_delta_step <- tune_extract_logloss_df(all_results_max_delta_step, param_name = "max_delta_step")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
tune_plot_logloss_curve(logloss_max_delta_step, type = "test", show_mean = TRUE, param_name = "max_delta_step")
tune_plot_logloss_curve(logloss_max_delta_step, type = "train", show_mean = FALSE, param_name = "max_delta_step")
tune_plot_logloss_curve(logloss_max_delta_step, type = "both", show_mean = FALSE, param_name = "max_delta_step")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_max_delta_step <- tune_logloss_gap(logloss_max_delta_step, param_name = "max_delta_step")

# plot logloss gap
# uses logloss_gap_df
tune_plot_logloss_gap(logloss_gap_max_delta_step, param_name = "max_delta_step") 


##################################################################################################
###   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF BAYES OPT   ###
##################################################################################################

# data to be used in the model
str(metabo_df)

scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step) {
  
  set.seed(1234)
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(metabo_df[, all_feat_cols]), 
                        label = metabo_df$condition_numeric)
  
  params <- list(objective = "binary:logistic",
                 eval_metric = "logloss",
                 eta = eta,
                 max_depth = as.integer(max_depth),
                 min_child_weight = as.integer(min_child_weight),
                 subsample = subsample,
                 colsample_bytree = colsample_bytree,
                 colsample_bynode = colsample_bynode,
                 gamma = gamma,
                 lambda = lambda,
                 alpha = alpha,
                 scale_pos_weight = scale_pos_weight,
                 max_delta_step = as.integer(max_delta_step))
  
  repeats <- 10
  repeat_logloss <- numeric(repeats)
  
  # loop over each repeat
  for (r in 1:repeats) {
    
    #  create 5-folds for cross-validation (stratified on condition)
    folds <- caret::createFolds(metabo_df$condition, k = 5, list = TRUE, returnTrain = FALSE)
    
    xgb_cv <- tryCatch({
      xgboost::xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 5000,
        folds = folds,
        stratified = TRUE,
        showsd = FALSE,
        early_stopping_rounds = 50,
        verbose = 0)
    }, error = function(e) return(NULL))
    
    if (is.null(xgb_cv)) {
      repeat_logloss[r] <- Inf
    } else {
      
      # take the best iteration of logloss for the repeat
      repeat_logloss[r] <- xgb_cv$evaluation_log$test_logloss_mean[xgb_cv$best_iteration]
    }
  }
  
  # average across repeats
  mean_logloss <- mean(repeat_logloss)
  
  # negative because Bayesian optimization maximizes the Score
  return(list(Score = -mean_logloss))
}

# used with metagenomics
bounds <- list(eta = c(0.0005, 0.01),
               max_depth = c(2L, 6L),
               min_child_weight = c(1L, 2L),
               subsample = c(0.9, 1.0),
               colsample_bytree = c(0.6, 1.0),
               colsample_bynode = c(0.6, 1.0),
               gamma = c(0, 4),
               lambda = c(0, 10),
               alpha = c(0, 4),
               scale_pos_weight = c(0.3, 2.0),
               max_delta_step = c(0L, 10L))


# resister back end
# parallelize evaluations of the scoring function
# each hyperparameter point being evaluated by Bayesian optimization can run on a separate core
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 44,
                   iters.n = 10,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ() 

# view results
stopstatus = optObj$stopStatus
score_sum = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
best_param_values = getBestPars(optObj)


##############################################################################
###   XGBOOST MODEL - EVALUATION OF MODEL WITH OPTIMIZED HYPERPARAMETERS   ###
##############################################################################

### use perf_xgb_evaluation function and associated analysis functions from ML_xgboost

# data to be used in the model (data.frame of all labels and features)
str(metabo_df)

# list of best hyperparameters to use (based on initial Bayesian optimization of hyperparameters)
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = best_param_values$eta,
                    scale_pos_weight = best_param_values$scale_pos_weight,
                    max_depth = best_param_values$max_depth,
                    min_child_weight = best_param_values$min_child_weight,
                    subsample = best_param_values$subsample,
                    colsample_bytree = best_param_values$colsample_bytree,
                    colsample_bynode = best_param_values$colsample_bynode,
                    lambda = best_param_values$lambda,
                    alpha = best_param_values$alpha,
                    gamma = best_param_values$gamma,
                    max_delta_step = best_param_values$max_delta_step)


### evaluate model using optimized hyperparameter values
perf_opt_params <- perf_xgb_evaluation(best_params = best_params,
                                       metagen = metabo_df,
                                       all_feat_cols = all_feat_cols,
                                       target_var = "condition",
                                       target_var_numeric = "condition_numeric",
                                       n_repeats = 50,
                                       n_folds = 5)


### functions to analyze the output of perf_xgb_evaluation
### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
# uses output of perf_xgb_evaluation
perf_summarize_performance(perf_opt_params)


### feature importance (importance of features across parameter settings)
# uses output of perf_xgb_evaluation, creates the all_importances_df
feature_freq_opt_params <- perf_feature_importance(perf_opt_params)
feature_freq_opt_params

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
perf_plot_feature_stability(feature_freq_opt_params, x = "freq_selected", y = "mean_gain")
perf_plot_feature_stability(feature_freq_opt_params, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
perf_feature_stability_table(feature_freq_opt_params, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
perf_mean_feature_frequency(feature_freq_opt_params)


### logloss and overfitting analysis
# uses output of perf_xgb_evaluation, creates logloss_df
logloss_opt_params <- perf_extract_logloss_df(perf_opt_params)
logloss_opt_params

# plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
# uses logloss_df
perf_plot_logloss_curve(logloss_opt_params, type = "test", show_mean = TRUE)
perf_plot_logloss_curve(logloss_opt_params, type = "train", show_mean = FALSE)
perf_plot_logloss_curve(logloss_opt_params, type = "both", show_mean = FALSE)

# prepare logloss gap (calculates the generalization gap (test loss - train loss) over boosting rounds)
# uses logloss_df, creates logloss_gap_df
logloss_gap_opt_params <- perf_prepare_logloss_gap(logloss_opt_params)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
perf_plot_logloss_gap(logloss_gap_opt_params)


#####################################################################
########   XGBOOST MODEL - FEATURE SELECTION USING TREE SHAP  #######
#####################################################################

### use xgb_shap_evaluation function and associated analysis functions from ML_xgboost

# data to be used in the model
str(metabo_df)

# list of best hyperparameters to use (from 50 repeats of 5-fold cross-validation parallelized on the cv folds)
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = best_param_values$eta,
                    scale_pos_weight = best_param_values$scale_pos_weight,
                    max_depth = best_param_values$max_depth,
                    min_child_weight = best_param_values$min_child_weight,
                    subsample = best_param_values$subsample,
                    colsample_bytree = best_param_values$colsample_bytree,
                    colsample_bynode = best_param_values$colsample_bynode,
                    lambda = best_param_values$lambda,
                    alpha = best_param_values$alpha,
                    gamma = best_param_values$gamma,
                    max_delta_step = best_param_values$max_delta_step)


### evaluate model and determine SHAP values using optimal hyperparameter values
shap_feature_selection <- xgb_shap_evaluation(best_params = best_params,
                                              metagen = metabo_df,
                                              all_feat_cols = all_feat_cols,
                                              target_var = "condition",
                                              target_var_numeric = "condition_numeric",
                                              n_repeats = 50,
                                              n_folds = 5)


### functions to analyze the output of xgb_shap_evaluation
### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
# uses output of xgb_shap_evaluation
perf_summarize_performance(shap_feature_selection)


### feature importance (importance of features across parameter settings)
# uses output of xgb_shap_evaluation, creates the all_importances_df
feature_freq_opt_params <- perf_feature_importance(shap_feature_selection) %>%
  arrange(desc(freq_selected))
feature_freq_opt_params

# plot feature frequency/selection by importance (mean gain or mean cover)
# uses the all_importances_df
perf_plot_feature_stability(feature_freq_opt_params, x = "freq_selected", y = "mean_gain")
perf_plot_feature_stability(feature_freq_opt_params, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
perf_feature_stability_table(feature_freq_opt_params, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
perf_mean_feature_frequency(feature_freq_opt_params)


### logloss and overfitting analysis
# uses output of xgb_shap_evaluation, creates logloss_df
logloss_opt_params <- perf_extract_logloss_df(shap_feature_selection)
logloss_opt_params

# plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
# uses logloss_df
perf_plot_logloss_curve(logloss_opt_params, type = "test", show_mean = TRUE)
perf_plot_logloss_curve(logloss_opt_params, type = "train", show_mean = FALSE)
perf_plot_logloss_curve(logloss_opt_params, type = "both", show_mean = FALSE)

# prepare logloss gap (calculates the generalization gap (test loss - train loss) over boosting rounds)
# uses logloss_df, creates logloss_gap_df
logloss_gap_opt_params <- perf_prepare_logloss_gap(logloss_opt_params)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
perf_plot_logloss_gap(logloss_gap_opt_params)


### analyze SHAP feature values and frequency of selection (similar to Boruta)
# aggregate median importance of SHAP values
shap_df <- do.call(rbind, shap_feature_selection$final_evaluation$shap_metrics)

# summarize across repeats
summary_shap_df <- shap_df %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_meanSHAP = mean(mean_abs_shap, na.rm = TRUE),
                   sd_meanSHAP = sd(mean_abs_shap, na.rm = TRUE),
                   count = n()) %>%
  arrange(desc(mean_meanSHAP)) 


### plot frequency of selection and SHAP values of features
n_repeats <- 50 # number of repeats used in function

# plot selection frequency of features
ggplot(summary_shap_df, aes(x = reorder(feature, count), y = count)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Feature selection frequency across repeats",
       x = "Feature", y = paste("Number of runs selected (out of", n_repeats, ")"))  


# plot average mean SHAP values for features
summary_shap_df$frequency <- ifelse(summary_shap_df$count >= n_repeats*0.80, "in_atleast_80%", "other") # add column: in all repeats or not

# plot features selected in half of repeats and highlight features selected in all repeats
ggplot(summary_shap_df[summary_shap_df$count >= n_repeats*0.60, ], 
       aes(x = reorder(feature, mean_meanSHAP), y = mean_meanSHAP, fill = frequency)) +
  scale_fill_manual(values = c("in_atleast_80%" = "indianred3", "other" = "steelblue")) +
  geom_col() + coord_flip() + theme_minimal(base_size = 12) +
  labs(title = "Average mean SHAP value of features",
       x = "Feature", y = "Average mean SHAP value", fill = "Feature selection frequency")


### features to keep
selected_features_df <- summary_shap_df %>% 
  filter(count >= n_repeats*0.60) %>% # frequency of feature selection
  head(10) # top mean_meanSHAP 

selected_features <- selected_features_df$feature # features to keep

# subset metabo_df to just confirmed features
shap_metabo_df <- metabo_df[, c("condition", "condition_numeric", selected_features)]
str(shap_metabo_df)


##############################################################################################################
########   XGBOOST MODEL - EVALUATION OF MODEL WITH BEST HYPERPARAMETERS USING SHAP-SELECTED FEATURES  #######
##############################################################################################################

### use perf_xgb_evaluation function and associated analysis functions from ML_xgboost

# data to be used in the model (data.frame of all labels and features)
str(shap_metabo_df)

# set predictor columns
subset_feat_cols <- setdiff(colnames(shap_metabo_df), c("condition", "condition_numeric"))

# list of best hyperparameters to use (based on initial Bayesian optimization of hyperparameters)
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = best_param_values$eta,
                    scale_pos_weight = best_param_values$scale_pos_weight,
                    max_depth = best_param_values$max_depth,
                    min_child_weight = best_param_values$min_child_weight,
                    subsample = best_param_values$subsample,
                    colsample_bytree = best_param_values$colsample_bytree,
                    colsample_bynode = best_param_values$colsample_bynode,
                    lambda = best_param_values$lambda,
                    alpha = best_param_values$alpha,
                    gamma = best_param_values$gamma,
                    max_delta_step = best_param_values$max_delta_step)


### evaluate model using optimized hyperparameter values
subset_performance <- perf_xgb_evaluation(best_params = best_params,
                                          metagen = shap_metabo_df, # subset based on selected features
                                          all_feat_cols = subset_feat_cols,
                                          target_var = "condition",
                                          target_var_numeric = "condition_numeric",
                                          n_repeats = 50,
                                          n_folds = 5)


### functions to analyze the output of perf_xgb_evaluation
### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
# uses output of perf_xgb_evaluation
perf_summarize_performance(subset_performance)


### feature importance (importance of features across parameter settings)
# uses output of perf_xgb_evaluation, creates the all_importances_df
feature_freq_subset <- perf_feature_importance(subset_performance)
feature_freq_subset

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
perf_plot_feature_stability(feature_freq_subset, x = "freq_selected", y = "mean_gain")
perf_plot_feature_stability(feature_freq_subset, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
perf_feature_stability_table(feature_freq_subset, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
perf_mean_feature_frequency(feature_freq_subset)


### logloss and overfitting analysis
# uses output of perf_xgb_evaluation, creates logloss_df
logloss_subset <- perf_extract_logloss_df(subset_performance)
logloss_subset

# plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
# uses logloss_df
perf_plot_logloss_curve(logloss_subset, type = "test", show_mean = TRUE)
perf_plot_logloss_curve(logloss_subset, type = "train", show_mean = FALSE)
perf_plot_logloss_curve(logloss_subset, type = "both", show_mean = FALSE)

# prepare logloss gap (calculates the generalization gap (test loss - train loss) over boosting rounds)
# uses logloss_df, creates logloss_gap_df
logloss_gap_subset <- perf_prepare_logloss_gap(logloss_subset)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
perf_plot_logloss_gap(logloss_gap_subset)


#########################################################################################################
########   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS USING SHAP-SELECTED FEATURES  #######
#########################################################################################################

# data to be used in the model (data.frame of all labels and features)
str(shap_metabo_df)

# set predictor columns
subset_feat_cols <- setdiff(colnames(shap_metabo_df), c("condition", "condition_numeric"))


scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step) {
  
  set.seed(1234)
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(shap_metabo_df[, subset_feat_cols]), 
                        label = shap_metabo_df$condition_numeric)
  
  params <- list(objective = "binary:logistic",
                 eval_metric = "logloss",
                 eta = eta,
                 max_depth = as.integer(max_depth),
                 min_child_weight = as.integer(min_child_weight),
                 subsample = subsample,
                 colsample_bytree = colsample_bytree,
                 colsample_bynode = colsample_bynode,
                 gamma = gamma,
                 lambda = lambda,
                 alpha = alpha,
                 scale_pos_weight = scale_pos_weight,
                 max_delta_step = as.integer(max_delta_step))
  
  repeats <- 10
  repeat_logloss <- numeric(repeats)
  
  # loop over each repeat
  for (r in 1:repeats) {
    
    #  create 5-folds for cross-validation (stratified on condition)
    folds <- caret::createFolds(shap_metabo_df$condition, k = 5, list = TRUE, returnTrain = FALSE)
    
    xgb_cv <- tryCatch({
      xgboost::xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 5000,
        folds = folds,
        stratified = TRUE,
        showsd = FALSE,
        early_stopping_rounds = 50,
        verbose = 0)
    }, error = function(e) return(NULL))
    
    if (is.null(xgb_cv)) {
      repeat_logloss[r] <- Inf
    } else {
      
      # take the best iteration of logloss for the repeat
      repeat_logloss[r] <- xgb_cv$evaluation_log$test_logloss_mean[xgb_cv$best_iteration]
    }
  }
  
  # average across repeats
  mean_logloss <- mean(repeat_logloss)
  
  # negative because Bayesian optimization maximizes the Score
  return(list(Score = -mean_logloss))
}

# define parameter bounds
bounds <- list(eta = c(0.0005, 0.01),
               max_depth = c(2L, 6L),
               min_child_weight = c(1L, 2L),
               subsample = c(0.9, 1.0),
               colsample_bytree = c(0.6, 1.0),
               colsample_bynode = c(0.6, 1.0),
               gamma = c(0, 4),
               lambda = c(0, 10),
               alpha = c(0, 4),
               scale_pos_weight = c(0.3, 2.0),
               max_delta_step = c(0L, 10L))

# resister back end
# parallelize evaluations of the scoring function
# each hyperparameter point being evaluated by Bayesian optimization can run on a separate core
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 44,
                   iters.n = 10,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ() 

# view results
stopstatus = optObj$stopStatus
score_sum = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
best_param_values_shap = getBestPars(optObj)


###################################################################################################################
########   XGBOOST MODEL - EVALUATION OF MODEL WITH BEST SHAP HYPERPARAMETERS USING SHAP-SELECTED FEATURES  #######
###################################################################################################################

### use perf_xgb_evaluation function and associated analysis functions from ML_xgboost

# data to be used in the model (data.frame of all labels and features)
str(shap_metabo_df)

# set predictor columns
subset_feat_cols <- setdiff(colnames(shap_metabo_df), c("condition", "condition_numeric"))

# list of best hyperparameters to use (based on initial Bayesian optimization of hyperparameters)
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = best_param_values_shap$eta,
                    scale_pos_weight = best_param_values_shap$scale_pos_weight,
                    max_depth = best_param_values_shap$max_depth,
                    min_child_weight = best_param_values_shap$min_child_weight,
                    subsample = best_param_values_shap$subsample,
                    colsample_bytree = best_param_values_shap$colsample_bytree,
                    colsample_bynode = best_param_values_shap$colsample_bynode,
                    lambda = best_param_values_shap$lambda,
                    alpha = best_param_values_shap$alpha,
                    gamma = best_param_values_shap$gamma,
                    max_delta_step = best_param_values_shap$max_delta_step)


### evaluate model using optimized hyperparameter values
subset_params_performance <- perf_xgb_evaluation(best_params = best_params,
                                                 metagen = shap_metabo_df, # subset based on selected features
                                                 all_feat_cols = subset_feat_cols,
                                                 target_var = "condition",
                                                 target_var_numeric = "condition_numeric",
                                                 n_repeats = 50,
                                                 n_folds = 5)


### functions to analyze the output of perf_xgb_evaluation
### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
# uses output of perf_xgb_evaluation
perf_summarize_performance(subset_params_performance)


### feature importance (importance of features across parameter settings)
# uses output of perf_xgb_evaluation, creates the all_importances_df
feature_freq_subset <- perf_feature_importance(subset_params_performance)
feature_freq_subset

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
perf_plot_feature_stability(feature_freq_subset, x = "freq_selected", y = "mean_gain")
perf_plot_feature_stability(feature_freq_subset, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
perf_feature_stability_table(feature_freq_subset, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
perf_mean_feature_frequency(feature_freq_subset)


### logloss and overfitting analysis
# uses output of perf_xgb_evaluation, creates logloss_df
logloss_subset <- perf_extract_logloss_df(subset_params_performance)
logloss_subset

# plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
# uses logloss_df
perf_plot_logloss_curve(logloss_subset, type = "test", show_mean = TRUE)
perf_plot_logloss_curve(logloss_subset, type = "train", show_mean = FALSE)
perf_plot_logloss_curve(logloss_subset, type = "both", show_mean = FALSE)

# prepare logloss gap (calculates the generalization gap (test loss - train loss) over boosting rounds)
# uses logloss_df, creates logloss_gap_df
logloss_gap_subset <- perf_prepare_logloss_gap(logloss_subset)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
perf_plot_logloss_gap(logloss_gap_subset)


###################################################################################################
###   XGBOOST LOGLOSS MODEL - TRAIN FINAL MODEL WITH SHAP SUBSET DATA AND BEST HYPERPARAMETERS  ###
###################################################################################################

# data to be used in the model (data.frame of all labels and features)
str(shap_metabo_df)

# set predictor columns
subset_feat_cols <- setdiff(colnames(shap_metabo_df), c("condition", "condition_numeric"))

# optimum number of nrounds chosen during final evaluation
best_nrounds <- round(mean(unlist(subset_params_performance$final_evaluation$best_nrounds)))

# best hyperparameter values determined by Bayesian optimization
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = best_param_values_shap$eta,
                    scale_pos_weight = best_param_values_shap$scale_pos_weight,
                    max_depth = best_param_values_shap$max_depth,
                    min_child_weight = best_param_values_shap$min_child_weight,
                    subsample = best_param_values_shap$subsample,
                    colsample_bytree = best_param_values_shap$colsample_bytree,
                    colsample_bynode = best_param_values_shap$colsample_bynode,
                    lambda = best_param_values_shap$lambda,
                    alpha = best_param_values_shap$alpha,
                    gamma = best_param_values_shap$gamma,
                    max_delta_step = best_param_values_shap$max_delta_step)

# train the model on the full dataset
dtrain_full <- xgboost::xgb.DMatrix(data = as.matrix(shap_metabo_df[, subset_feat_cols]),
                                    label = shap_metabo_df$condition_numeric)

final_model <- xgboost::xgb.train(params = best_params,
                                  data = dtrain_full,
                                  nrounds = best_nrounds,
                                  verbose = 1)


#################################################################################
###   XGBOOST LOGLOSS MODEL - SHAP VALUES - DEPENDENCE AND INTERACTION PLOTS  ###
#################################################################################

selected_features

# compute tree SHAP values
shap_values <- predict(final_model, newdata = dtrain_full, predcontrib = TRUE)
shap_df <- as.data.frame(shap_values)
shap_df$BIAS <- NULL  # remove bias term


### SHAP summary plot
shap_long <- shap.prep(xgb_model = final_model, X_train = as.matrix(shap_metabo_df[, subset_feat_cols]))
shap.plot.summary(shap_long)


### mean absolute SHAP value per feature
shap_mean_abs <- sort(colMeans(abs(shap_df)), decreasing = TRUE)
shap_mean_abs <- as.data.frame(shap_mean_abs) # covert to data.frame
shap_mean_abs$feature <- rownames(shap_mean_abs)
shap_mean_abs

# plot mean absolute SHAP value per feature
ggplot(shap_mean_abs, aes(x = reorder(feature, shap_mean_abs), y = shap_mean_abs)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Mean absolute SHAP value", 
       x = "Feature", y = "Mean absolute SHAP value")


### SHAP dependence plots (how SHAP values for a given feature vary as the input values for the feature vary)
# plot dependence plot for feature of interest (shap.plot.dependence and ggplot)
shap.plot.dependence(data_long = shap_long, x = "beta.alanine", y = NULL)
shap.plot.dependence(data_long = shap_long, x = "myristate..14.0.", y = NULL)

# wide table of CLR-transformed relative abundance
rfvalue_wide <- shap_long %>%
  select(ID, variable, rfvalue) %>%
  pivot_wider(names_from = variable, values_from = rfvalue)


feature_name <- "beta.alanine"
shap_dep <- shap_long %>% filter(variable == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "X7.ketodeoxycholate"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- CLR abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- CLR abun"))


feature_name <- "myristate..14.0."
shap_dep <- shap_long %>% filter(variable == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "X2.dimethylaminoethanol"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- CLR abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- CLR abun"))


### SHAP interaction values (how pairs of features interact in affecting the prediction)
interaction_values <- predict(final_model,
                              newdata = as.matrix(shap_metabo_df[, subset_feat_cols]),
                              predinteraction = TRUE)
interaction_values <- interaction_values[, -ncol(interaction_values), -ncol(interaction_values)] # remove BIAS term
mean_interactions <- apply(abs(interaction_values), c(2, 3), mean) # average absolute interaction strengths

# set row and column names
feature_names <- colnames(shap_metabo_df[, subset_feat_cols])
rownames(mean_interactions) <- feature_names
colnames(mean_interactions) <- feature_names

# convert to long format for plotting
interaction_long <- as.data.frame(mean_interactions) %>%
  rownames_to_column("Feature1") %>%
  pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "InteractionStrength")

# remove self-interactions (diagonal) - to see other interactions more clearly
interaction_long <- interaction_long %>%
  filter(Feature1 != Feature2)

# plot interaction heatmap
ggplot(interaction_long, aes(x = Feature1, y = Feature2, fill = InteractionStrength)) +
  geom_tile() + scale_fill_viridis_c() + theme_minimal() + coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SHAP feature interactions", fill = "Mean absolute interaction")


### feature-specific interactions
target_feat <- "beta.alanine" # feature of interest
target_feat <- "X.24683" # feature of interest

# data.frame with average feature specific interaction for target_feat 
target_interact <- data.frame(feature = feature_names,
                              interaction = mean_interactions[which(feature_names == target_feat), ])
target_interact <- target_interact[target_interact$feature != target_feat, ] # remove self-interaction

# plot feature-specific interactions
ggplot(target_interact, aes(x = reorder(feature, interaction), y = interaction)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = paste("SHAP interactions with", target_feat),
       x = "Interacting feature", y = "Mean absolute interaction")


### mean absolute SHAP value per class per feature
shap_df$condition <- shap_metabo_df$condition
abs_shap_by_class <- shap_df %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

# plot mean absolute SHAP value per class per feature
ggplot(abs_shap_by_class, aes(x = reorder(name, mean_abs_shap), y = mean_abs_shap, fill = condition)) +
  geom_col(position = "dodge") +
  coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean absolute SHAP values",
       x = "Feature", y = "Mean abs SHAP", fill = "Condition")


### mean SHAP value per class per feature
shap_df$condition <- shap_metabo_df$condition
shap_by_class <- shap_df %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_shap = mean(value), .groups = "drop")

# plot mean SHAP value per class per feature
ggplot(shap_by_class, aes(x = reorder(name, mean_shap), y = mean_shap, fill = condition)) +
  geom_col(position = "dodge") +
  coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean SHAP values",
       x = "Feature", y = "Mean SHAP", fill = "Condition")


### how the impact of a driven feature on model prediction varies between healthy and disease samples
# histogram of SHAP values
ggplot(shap_df, aes(x = beta.alanine, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - beta.alanine",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_df, aes(x = condition, y = beta.alanine, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for beta.alanine by condition",
       y = "SHAP value", x = "Condition")

# box plot of CLR-abundance
ggplot(shap_metabo_df, aes(x = condition, y = beta.alanine, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "CLR abundance of beta.alanine by condition",
       y = "CLR Abundance", x = "Condition")


# histogram of SHAP values
ggplot(shap_df, aes(x = X.24683, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - X.24683",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_df, aes(x = condition, y = X.24683, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for X.24683 by condition",
       y = "SHAP value", x = "Condition")

# box plot of CLR-abundance
ggplot(shap_metabo_df, aes(x = condition, y = X.24683, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "CLR abundance of X.24683 by condition",
       y = "CLR abundance", x = "Condition")


### SHAP values versus predicted probabilities
shap_df$ID <- rownames(shap_df)
shap_df$pred_prob <- predict(final_model, newdata = dtrain_full)

# X.24683
plot_df <- shap_df %>%
  select(ID, X.24683, condition, pred_prob)

ggplot(plot_df, aes(x = X.24683, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "X.24683 SHAP value versus prediction probability",
       x = "SHAP value for X.24683", y = "Predicted probability (disease)")

# beta.alanine
plot_df <- shap_df %>%
  select(ID, beta.alanine, condition, pred_prob)

ggplot(plot_df, aes(x = beta.alanine, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "beta.alanine SHAP value versus prediction probability",
       x = "SHAP value for beta.alanine", y = "Predicted probability (disease)")


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
#   [1] DiceKriging_1.6.0             SHAPforxgboost_0.1.3          doParallel_1.0.17            
# [4] iterators_1.0.14              foreach_1.5.2                 ParBayesianOptimization_1.2.6
# [7] Matrix_1.7-3                  MLmetrics_1.1.3               pROC_1.19.0.1                
# [10] caret_7.0-1                   lattice_0.22-7                xgboost_1.7.11.1             
# [13] impute_1.82.0                 readxl_1.4.5                  lubridate_1.9.4              
# [16] forcats_1.0.0                 stringr_1.5.1                 dplyr_1.1.4                  
# [19] purrr_1.1.0                   readr_2.1.5                   tidyr_1.3.1                  
# [22] tibble_3.3.0                  tidyverse_2.0.0               ggplot2_3.5.2                
# 
# loaded via a namespace (and not attached):
#   [1] rlang_1.1.6          magrittr_2.0.3       e1071_1.7-16         compiler_4.5.0      
# [5] vctrs_0.6.5          reshape2_1.4.4       lhs_1.2.0            pkgconfig_2.0.3     
# [9] crayon_1.5.3         backports_1.5.0      labeling_0.4.3       utf8_1.2.6          
# [13] prodlim_2025.04.28   tzdb_0.5.0           jsonlite_2.0.0       recipes_1.3.1       
# [17] tweenr_2.0.3         broom_1.0.9          R6_2.6.1             stringi_1.8.7       
# [21] RColorBrewer_1.1-3   parallelly_1.45.1    car_3.1-3            rpart_4.1.24        
# [25] cellranger_1.1.0     Rcpp_1.1.0           future.apply_1.20.0  splines_4.5.0       
# [29] nnet_7.3-20          timechange_0.3.0     tidyselect_1.2.1     rstudioapi_0.17.1   
# [33] dichromat_2.0-0.1    abind_1.4-8          timeDate_4041.110    codetools_0.2-20    
# [37] listenv_0.9.1        plyr_1.8.9           withr_3.0.2          future_1.67.0       
# [41] survival_3.8-3       proxy_0.4-27         polyclip_1.10-7      pillar_1.11.0       
# [45] ggpubr_0.6.1         carData_3.0-5        checkmate_2.3.2      stats4_4.5.0        
# [49] generics_0.1.4       dbscan_1.2.2         hms_1.1.3            scales_1.4.0        
# [53] globals_0.18.0       class_7.3-23         glue_1.8.0           tools_4.5.0         
# [57] data.table_1.17.8    ModelMetrics_1.2.2.2 gower_1.0.2          ggsignif_0.6.4      
# [61] grid_4.5.0           ipred_0.9-15         nlme_3.1-168         BBmisc_1.13         
# [65] ggforce_0.5.0        Formula_1.2-5        cli_3.6.5            lava_1.8.1          
# [69] gtable_0.3.6         rstatix_0.7.2        digest_0.6.37        farver_2.1.2        
# [73] lifecycle_1.0.4      hardhat_1.4.1        MASS_7.3-65

