# XGBoost Workflow - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(tidyverse)
library(compositions)
library(xgboost)
library(caret)
library(pROC)
library(MLmetrics)
library(Matrix)
library(ParBayesianOptimization)
library(doParallel)
library(foreach)


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

### merge metadata and feature table
metagen <- merge(meta, feat, by = "sample_name", all.x = TRUE)
metagen <- metagen[,-1] # remove sample_name

# make sure names are syntactically valid 
colnames(metagen) <- make.names(colnames(metagen))


### for xgboost analysis
# convert condition to numeric (0 = healthy, 1 = disease)
metagen$condition_numeric <- ifelse(metagen$condition == "disease", 1, 0)
# set factor labels (healthy = negative, disease = positive)
metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease"))

# set predictor columns
all_feat_cols <- setdiff(colnames(metagen), c("condition", "condition_numeric"))


#####################################################################
########   BASELINE XGBOOST MODEL - DEFAULT HYPERPARAMETERS   #######
#####################################################################


# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation (250 models)
folds <- createMultiFolds(metagen$condition, k = 5, times = 50)

# create list to store performance metrics
xgb_metrics <- list() # list to store performance metrics
xgb_importances <- list() # list to store feature importances

# loop through the folds
for (key in names(folds)) {

  # splits the dataset into training and testing sets
  train_idx <- folds[[key]] # train indices
  train_data <- metagen[train_idx, ] # training data
  test_data <- metagen[-train_idx, ] # testing data
  
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


######################################################################################
########   BASELINE XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING   #######
######################################################################################


# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation (250 models)
folds <- createMultiFolds(metagen$condition, k = 5, times = 50)

# create list to store performance metrics
xgb_metrics <- list() # list to store performance metrics
xgb_importances <- list() # list to store feature importances
best_nrounds_list <- list() # list to store best nrounds

# loop through the folds
for (key in names(folds)) {
  
  # splits the dataset into training and testing sets
  train_idx <- folds[[key]] # train indices
  train_data <- metagen[train_idx, ] # training data
  test_data <- metagen[-train_idx, ] # testing data
  
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


### best_nrounds is very low -> model overfitting very quickly
### default hyperparameters too aggressive for microbiome data (high-dimensional, low signal-to-noise ratio, small sample sizes)


#############################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - LEARNING RATE   #######
#############################################################################################


# list of baseline hyperparameters to use (excluding eta)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    max_depth = 6,
                    subsample = 1,
                    colsample_bytree = 1,
                    lambda = 1,
                    alpha = 0,
                    min_child_weight = 1)

# grid of eta values (learning rate)
eta_grid <- c(0.005, 0.01, 0.05, 0.1, 0.3)

# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation
folds <- createMultiFolds(metagen$condition, k = 5, times = 50)

# create list to store results (metrics, feature importances, best nrounds) for each eta setting
all_eta_results <- list()

# loop over eta values
for (e in eta_grid) {
  cat("Running with eta =", e, "\n")
  
  # update params with current eta to test
  params <- base_params
  params$eta <- e
  
  # create lists to store results for this setting
  xgb_metrics <- list()
  xgb_importances <- list()
  best_nrounds_list <- list()
  eval_logs <- list()
  
  # loop over each of the 250 cross-validation folds
  for (key in names(folds)) {
    
    # splits the dataset into training and testing sets for the current fold
    train_idx <- folds[[key]] # train indices 
    train_data <- metagen[train_idx, ] # training data
    test_data  <- metagen[-train_idx, ] # testing data
    
    # convert to DMatrix
    dtrain <- xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), label = train_data$condition_numeric)
    dtest  <- xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), label = test_data$condition_numeric)
    
    # run internal 5-fold cross-validation with early stopping
    xgb_cv <- xgb.cv(data = dtrain,
                     params = params,
                     nrounds = 5000,
                     nfold = 5,
                     early_stopping_rounds = 50,
                     maximize = FALSE,
                     stratified = TRUE,
                     showsd = FALSE,
                     verbose = 0)
    
    # best number of boosting rounds found in the internal 5-fold cross-validation above
    best_nrounds <- xgb_cv$best_iteration
    
    # store evaluation_log from xgb.cv()
    eval_log <- xgb_cv$evaluation_log
    eval_log$fold_id <- key 
    eval_logs[[key]] <- eval_log
    
    # train final xgboost model with best nrounds
    xgb_model <- xgboost(data = dtrain,
                         params = params,
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
    
    # calculate test logloss
    eps <- 1e-15  # avoid log(0)
    prob_clipped <- pmin(pmax(preds_prob, eps), 1 - eps)
    logloss <- -mean(test_data$condition_numeric * log(prob_clipped) +
                       (1 - test_data$condition_numeric) * log(1 - prob_clipped))
    
    # store metrics, importances and best nrounds
    xgb_metrics[[key]] <- list(cm = cm, auc = auc_val, logloss = logloss)
    imp <- xgb.importance(feature_names = all_feat_cols, model = xgb_model)
    imp$Repeat_Fold <- key
    xgb_importances[[key]] <- imp
    best_nrounds_list[[key]] <- best_nrounds
  }
  
  ### summarize performance
  # create vectors to store metrics
  balanced_accuracy <- numeric()
  f1_score <- numeric()
  precision <- numeric()
  sensitivity <- numeric()
  specificity <- numeric()
  auc_values <- numeric()
  logloss_values <- numeric() 
  nrounds_values <- unlist(best_nrounds_list)
  
  # loop through stored lists to extract metrics
  for (key in names(xgb_metrics)) {
    cm <- xgb_metrics[[key]]$cm
    auc_val <- xgb_metrics[[key]]$auc
    logloss_val <- xgb_metrics[[key]]$logloss
    
    balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
    f1_score <- c(f1_score, cm$byClass["F1"])
    precision <- c(precision, cm$byClass["Precision"])
    sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
    specificity <- c(specificity, cm$byClass["Specificity"])
    auc_values <- c(auc_values, auc_val)
    logloss_values <- c(logloss_values, logloss_val)
  }
  
  # aggregate stored metrics into a data.frame
  xgb_summary <- data.frame(eta = e,
                            mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                            sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                            mean_f1 = mean(f1_score, na.rm = TRUE),
                            sd_f1 = sd(f1_score, na.rm = TRUE),
                            mean_precision = mean(precision, na.rm = TRUE),
                            sd_precision = sd(precision, na.rm = TRUE),
                            mean_sens = mean(sensitivity, na.rm = TRUE),
                            sd_sens = sd(sensitivity, na.rm = TRUE),
                            mean_spec = mean(specificity, na.rm = TRUE),
                            sd_spec = sd(specificity, na.rm = TRUE),
                            mean_auc = mean(auc_values, na.rm = TRUE),
                            sd_auc = sd(auc_values, na.rm = TRUE),
                            mean_logloss = mean(logloss_values, na.rm = TRUE),
                            sd_logloss = sd(logloss_values, na.rm = TRUE),
                            mean_nrounds = mean(nrounds_values, na.rm = TRUE),
                            sd_nrounds = sd(nrounds_values, na.rm = TRUE))
  
  # aggregate feature importances
  all_importances <- bind_rows(xgb_importances)
  mean_importance <- all_importances %>%
    group_by(Feature = Feature) %>%
    summarise(mean_gain = mean(Gain, na.rm = TRUE),
              mean_cover = mean(Cover, na.rm = TRUE),
              mean_freq = mean(Frequency, na.rm = TRUE),
              freq_selected = n()) %>%
    arrange(desc(mean_gain))
  
  # store all values under current eta value
  all_eta_results[[paste0("eta_", e)]] <- list(summary = xgb_summary,
                                               feature_importance = mean_importance,
                                               raw_metrics = xgb_metrics,
                                               best_nrounds = best_nrounds_list,
                                               eval_logs = eval_logs)
  
}


### combine all summaries (metrics and SD for all metrics per eta) for comparison
performance_summary <- bind_rows(lapply(all_eta_results, `[[`, "summary"))
performance_summary


### feature importances
# feature importances by eta
head(all_eta_results[["eta_0.3"]]$feature_importance, 20)

# combine feature importances from all eta values into one tibble
all_importances <- lapply(names(all_eta_results), function(name) {
  df <- all_eta_results[[name]]$feature_importance
  df$eta <- sub("eta_", "", name)  # add eta as a column
  df
})
all_importances_df <- bind_rows(all_importances)
all_importances_df <- all_importances_df %>%
  arrange(desc(freq_selected)) # arrange by freq_selected
print(all_importances_df, n = 20)

all_importances_df <- all_importances_df %>%
  mutate(eta_clean = as.numeric(gsub("eta_", "", eta))) # convert eta string to numeric

# plot frequency of selection by eta value
top_feat <- "Acutalibacter_muris"
ggplot(filter(all_importances_df, Feature == top_feat), aes(x = eta_clean, y = freq_selected)) +
  geom_line() + geom_point() +
  ggtitle(paste("Frequency of selection vs eta for", top_feat)) + xlab("Eta") + ylab("Frequency of selection")

# number of features selected more than 75/250 models
feature_stability <- all_importances_df %>%
  group_by(eta) %>%
  summarise(n_features_selected_75plus = sum(freq_selected >= 75))
feature_stability

# mean of frequency of feature selection
feature_mean_freq <- all_importances_df %>%
  group_by(eta) %>%
  summarise(mean_frequency_selection = mean(freq_selected))
feature_mean_freq

# feature stability by mean gain
ggplot(all_importances_df, aes(x = freq_selected, y = mean_gain, color = as.factor(eta))) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Feature importance stability across eta values", x = "Frequency selected", 
       y = "Mean gain", color = "eta")
  
# feature stability by mean gain
ggplot(all_importances_df, aes(x = freq_selected, y = mean_cover, color = as.factor(eta))) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Feature importance stability across eta values", x = "Frequency selected", 
       y = "Mean cover", color = "eta")


### overfitting analysis - logloss over boosting rounds by eta values (stability of learning process across folds)
# loop through eta values to extract eval_logs
logloss_all <- lapply(names(all_eta_results), function(eta_name) {
  logs <- all_eta_results[[eta_name]]$eval_logs
  bind_rows(logs, .id = "Fold") %>%
    mutate(eta = eta_name)
}) %>% bind_rows()

# add max boosting rounds per eta to logloss
max_iter_by_eta <- sapply(all_eta_results, function(x) max(unlist(x$best_nrounds), na.rm = TRUE))
logloss_all <- logloss_all %>%
  mutate(max_iter = max_iter_by_eta[eta]) %>%
  filter(iter <= max_iter) # truncate to valid range of boosting rounds per eta

# compute the mean test logloss per eta and iteration
mean_logloss <- logloss_all %>%
  group_by(eta, iter) %>%
  summarise(mean_test_logloss = mean(test_logloss_mean, na.rm = TRUE), .groups = "drop")

# plot each fold's test logloss curve for each eta (250 folds each)
ggplot(logloss_all, aes(x = iter, y = test_logloss_mean, group = Fold)) +
  geom_line(alpha = 0.2, color = "black") +
  geom_line(data = mean_logloss, aes(x = iter, y = mean_test_logloss, group = NULL), color = "blue", linewidth = 1) + # mean curve
  facet_wrap(~ eta, scales = "free_x") + # use valid range of boosting rounds per eta
  labs(title = "Test logloss across folds for each eta",
       x = "Boosting round", y = "Test logloss")

# plot each fold's train logloss curve for each eta (250 folds each)
ggplot(logloss_all, aes(x = iter, y = train_logloss_mean, group = Fold)) +
  geom_line(alpha = 0.2, color = "red") +
  facet_wrap(~ eta, scales = "free_x") + # use valid range of boosting rounds per eta
  labs(title = "Train logloss across folds for each eta",
       x = "Boosting round", y = "Train logloss")

# plot each fold's test and train logloss curve for each eta 
# pivot to long format for logloss
logloss_long <- logloss_all %>%
  select(iter, Fold, eta, train_logloss_mean, test_logloss_mean) %>%
  pivot_longer(cols = c(train_logloss_mean, test_logloss_mean), names_to = "set", values_to = "logloss") %>%
  mutate(set = ifelse(set == "train_logloss_mean", "Train", "Test"))

# plot each fold's test and train logloss curve for each eta on same panel
ggplot(logloss_long, aes(x = iter, y = logloss, group = interaction(Fold, set), color = set)) +
  geom_line(alpha = 0.2) +
  facet_wrap(~ eta, scales = "free_x") +
  scale_color_manual(values = c("Train" = "red", "Test" = "black")) +
  labs(title = "Train and Test logloss across folds for each eta",
       x = "Boosting round", y = "Logloss", color = "Data")

# mean test versus train gap plot
# quantify generalization error over time per eta
logloss_gap <- logloss_long %>%
  pivot_wider(names_from = set, values_from = logloss) %>%
  mutate(gap = Test - Train) %>%
  group_by(eta, iter) %>%
  summarise(mean_gap = mean(gap, na.rm = TRUE), .groups = "drop")

ggplot(logloss_gap, aes(x = iter, y = mean_gap)) +
  geom_line() +
  facet_wrap(~ eta, scales = "free_x") +
  labs(title = "Mean logloss gap (test logloss - train logloss) across boosting rounds",
       y = "Gap", x = "Boosting round")


############################################################################
########   XGBOOST MODEL - FUNCTION TO EXPLORE HYPERPARAMETER SPACE  #######
############################################################################

### function to explore hyperparameter space with xgboost

tune_xgb_param <- function(param_grid_name,
                           param_grid_values,
                           base_params,
                           metagen,
                           all_feat_cols,
                           target_var,
                           target_var_numeric,
                           n_repeats = 50) {
  
  # set seed
  set.seed(1234)
  
  # n_repeats (default 50) repeats of stratified 5-fold cross-validation
  folds <- caret::createMultiFolds(metagen[[target_var]], k = 5, times = n_repeats)
  
  # list to store all results (metrics, feature importance, best nrounds, logloss) for each param grid value
  all_results <- list()
  
  # loop over val in param_grid_values
  for (val in param_grid_values) {
    cat("Running with", param_grid_name, "=", val, "\n")
    
    # update parameters with current val
    params <- base_params
    params[[param_grid_name]] <- val
    
    # create lists to store results 
    xgb_metrics <- list()
    xgb_importances <- list()
    best_nrounds_list <- list()
    eval_logs <- list()
    
    # loop over each of the cross-validation folds
    for (key in names(folds)) {
      
      # splits the dataset into training and testing sets for the current fold
      train_idx <- folds[[key]] # train indices 
      train_data <- metagen[train_idx, ] # training data
      test_data  <- metagen[-train_idx, ] # testing data
      
      # convert to DMatrix
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), 
                                     label = train_data[[target_var_numeric]])
      dtest  <- xgboost::xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), 
                                     label = test_data[[target_var_numeric]])
      
      # run internal 5-fold cross-validation with early stopping
      xgb_cv <- xgboost::xgb.cv(data = dtrain,
                                params = params,
                                nrounds = 5000,
                                nfold = 5,
                                early_stopping_rounds = 50,
                                maximize = FALSE,
                                stratified = TRUE,
                                showsd = FALSE,
                                verbose = 0)
      
      # best number of boosting rounds found by the internal 5-fold cross-validation above
      best_nrounds <- xgb_cv$best_iteration
      
      # store evaluation_og from xgb.cv()
      eval_log <- xgb_cv$evaluation_log
      eval_log$fold_id <- key
      eval_logs[[key]] <- eval_log
      
      # train final xgboost model with best nrounds
      xgb_model <- xgboost::xgboost(data = dtrain,
                                    params = params,
                                    nrounds = best_nrounds,
                                    verbose = 0)
      
      # evaluate on test set
      preds_prob <- predict(xgb_model, dtest)
      preds_label <- ifelse(preds_prob > 0.5, "disease", "healthy")
      preds_label <- factor(preds_label, levels = c("healthy", "disease"))
      
      # calculate AUC
      auc_val <- pROC::auc(response = test_data[[target_var]],
                           predictor = preds_prob,
                           levels = c("healthy", "disease"),
                           direction = "<")
      
      # generate confusion matrix
      cm <- caret::confusionMatrix(preds_label, test_data[[target_var]], positive = "disease")
      
      # calcualte logloss
      eps <- 1e-15 # avoid log(0)
      prob_clipped <- pmin(pmax(preds_prob, eps), 1 - eps)
      logloss <- -mean(test_data[[target_var_numeric]] * log(prob_clipped) +
                         (1 - test_data[[target_var_numeric]]) * log(1 - prob_clipped))
      
      # store performance metrics, importances, best nrounds and logloss
      xgb_metrics[[key]] <- list(cm = cm, auc = auc_val, logloss = logloss)
      imp <- xgboost::xgb.importance(feature_names = all_feat_cols, model = xgb_model)
      imp$Repeat_Fold <- key
      xgb_importances[[key]] <- imp
      best_nrounds_list[[key]] <- best_nrounds
    }
    
    ### summarize performance metrics
    # create vectors to store metrics
    balanced_accuracy <- f1_score <- precision <- sensitivity <- specificity <- auc_values <- logloss_values <- numeric()
    nrounds_values <- unlist(best_nrounds_list)
    
    # loop through stored lists to extract metrics
    for (key in names(xgb_metrics)) {
      cm <- xgb_metrics[[key]]$cm
      auc_val <- xgb_metrics[[key]]$auc
      logloss_val <- xgb_metrics[[key]]$logloss
      
      balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
      f1_score <- c(f1_score, cm$byClass["F1"])
      precision <- c(precision, cm$byClass["Precision"])
      sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
      specificity <- c(specificity, cm$byClass["Specificity"])
      auc_values <- c(auc_values, auc_val)
      logloss_values <- c(logloss_values, logloss_val)
    }
    
    # aggregate stored metrics into a data.frame
    summary_df <- data.frame(param = val,
                             mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(f1_score, na.rm = TRUE),
                             sd_f1 = sd(f1_score, na.rm = TRUE),
                             mean_precision = mean(precision, na.rm = TRUE),
                             sd_precision = sd(precision, na.rm = TRUE),
                             mean_sens = mean(sensitivity, na.rm = TRUE),
                             sd_sens = sd(sensitivity, na.rm = TRUE),
                             mean_spec = mean(specificity, na.rm = TRUE),
                             sd_spec = sd(specificity, na.rm = TRUE),
                             mean_auc = mean(auc_values, na.rm = TRUE),
                             sd_auc = sd(auc_values, na.rm = TRUE),
                             mean_logloss = mean(logloss_values, na.rm = TRUE),
                             sd_logloss = sd(logloss_values, na.rm = TRUE),
                             mean_nrounds = mean(nrounds_values, na.rm = TRUE),
                             sd_nrounds = sd(nrounds_values, na.rm = TRUE))
    
    # aggregate feature importances
    all_importances <- dplyr::bind_rows(xgb_importances)
    mean_importance <- all_importances %>%
      dplyr::group_by(Feature = Feature) %>%
      dplyr::summarise(mean_gain = mean(Gain, na.rm = TRUE),
                       mean_cover = mean(Cover, na.rm = TRUE),
                       mean_freq = mean(Frequency, na.rm = TRUE),
                       freq_selected = dplyr::n()) %>%
      dplyr::arrange(desc(mean_gain))
    
    # store all values under current param value name
    all_results[[paste0(param_grid_name, "_", val)]] <- list(summary = summary_df,
                                                             feature_importance = mean_importance,
                                                             raw_metrics = xgb_metrics,
                                                             best_nrounds = best_nrounds_list,
                                                             eval_logs = eval_logs)
  }
  
  return(all_results)
}


### arguments to include in the function
# name of parameter to tune
param_grid_name <- "eta"

# grid of param values to test
param_grid_values <- c(0.005, 0.01, 0.05, 0.1, 0.3)

# list of baseline hyperparameters to use (excluding param to test)
base_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = 0.3,
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

# metagen = data.frame of all features and labels (samples = rows, features + labels = columns)
metagen <- metagen

# predictor/feature columns (e.g., setdiff(colnames(metagen), c("condition", "condition_numeric")))
all_feat_cols <- all_feat_cols

# the label/variable as a factor, need to specify the order) (e.g., metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease")))
target_var <- "condition"

# the label/variable as a numeric (e.g., metagen$condition_numeric <- as.numeric(metagen$condition) - 1)
target_var_numeric <- "condition_numeric"

# number of repeats to use for cross-validation (default is 50)
n_repeats <- 50


### repeat eta exploration using function
param_grid_values <- c(0.005, 0.01, 0.05, 0.1, 0.3)
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
                                  metagen = metagen,
                                  all_feat_cols = all_feat_cols,
                                  target_var = "condition",
                                  target_var_numeric = "condition_numeric",
                                  n_repeats = 50)

### structure of output of tune_xgb_param
# all_results_param_grid_name
# ├── param_grid_name_param_grid_value
# │   ├── summary (performance metrics across folds)
# │   ├── feature_importance (average feature importances)
# │   ├── raw_metrics (raw confusion matrices, AUCs, logloss per fold)
# │   ├── best_nrounds (best round for each outer CV fold)
# │   └── eval_logs (xgboost evaluation logs for diagnostics)
# ├── param_grid_name_param_grid_value
# │   ├── summary
# │   ├── feature_importance
# │   ...


### functions to analyze the output of tune_xgb_param

### summarize performance metrics 
# uses ouput of tune_xgb_param
# mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds
summarize_performance <- function(results_list) {
  dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
}
summarize_performance(all_results_eta)

 
### feature importance (importance of features across parameter settings)
# uses ouput of tune_xgb_param, creates the all_importances_df
# combine feature importances across parameter values into a long-format table
# extracts numeric parameter value from list name (param_grid_values) and adds as column to each row
combine_feature_importance <- function(results_list, param_name = "eta") {
  purrr::imap_dfr(results_list, function(res, name) {
    df <- res$feature_importance
    df[[param_name]] <- as.numeric(gsub(paste0(param_name, "_"), "", name))
    df
  })
}
feature_freq_eta <- combine_feature_importance(all_results_eta, param_name = "eta")

# plot feature selection by frequency/selection across parameter values
# uses the all_importances_df
plot_feature_selection_frequency <- function(all_importances_df, feature, param_name = "eta") {
  filtered_df <- dplyr::filter(all_importances_df, Feature == feature)
  ggplot(filtered_df, aes(x = .data[[param_name]], y = freq_selected)) +
    geom_line() + geom_point() +
    labs(title = paste("Frequency of selection vs", param_name, "for", feature),
      x = param_name, y = "Frequency of selection")
}
plot_feature_selection_frequency(feature_freq_eta, feature = "Acutalibacter_muris", param_name = "eta")

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability <- function(all_importances_df, x = "freq_selected", y = "mean_gain", color_by = "eta") {
  ggplot(all_importances_df, aes(x = .data[[x]], y = .data[[y]], color = as.factor(.data[[color_by]]))) +
    geom_point(alpha = 0.6) + theme_minimal() +
    labs(title = "Feature importance stability", x = x, y = y, color = color_by)
}
plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_gain", color_by = "eta")
plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_cover", color_by = "eta")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table <- function(all_importances_df, threshold_frac = 0.4, n_repeats = 50, param_name = "eta") {
  threshold <- threshold_frac * n_repeats * 5 # 5 folds per repeat
  all_importances_df %>%
    group_by_at(param_name) %>%
    summarise(n_features_selected = sum(freq_selected >= threshold), .groups = "drop")
}
get_feature_stability_table(feature_freq_eta, threshold_frac = 0.4, n_repeats = 50, param_name = "eta")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency <- function(all_importances_df, param_name = "eta") {
  all_importances_df %>%
    group_by_at(param_name) %>%
    summarise(mean_frequency_selection = mean(freq_selected), .groups = "drop")
}
get_mean_feature_frequency(feature_freq_eta, param_name = "eta")


### logloss and overfitting analysis
# extract evaluation_logs (training/testing loss per boosting round) for all param_name values
# combines values into one table with Fold and param_name labeled
# uses output of tune_xgb_param, creates logloss_df
extract_logloss_df <- function(results_list, param_name = "eta") {
  purrr::imap_dfr(results_list, function(res, name) {
    logs <- res$eval_logs
    df <- dplyr::bind_rows(logs, .id = "Fold")
    df[[param_name]] <- name
    df
  })
}
logloss_eta <- extract_logloss_df(all_results_eta, param_name = "eta")

# plot logloss
# uses logloss_df
# plots change of logloss over boosting rounds (type = "train", "test", or "both")
plot_logloss_curve <- function(logloss_df, type = "test", show_mean = TRUE, param_name = "eta") {
  df <- logloss_df %>%
    mutate(!!param_name := as.numeric(gsub(paste0(param_name, "_"), "", .data[[param_name]])))
  
  p <- ggplot(df, aes(x = iter)) +
    facet_wrap(as.formula(paste("~", param_name)), scales = "free_x")
  
  if (type == "test") {
    p <- p + geom_line(aes(y = test_logloss_mean, group = Fold), alpha = 0.2, color = "black")
    if (show_mean) {
      mean_df <- df %>% group_by_at(c(param_name, "iter")) %>%
        summarise(mean_logloss = mean(test_logloss_mean, na.rm = TRUE), .groups = "drop")
      p <- p + geom_line(data = mean_df, aes(x = iter, y = mean_logloss), color = "blue", linewidth = 1)
    }
  } else if (type == "train") {
    p <- p + geom_line(aes(y = train_logloss_mean, group = Fold), alpha = 0.2, color = "red")
  } else if (type == "both") {
    df_long <- df %>%
      select(iter, Fold, all_of(param_name), train_logloss_mean, test_logloss_mean) %>%
      pivot_longer(cols = c(train_logloss_mean, test_logloss_mean), names_to = "set", values_to = "logloss")
    
    p <- ggplot(df_long, aes(x = iter, y = logloss, color = set, group = interaction(Fold, set))) +
      geom_line(alpha = 0.2) +
      scale_color_manual(values = c("train_logloss_mean" = "red", "test_logloss_mean" = "black")) +
      facet_wrap(as.formula(paste("~", param_name)), scales = "free_x")
  }
  
  p + labs(title = "Logloss curves", x = "Boosting round", y = "Logloss")
}
plot_logloss_curve(logloss_eta, type = "test", show_mean = TRUE, param_name = "eta")
plot_logloss_curve(logloss_eta, type = "train", show_mean = FALSE, param_name = "eta")
plot_logloss_curve(logloss_eta, type = "both", show_mean = FALSE, param_name = "eta")

# prepare logloss gap
# uses logloss_df, creates logloss_gap_df
# calculates the generalization gap (test loss - train loss) over boosting rounds and parameter values
prepare_logloss_gap <- function(logloss_df, param_name = "eta") {
  logloss_df %>%
    mutate(!!param_name := as.numeric(gsub(paste0(param_name, "_"), "", .data[[param_name]]))) %>%
    group_by_at(c(param_name, "iter")) %>%
    summarise(mean_train = mean(train_logloss_mean, na.rm = TRUE),
              mean_test = mean(test_logloss_mean, na.rm = TRUE),
              gap = mean(test_logloss_mean - train_logloss_mean, na.rm = TRUE),
              .groups = "drop")
}
logloss_gap_eta <- prepare_logloss_gap(logloss_eta, param_name = "eta")

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
plot_logloss_gap <- function(logloss_gap_df, param_name = "eta") {
  ggplot(logloss_gap_df, aes(x = iter, y = gap)) +
    geom_line() + theme_minimal() +
    facet_wrap(as.formula(paste("~", param_name)), scales = "free_x") +
    labs(title = paste("Mean logloss gap (test - train) across boosting rounds"),
      x = "Boosting round", y = "Gap")
}
plot_logloss_gap(logloss_gap_eta, param_name = "eta")


#############################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - CLASS WEIGHTS   #######
#############################################################################################


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
                                           metagen = metagen, 
                                           all_feat_cols = all_feat_cols,
                                           target_var = "condition",
                                           target_var_numeric = "condition_numeric",
                                           n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_class_weight)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_class_weight <- combine_feature_importance(all_results_class_weight, param_name = "scale_pos_weight")
feature_import_class_weight %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_class_weight, feature = "Acutalibacter_muris", param_name = "scale_pos_weight")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_class_weight, x = "freq_selected", y = "mean_gain", color_by = "scale_pos_weight")
plot_feature_stability(feature_import_class_weight, x = "freq_selected", y = "mean_cover", color_by = "scale_pos_weight")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_class_weight, threshold_frac = 0.4, n_repeats = 10, param_name = "scale_pos_weight")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_class_weight, param_name = "scale_pos_weight")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_class_weight <- extract_logloss_df(all_results_class_weight, param_name = "scale_pos_weight")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_class_weight, type = "test", show_mean = TRUE, param_name = "scale_pos_weight")
plot_logloss_curve(logloss_class_weight, type = "train", show_mean = FALSE, param_name = "scale_pos_weight")
plot_logloss_curve(logloss_class_weight, type = "both", show_mean = FALSE, param_name = "scale_pos_weight")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_class_weight <- prepare_logloss_gap(logloss_class_weight, param_name = "scale_pos_weight")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_class_weight, param_name = "scale_pos_weight")


#########################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MAX DEPTH   #######
#########################################################################################


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
                                        metagen = metagen, 
                                        all_feat_cols = all_feat_cols,
                                        target_var = "condition",
                                        target_var_numeric = "condition_numeric",
                                        n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_max_depth)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_max_depth <- combine_feature_importance(all_results_max_depth, param_name = "max_depth")
feature_import_max_depth %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_max_depth, feature = "Lachnoclostridium_sp._YL32", param_name = "max_depth")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_max_depth, x = "freq_selected", y = "mean_gain", color_by = "max_depth")
plot_feature_stability(feature_import_max_depth, x = "freq_selected", y = "mean_cover", color_by = "max_depth")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_max_depth, threshold_frac = 0.4, n_repeats = 10, param_name = "max_depth")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_max_depth, param_name = "max_depth")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_max_depth <- extract_logloss_df(all_results_max_depth, param_name = "max_depth")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_max_depth, type = "test", show_mean = TRUE, param_name = "max_depth")
plot_logloss_curve(logloss_max_depth, type = "train", show_mean = FALSE, param_name = "max_depth")
plot_logloss_curve(logloss_max_depth, type = "both", show_mean = FALSE, param_name = "max_depth")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_max_depth <- prepare_logloss_gap(logloss_max_depth, param_name = "max_depth")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_max_depth, param_name = "max_depth")


################################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MIN CHILD WEIGHT   #######
################################################################################################


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
                                               metagen = metagen, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_min_child_weight)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_min_child_weight <- combine_feature_importance(all_results_min_child_weight, param_name = "min_child_weight")
feature_import_min_child_weight %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_min_child_weight, feature = "Lachnoclostridium_sp._YL32", param_name = "min_child_weight")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_min_child_weight, x = "freq_selected", y = "mean_gain", color_by = "min_child_weight")
plot_feature_stability(feature_import_min_child_weight, x = "freq_selected", y = "mean_cover", color_by = "min_child_weight")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_min_child_weight, threshold_frac = 0.4, n_repeats = 10, param_name = "min_child_weight")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_min_child_weight, param_name = "min_child_weight")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_min_child_weight <- extract_logloss_df(all_results_min_child_weight, param_name = "min_child_weight")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_min_child_weight, type = "test", show_mean = FALSE, param_name = "min_child_weight")
plot_logloss_curve(logloss_min_child_weight, type = "train", show_mean = FALSE, param_name = "min_child_weight")
plot_logloss_curve(logloss_min_child_weight, type = "both", show_mean = FALSE, param_name = "min_child_weight")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_min_child_weight <- prepare_logloss_gap(logloss_min_child_weight, param_name = "min_child_weight")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_min_child_weight, param_name = "min_child_weight")


#########################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - SUBSAMPLE   #######
#########################################################################################


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
                                        metagen = metagen, 
                                        all_feat_cols = all_feat_cols,
                                        target_var = "condition",
                                        target_var_numeric = "condition_numeric",
                                        n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_subsample)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_subsample <- combine_feature_importance(all_results_subsample, param_name = "subsample")
feature_import_subsample %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_subsample, feature = "Lachnoclostridium_sp._YL32", param_name = "subsample")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_subsample, x = "freq_selected", y = "mean_gain", color_by = "subsample")
plot_feature_stability(feature_import_subsample, x = "freq_selected", y = "mean_cover", color_by = "subsample")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_subsample, threshold_frac = 0.4, n_repeats = 10, param_name = "subsample")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_subsample, param_name = "subsample")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_subsample <- extract_logloss_df(all_results_subsample, param_name = "subsample")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_subsample, type = "test", show_mean = TRUE, param_name = "subsample")
plot_logloss_curve(logloss_subsample, type = "train", show_mean = FALSE, param_name = "subsample")
plot_logloss_curve(logloss_subsample, type = "both", show_mean = FALSE, param_name = "subsample")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_subsample <- prepare_logloss_gap(logloss_subsample, param_name = "subsample")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_subsample, param_name = "subsample")


################################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - COLSAMPLE_BYTREE   #######
################################################################################################


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
                                               metagen = metagen, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_colsample_bytree)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_colsample_bytree <- combine_feature_importance(all_results_colsample_bytree, param_name = "colsample_bytree")
feature_import_colsample_bytree %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_colsample_bytree, feature = "Lachnoclostridium_sp._YL32", param_name = "colsample_bytree")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_colsample_bytree, x = "freq_selected", y = "mean_gain", color_by = "colsample_bytree")
plot_feature_stability(feature_import_colsample_bytree, x = "freq_selected", y = "mean_cover", color_by = "colsample_bytree")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_colsample_bytree, threshold_frac = 0.4, n_repeats = 10, param_name = "colsample_bytree")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_colsample_bytree, param_name = "colsample_bytree")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_colsample_bytree <- extract_logloss_df(all_results_colsample_bytree, param_name = "colsample_bytree")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_colsample_bytree, type = "test", show_mean = TRUE, param_name = "colsample_bytree")
plot_logloss_curve(logloss_colsample_bytree, type = "train", show_mean = FALSE, param_name = "colsample_bytree")
plot_logloss_curve(logloss_colsample_bytree, type = "both", show_mean = FALSE, param_name = "colsample_bytree")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_colsample_bytree <- prepare_logloss_gap(logloss_colsample_bytree, param_name = "colsample_bytree")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_colsample_bytree, param_name = "colsample_bytree")


################################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - COLSAMPLE_BYNODE   #######
################################################################################################


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
                                               metagen = metagen, 
                                               all_feat_cols = all_feat_cols,
                                               target_var = "condition",
                                               target_var_numeric = "condition_numeric",
                                               n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_colsample_bynode)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_colsample_bynode <- combine_feature_importance(all_results_colsample_bynode, param_name = "colsample_bynode")
feature_import_colsample_bynode %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_colsample_bynode, feature = "Lachnoclostridium_sp._YL32", param_name = "colsample_bynode")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_colsample_bynode, x = "freq_selected", y = "mean_gain", color_by = "colsample_bynode")
plot_feature_stability(feature_import_colsample_bynode, x = "freq_selected", y = "mean_cover", color_by = "colsample_bynode")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_colsample_bynode, threshold_frac = 0.4, n_repeats = 10, param_name = "colsample_bynode")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_colsample_bynode, param_name = "colsample_bynode")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_colsample_bynode <- extract_logloss_df(all_results_colsample_bynode, param_name = "colsample_bynode")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_colsample_bynode, type = "test", show_mean = TRUE, param_name = "colsample_bynode")
plot_logloss_curve(logloss_colsample_bynode, type = "train", show_mean = FALSE, param_name = "colsample_bynode")
plot_logloss_curve(logloss_colsample_bynode, type = "both", show_mean = FALSE, param_name = "colsample_bynode")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_colsample_bynode <- prepare_logloss_gap(logloss_colsample_bynode, param_name = "colsample_bynode")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_colsample_bynode, param_name = "colsample_bynode")


######################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - LAMBDA   #######
######################################################################################


### run xgboost model
param_grid_values <- c(0.1, 0.5, 1, 1.5, 5, 10)
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
                                    metagen = metagen, 
                                    all_feat_cols = all_feat_cols,
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric",
                                    n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_lambda)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_lambda <- combine_feature_importance(all_results_lambda, param_name = "lambda")
feature_import_lambda %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_lambda, feature = "Lachnoclostridium_sp._YL32", param_name = "lambda")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_lambda, x = "freq_selected", y = "mean_gain", color_by = "lambda")
plot_feature_stability(feature_import_lambda, x = "freq_selected", y = "mean_cover", color_by = "lambda")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_lambda, threshold_frac = 0.4, n_repeats = 10, param_name = "lambda")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_lambda, param_name = "lambda")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_lambda <- extract_logloss_df(all_results_lambda, param_name = "lambda")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_lambda, type = "test", show_mean = TRUE, param_name = "lambda")
plot_logloss_curve(logloss_lambda, type = "train", show_mean = FALSE, param_name = "lambda")
plot_logloss_curve(logloss_lambda, type = "both", show_mean = FALSE, param_name = "lambda")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_lambda <- prepare_logloss_gap(logloss_lambda, param_name = "lambda")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_lambda, param_name = "lambda") 


#####################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - ALPHA   #######
#####################################################################################


### run xgboost model
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 5)
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
                                    metagen = metagen, 
                                    all_feat_cols = all_feat_cols,
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric",
                                    n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_alpha)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_alpha <- combine_feature_importance(all_results_alpha, param_name = "alpha")
feature_import_alpha %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_alpha, feature = "Lachnoclostridium_sp._YL32", param_name = "alpha")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_alpha, x = "freq_selected", y = "mean_gain", color_by = "alpha")
plot_feature_stability(feature_import_alpha, x = "freq_selected", y = "mean_cover", color_by = "alpha")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_alpha, threshold_frac = 0.4, n_repeats = 10, param_name = "alpha")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_alpha, param_name = "alpha")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_alpha <- extract_logloss_df(all_results_alpha, param_name = "alpha")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_alpha, type = "test", show_mean = TRUE, param_name = "alpha")
plot_logloss_curve(logloss_alpha, type = "train", show_mean = FALSE, param_name = "alpha")
plot_logloss_curve(logloss_alpha, type = "both", show_mean = FALSE, param_name = "alpha")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_alpha <- prepare_logloss_gap(logloss_alpha, param_name = "alpha")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_alpha, param_name = "alpha") 


#####################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - GAMMA   #######
#####################################################################################


### run xgboost model
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 5)
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
                                    metagen = metagen, 
                                    all_feat_cols = all_feat_cols,
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric",
                                    n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_gamma)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_gamma <- combine_feature_importance(all_results_gamma, param_name = "gamma")
feature_import_gamma %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_gamma, feature = "Lachnoclostridium_sp._YL32", param_name = "gamma")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_gamma, x = "freq_selected", y = "mean_gain", color_by = "gamma")
plot_feature_stability(feature_import_gamma, x = "freq_selected", y = "mean_cover", color_by = "gamma")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_gamma, threshold_frac = 0.4, n_repeats = 10, param_name = "gamma")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_gamma, param_name = "gamma")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_gamma <- extract_logloss_df(all_results_gamma, param_name = "gamma")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_gamma, type = "test", show_mean = TRUE, param_name = "gamma")
plot_logloss_curve(logloss_gamma, type = "train", show_mean = FALSE, param_name = "gamma")
plot_logloss_curve(logloss_gamma, type = "both", show_mean = FALSE, param_name = "gamma")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_gamma <- prepare_logloss_gap(logloss_gamma, param_name = "gamma")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_gamma, param_name = "gamma") 


##############################################################################################
########   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MAX_DELTA_STEP   #######
##############################################################################################


### run xgboost model
param_grid_values <- c(0, 1, 2, 3, 4, 5)
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
                                             metagen = metagen, 
                                             all_feat_cols = all_feat_cols,
                                             target_var = "condition",
                                             target_var_numeric = "condition_numeric",
                                             n_repeats = 10)


### summarize performance metrics
# uses ouput of tune_xgb_param
summarize_performance(all_results_max_delta_step)


### feature importance
# combine feature importances
# uses ouput of tune_xgb_param, creates the all_importances_df
feature_import_max_delta_step <- combine_feature_importance(all_results_max_delta_step, param_name = "max_delta_step")
feature_import_max_delta_step %>% arrange(desc(freq_selected))

# plot specific feature by frequency across parameter values
# uses the all_importances_df
plot_feature_selection_frequency(feature_import_max_delta_step, feature = "Lachnoclostridium_sp._YL32", param_name = "max_delta_step")

# plot feature frequency importance (mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability(feature_import_max_delta_step, x = "freq_selected", y = "mean_gain", color_by = "max_delta_step")
plot_feature_stability(feature_import_max_delta_step, x = "freq_selected", y = "mean_cover", color_by = "max_delta_step")

# number of features selected in more than threshold_frac models
# uses all_importances_df
get_feature_stability_table(feature_import_max_delta_step, threshold_frac = 0.4, n_repeats = 10, param_name = "max_delta_step")

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency(feature_import_max_delta_step, param_name = "max_delta_step")


### logloss and overfitting analysis
# extract evaluation_logs and combine into one table
# uses tune_xgb_param, creates logloss_df
logloss_max_delta_step <- extract_logloss_df(all_results_max_delta_step, param_name = "max_delta_step")

# plot logloss over boosting rounds (type = "train", "test", or "both")
# uses logloss_df
plot_logloss_curve(logloss_max_delta_step, type = "test", show_mean = TRUE, param_name = "max_delta_step")
plot_logloss_curve(logloss_max_delta_step, type = "train", show_mean = FALSE, param_name = "max_delta_step")
plot_logloss_curve(logloss_max_delta_step, type = "both", show_mean = FALSE, param_name = "max_delta_step")

# calculate logloss gap (test loss - train loss) over boosting rounds
# uses logloss_df, creates logloss_gap_df
logloss_gap_max_delta_step <- prepare_logloss_gap(logloss_max_delta_step, param_name = "max_delta_step")

# plot logloss gap
# uses logloss_gap_df
plot_logloss_gap(logloss_gap_max_delta_step, param_name = "max_delta_step") 


###########################################################################################################
########   XGBOOST LOGLOSS MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS WITHOUT PARALLELIZATON  #######
###########################################################################################################

scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step) {
  
  set.seed(1234)
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(metagen[, all_feat_cols]), 
                        label = metagen$condition_numeric)
  
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
    folds <- caret::createFolds(metagen$condition, k = 5, list = TRUE, returnTrain = FALSE)
    
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
bounds <- list(eta = c(0.001, 0.02),
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

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 12,
                   iters.n = 10,
                   acq = "ei",
                   parallel = FALSE,
                   verbose = 1)

# view results
optObj$stopStatus
head(optObj$scoreSummary[order(-Score), ])
getBestPars(optObj)


###################################################################################################################
########   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF CV FOLDS - LOGLOSS  #######
###################################################################################################################

library(future.apply) # allows parallel versions of *apply functions
plan(multisession, workers = 5)  # each worker in a separate R session, one worker per 2 repeats

# data and feat_cols explicitly passed so each worker has a complete copy of data (avoid scope issue in parallel execution)
scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step,
                             data, feat_cols) {
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(data[, feat_cols]), 
                        label = data$condition_numeric)
  
  params <- list(objective = "binary:logistic",
                 eval_metric = "logloss",
                 eta = eta,
                 max_depth = as.integer(max_depth), # integer
                 min_child_weight = as.integer(min_child_weight), # integer
                 subsample = subsample,
                 colsample_bytree = colsample_bytree,
                 colsample_bynode = colsample_bynode,
                 gamma = gamma,
                 lambda = lambda,
                 alpha = alpha,
                 scale_pos_weight = scale_pos_weight,
                 max_delta_step = as.integer(max_delta_step)) # integer
  
  # number of repeats of stratified 5-fold cross-validation
  repeats <- 10 
  
  # runs the anonymous function in parallel once per repeat, parallelized across the workers
  # output collected by repeat_logloss
  # data and feat_cols explicitly passed to each worker
  repeat_logloss <- future_sapply(1:repeats,
                                  function(r, data_local, feat_cols_local) {
                                    
                                    # gives each repeat a different random seed for fold generation (prevents identical scores and zero variance for Bayesian optimizaiton)
                                    set.seed(sample.int(1e6, 1))
                                    
                                    # create 5-folds for cross-validation
                                    folds <- caret::createFolds(data_local$condition, k = 5, list = TRUE, returnTrain = FALSE)
                                    
                                    # ensures that any xgboost error returns NULL instead of crashing
                                    xgb_cv <- tryCatch({
                                      xgboost::xgb.cv(params = params,
                                                      data = xgb.DMatrix(data = as.matrix(data_local[, feat_cols_local]), 
                                                                         label = data_local$condition_numeric),
                                                      nrounds = 5000,
                                                      folds = folds,
                                                      stratified = TRUE,
                                                      showsd = FALSE,
                                                      early_stopping_rounds = 50,
                                                      verbose = 0)
                                      }, error = function(e) NULL)
                                    
                                    # replaces failures with a finite large number (prevents NA/Inf issues that will crash BayesOpt)
                                    if (is.null(xgb_cv)) return(Inf)
                                    
                                    # returns the best iteration of logloss for the repeat
                                    xgb_cv$evaluation_log$test_logloss_mean[xgb_cv$best_iteration]
                                    },
                                  
                                  future.seed = TRUE, # ensures random number generation is reproducible across parallel workers
                                  data_local = data,
                                  feat_cols_local = feat_cols) # pass dataset and feature columns from parent environment to each worker
  
  # average across repeats
  mean_logloss <- mean(repeat_logloss)
  
  # negative because Bayesian optimization maximizes Score (converts minimization into a maximization problem)
  return(list(Score = -mean_logloss))
}

# define parameter bounds
bounds <- list(eta = c(0.001, 0.02),
               max_depth = c(2L, 6L), # integer
               min_child_weight = c(1L, 2L), # integer
               subsample = c(0.9, 1.0),
               colsample_bytree = c(0.6, 1.0),
               colsample_bynode = c(0.6, 1.0),
               gamma = c(0, 4),
               lambda = c(0, 10),
               alpha = c(0, 4),
               scale_pos_weight = c(0.3, 2.0),
               max_delta_step = c(0L, 10L)) # integer

# run Bayesian optimization sequentially, letting repeats run in parallel
set.seed(1234) # reproducible search of initial points
optObj <- bayesOpt(FUN = function(...) scoring_function(..., data = metagen, feat_cols = all_feat_cols), # wrapper to pass dataset and feature columns to scoring_fucntion
                   bounds = bounds,
                   initPoints = 22,
                   iters.n = 50,
                   acq = "ei",
                   parallel = FALSE,
                   verbose = 1)

# view results
stopstatus_cv_log = optObj$stopStatus
scoresum_cv_log = optObj$scoreSummary[order(-Score), ]
head(scoresum_cv_log)
bestparams_cv_log = getBestPars(optObj)

plan(sequential) # resets future plan to normal single-threaded execution


###############################################################################################################
########   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF CV FOLDS - AUC  #######
###############################################################################################################

library(future.apply) # allows parallel versions of *apply functions
plan(multisession, workers = 5)  # each worker in a separate R session, one worker per 2 repeats

# data and feat_cols explicitly passed so each worker has a complete copy of data (avoid scope issue in parallel execution)
scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step,
                             data, feat_cols) {
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(data[, feat_cols]), 
                        label = data$condition_numeric)
  
  params <- list(objective = "binary:logistic",
                 eval_metric = "auc",
                 eta = eta,
                 max_depth = as.integer(max_depth), # integer
                 min_child_weight = as.integer(min_child_weight), # integer
                 subsample = subsample,
                 colsample_bytree = colsample_bytree,
                 colsample_bynode = colsample_bynode,
                 gamma = gamma,
                 lambda = lambda,
                 alpha = alpha,
                 scale_pos_weight = scale_pos_weight,
                 max_delta_step = as.integer(max_delta_step)) # integer
  
  # number of repeats of stratified 5-fold cross-validation
  repeats <- 10 
  
  # runs the anonymous function in parallel once per repeat, parallelized across the workers
  # output collected by repeat_logloss
  # data and feat_cols explicitly passed to each worker
  repeat_auc <- future_sapply(1:repeats,
                              function(r, data_local, feat_cols_local) {
                                
                                # gives each repeat a different random seed for fold generation (prevents identical scores and zero variance for Bayesian optimizaiton)
                                set.seed(sample.int(1e6, 1))
                                
                                # create 5-folds for cross-validation
                                folds <- caret::createFolds(data_local$condition, k = 5, list = TRUE, returnTrain = FALSE)
                                
                                # ensures that any xgboost error returns NULL instead of crashing
                                xgb_cv <- tryCatch({
                                  xgboost::xgb.cv(params = params,
                                                  data = xgb.DMatrix(data = as.matrix(data_local[, feat_cols_local]), 
                                                                     label = data_local$condition_numeric),
                                                  nrounds = 5000,
                                                  folds = folds,
                                                  stratified = TRUE,
                                                  showsd = FALSE,
                                                  early_stopping_rounds = 50,
                                                  verbose = 0)
                                }, error = function(e) NULL)
                                
                                # replaces failures with a finite large number (prevents NA/Inf issues that will crash BayesOpt)
                                if (is.null(xgb_cv)) return(0)
                                
                                # returns the best iteration of auc for the repeat
                                xgb_cv$evaluation_log$test_auc_mean[xgb_cv$best_iteration]
                              },
                              
                              future.seed = TRUE, # ensures random number generation is reproducible across parallel workers
                              data_local = data,
                              feat_cols_local = feat_cols) # pass dataset and feature columns from parent environment to each worker
  
  # average across repeats
  mean_auc <- mean(repeat_auc)
  
  # return mean auc values
  return(list(Score = mean_auc))
}

# define parameter bounds
bounds <- list(eta = c(0.001, 0.02),
               max_depth = c(2L, 6L), # integer
               min_child_weight = c(1L, 2L), # integer
               subsample = c(0.9, 1.0),
               colsample_bytree = c(0.6, 1.0),
               colsample_bynode = c(0.6, 1.0),
               gamma = c(0, 4),
               lambda = c(0, 10),
               alpha = c(0, 4),
               scale_pos_weight = c(0.3, 2.0),
               max_delta_step = c(0L, 10L)) # integer

# run Bayesian optimization sequentially, letting repeats run in parallel
set.seed(1234) # reproducible search of initial points
optObj <- bayesOpt(FUN = function(...) scoring_function(..., data = metagen, feat_cols = all_feat_cols), # wrapper to pass dataset and feature columns to scoring_fucntion
                   bounds = bounds,
                   initPoints = 22,
                   iters.n = 50,
                   acq = "ei",
                   parallel = FALSE,
                   verbose = 1)

# view results
stopstatus_cv_auc = optObj$stopStatus
scoresum_cv_auc = optObj$scoreSummary[order(-Score), ]
head(scoresum_cv_auc)
bestparams_cv_auc = getBestPars(optObj)

plan(sequential) # resets future plan to normal single-threaded execution


####################################################################################################################
########   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF BAYES OPT - LOGLOSS  #######
####################################################################################################################

scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step) {
  
  set.seed(1234)
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(metagen[, all_feat_cols]), 
                        label = metagen$condition_numeric)
  
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
    folds <- caret::createFolds(metagen$condition, k = 5, list = TRUE, returnTrain = FALSE)
    
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
bounds <- list(eta = c(0.001, 0.02),
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
# each hyperparameter point being evaluated by Bayesain optimization can run on a separate core
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 22,
                   iters.n = 50,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ() 

# view results
stopstatus_bo_log = optObj$stopStatus
scoresum_bo_log = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
bestparams_bo_log = getBestPars(optObj)


################################################################################################################
########   XGBOOST MODEL - BAYESIAN OPTIMIZATION OF HYPERPARAMETERS - PARALLELIZATON OF BAYES OPT - AUC  #######
################################################################################################################

scoring_function <- function(eta, max_depth, min_child_weight, subsample,
                             colsample_bytree, colsample_bynode, gamma,
                             lambda, alpha, scale_pos_weight, max_delta_step) {
  
  set.seed(1234)
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(metagen[, all_feat_cols]), 
                        label = metagen$condition_numeric)
  
  params <- list(objective = "binary:logistic",
                 eval_metric = "auc",
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
  repeat_auc <- numeric(repeats)
  
  # loop over each repeat
  for (r in 1:repeats) {
    
    #  create 5-folds for cross-validation (stratified on condition)
    folds <- caret::createFolds(metagen$condition, k = 5, list = TRUE, returnTrain = FALSE)
    
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
      repeat_auc[r] <- 0
    } else {
      
      # take the best iteration of auc for the repeat
      repeat_auc[r] <- xgb_cv$evaluation_log$test_auc_mean[xgb_cv$best_iteration]
    }
  }
  
  # average across repeats
  mean_auc <- mean(repeat_auc)
  
  # return mean auc values
  return(list(Score = mean_auc))
}

# define parameter bounds
bounds <- list(eta = c(0.001, 0.02),
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
# each hyperparameter point being evaluated by Bayesain optimization can run on a separate core
doParallel::registerDoParallel(parallel::detectCores() - 1)

set.seed(1234)
optObj <- bayesOpt(FUN = scoring_function,
                   bounds = bounds,
                   initPoints = 22,
                   iters.n = 50,
                   acq = "ei",
                   parallel = TRUE,
                   verbose = 1)

# unregister the backend
registerDoSEQ() 

# view results
stopstatus_bo_auc = optObj$stopStatus
scoresum_bo_auc = optObj$scoreSummary[order(-Score), ]
head(optObj$scoreSummary[order(-Score), ])
bestparams_bo_auc = getBestPars(optObj)


############################################################################################################
########   XGBOOST MODEL - EVALUATION OF MODEL WITH BEST HYPERPARAMETERS USING LOGLOSS - BO LOGLOSS  #######
############################################################################################################

### function to evaluate model using hyperparameter values
final_xgb_evaluation <- function(best_params,
                                 metagen,
                                 all_feat_cols,
                                 target_var,
                                 target_var_numeric,
                                 n_repeats = 50,
                                 n_folds = 5) {
  
  # set seed for reproducibility
  set.seed(1234)
  
  # lists to store all results (metrics, feature importance, best nrounds, eval_logs)
  all_results <- list()
  
  # create lists to store results 
  xgb_metrics <- list()
  xgb_importances <- list()
  best_nrounds_list <- list()
  eval_logs <- list()
  
  # n_repeats of n_folds cross-validation
  for (r in 1:n_repeats) {
    cat("Repeat:", r, "\n")
    
    # different seed per repeat to get different folds
    set.seed(1234 + r)
    folds <- caret::createFolds(metagen[[target_var]], k = n_folds, list = TRUE, returnTrain = FALSE)
    
    # loop over each cross-validation fold
    for (key in names(folds)) {
      
      # splits the dataset into training and testing sets for the current fold
      test_idx <- folds[[key]] 
      test_data <- metagen[test_idx, ]
      train_data <- metagen[-test_idx, ]
      
      
      # convert to DMatrix
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), 
                                     label = train_data[[target_var_numeric]])
      dtest <- xgboost::xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), 
                                     label = test_data[[target_var_numeric]])
      
      # run internal 5-fold cross-validation to determine best nrounds
      xgb_cv <- xgboost::xgb.cv(data = dtrain,
                                params = best_params,
                                nrounds = 5000,
                                nfold = 5,
                                early_stopping_rounds = 50,
                                maximize = FALSE,
                                stratified = TRUE,
                                showsd = FALSE,
                                verbose = 0)
      
      # best number of boosting rounds found by the internal 5-fold cross-validation above
      best_nrounds <- xgb_cv$best_iteration
      
      # store evaluation_log from xgb.cv()
      eval_log <- xgb_cv$evaluation_log
      eval_log$fold_id <- key
      eval_logs[[paste0("R", r, "_F", key)]] <- eval_log
      
      # train final model on training data using best nrounds
      xgb_model <- xgboost::xgboost(data = dtrain,
                                    params = best_params,
                                    nrounds = best_nrounds,
                                    verbose = 0)
      
      # evaluate on test set
      preds_prob <- predict(xgb_model, dtest)
      preds_label <- ifelse(preds_prob > 0.5, "disease", "healthy")
      preds_label <- factor(preds_label, levels = c("healthy", "disease"))
      
      # calculate AUC
      auc_val <- pROC::auc(response = test_data[[target_var]],
                           predictor = preds_prob,
                           levels = c("healthy", "disease"),
                           direction = "<")
      
      # generate confusion matrix
      cm <- caret::confusionMatrix(preds_label, test_data[[target_var]], positive = "disease")
      
      # calculate logloss
      eps <- 1e-15 # avoid log(0)
      prob_clipped <- pmin(pmax(preds_prob, eps), 1 - eps)
      logloss <- -mean(test_data[[target_var_numeric]] * log(prob_clipped) +
                         (1 - test_data[[target_var_numeric]]) * log(1 - prob_clipped))
      
      # store performance metrics, importance, best_nrounds and logloss
      xgb_metrics[[paste0("R", r, "_F", key)]] <- list(cm = cm, auc = auc_val, logloss = logloss)
      imp <- xgboost::xgb.importance(feature_names = all_feat_cols, model = xgb_model)
      imp$Repeat_Fold <- paste0("R", r, "_F", key)
      xgb_importances[[paste0("R", r, "_F", key)]] <- imp
      best_nrounds_list[[paste0("R", r, "_F", key)]] <- best_nrounds
    }
  }
  
  ### summarize performance metrics
  # create vectors to store metrics
  balanced_accuracy <- f1_score <- precision <- sensitivity <- specificity <- auc_values <- logloss_values <- numeric()
  nrounds_values <- unlist(best_nrounds_list)
  
  # loop through stored lists to extract metrics
  for (key in names(xgb_metrics)) {
    cm <- xgb_metrics[[key]]$cm
    auc_val <- xgb_metrics[[key]]$auc
    logloss_val <- xgb_metrics[[key]]$logloss
    
    balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
    f1_score <- c(f1_score, cm$byClass["F1"])
    precision <- c(precision, cm$byClass["Precision"])
    sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
    specificity <- c(specificity, cm$byClass["Specificity"])
    auc_values <- c(auc_values, auc_val)
    logloss_values <- c(logloss_values, logloss_val)
  }
  
  # aggregate stored metrics into a data.frame
  summary_df <- data.frame(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
                           sd_bal_acc = sd(balanced_accuracy, na.rm = TRUE),
                           mean_f1 = mean(f1_score, na.rm = TRUE),
                           sd_f1 = sd(f1_score, na.rm = TRUE),
                           mean_precision = mean(precision, na.rm = TRUE),
                           sd_precision = sd(precision, na.rm = TRUE),
                           mean_sens = mean(sensitivity, na.rm = TRUE),
                           sd_sens = sd(sensitivity, na.rm = TRUE),
                           mean_spec = mean(specificity, na.rm = TRUE),
                           sd_spec = sd(specificity, na.rm = TRUE),
                           mean_auc = mean(auc_values, na.rm = TRUE),
                           sd_auc = sd(auc_values, na.rm = TRUE),
                           mean_logloss = mean(logloss_values, na.rm = TRUE),
                           sd_logloss = sd(logloss_values, na.rm = TRUE),
                           mean_nrounds = mean(nrounds_values, na.rm = TRUE),
                           sd_nrounds = sd(nrounds_values, na.rm = TRUE))
  
  # aggregate feature importances
  all_importances <- dplyr::bind_rows(xgb_importances)
  mean_importance <- all_importances %>%
    dplyr::group_by(Feature = Feature) %>%
    dplyr::summarise(mean_gain = mean(Gain, na.rm = TRUE),
                     mean_cover = mean(Cover, na.rm = TRUE),
                     mean_freq = mean(Frequency, na.rm = TRUE),
                     freq_selected = dplyr::n()) %>%
    dplyr::arrange(desc(mean_gain))
  
  # store all values under current param value name
  all_results[["final_evaluation"]] <- list(summary = summary_df,
                                            feature_importance = mean_importance,
                                            raw_metrics = xgb_metrics,
                                            best_nrounds = best_nrounds_list,
                                            eval_logs = eval_logs)
  
  return(all_results)
}

# list of best hyperparameters to use (from 50 repeats of 5-fold cross-validation parallelized on the cv folds)
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = bestparams_bo_log$eta,
                    scale_pos_weight = bestparams_bo_log$scale_pos_weight,
                    max_depth = bestparams_bo_log$max_depth,
                    min_child_weight = bestparams_bo_log$min_child_weight,
                    subsample = bestparams_bo_log$subsample,
                    colsample_bytree = bestparams_bo_log$colsample_bytree,
                    colsample_bynode = bestparams_bo_log$colsample_bynode,
                    lambda = bestparams_bo_log$lambda,
                    alpha = bestparams_bo_log$alpha,
                    gamma = bestparams_bo_log$alpha,
                    max_delta_step = bestparams_bo_log$max_delta_step)

# metagen = data.frame of all features and labels (samples = rows, features + labels = columns)
metagen <- metagen

# predictor/feature columns (e.g., setdiff(colnames(metagen), c("condition", "condition_numeric")))
all_feat_cols <- all_feat_cols

# the label/variable as a factor, need to specify the order) (e.g., metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease")))
target_var <- "condition"

# the label/variable as a numeric (e.g., metagen$condition_numeric <- as.numeric(metagen$condition) - 1)
target_var_numeric <- "condition_numeric"


### run final evaluation
final_results_log <- final_xgb_evaluation(best_params = best_params,
                                          metagen = metagen,
                                          all_feat_cols = all_feat_cols,
                                          target_var = "condition",
                                          target_var_numeric = "condition_numeric",
                                          n_repeats = 50,
                                          n_folds = 5)
  
### functions to analyze the output of final_xgb_evaluation

### summarize performance metrics (same as for tune_xgb_param)
# uses ouput of final_xgb_evaluation
# mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds
summarize_performance <- function(results_list) {
  dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
}
summarize_performance(final_results_log)


### feature importance (importance of features across parameter settings)
# uses ouput of final_xgb_evaluation, creates the all_importances_df
final_feature_importance <- function(results_list) {
  purrr::imap_dfr(results_list, function(res, name) {
    df <- res$feature_importance
    return(df)
  })
}
feature_freq_final_log <- final_feature_importance(final_results_log)
feature_freq_final_log

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability_final <- function(all_importances_df, x = "freq_selected", y = "mean_gain") {
  ggplot(all_importances_df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.6, color = "blue") + theme_minimal() +
    labs(title = "Feature importance stability", x = x, y = y)
} 
plot_feature_stability_final(feature_freq_final_log, x = "freq_selected", y = "mean_gain")
plot_feature_stability_final(feature_freq_final_log, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
get_feature_stability_table_final <- function(all_importances_df, threshold_frac = 0.4, n_repeats = 50) {
  threshold <- threshold_frac * n_repeats * 5 # 5 folds per repeat
  all_importances_df %>%
    summarise(n_features_selected = sum(freq_selected >= threshold), .groups = "drop")
}
get_feature_stability_table_final(feature_freq_final_log, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency_final <- function(all_importances_df) {
  all_importances_df %>%
    summarise(mean_frequency_selection = mean(freq_selected), .groups = "drop")
}
get_mean_feature_frequency_final(feature_freq_final_log)


### logloss and overfitting analysis
# extract evaluation_logs (training/testing loss per boosting round)
# uses output of final_xgb_evaluation, creates logloss_df
extract_logloss_df_final <- function(results_list) {
  purrr::imap_dfr(results_list, function(res, name) {
    logs <- res$eval_logs
    df <- dplyr::bind_rows(logs, .id = "Fold")
    return(df)
  })
}
logloss_final_log <- extract_logloss_df_final(final_results_log)
logloss_final_log

# plot logloss
# uses logloss_df
# plots change of logloss over boosting rounds (type = "train", "test", or "both")
plot_logloss_curve_final <- function(logloss_df, type = "test", show_mean = TRUE) {
  p <- ggplot(logloss_df, aes(x = iter))
  
  if (type == "test") {
    p <- p + geom_line(aes(y = test_logloss_mean, group = Fold), alpha = 0.2, color = "black")
    
    if (show_mean) {
      mean_df <- logloss_df %>%
        group_by(iter) %>%
        summarise(mean_logloss = mean(test_logloss_mean, na.rm = TRUE), .groups = "drop")
      
      p <- p + geom_line(data = mean_df, aes(x = iter, y = mean_logloss),
                         color = "blue", linewidth = 1)
    }
    
  } else if (type == "train") {
    p <- p + geom_line(aes(y = train_logloss_mean, group = Fold), alpha = 0.2, color = "red")
    
  } else if (type == "both") {
    df_long <- logloss_df %>%
      select(iter, Fold, train_logloss_mean, test_logloss_mean) %>%
      tidyr::pivot_longer(cols = c(train_logloss_mean, test_logloss_mean),
                          names_to = "set", values_to = "logloss")
    
    p <- ggplot(df_long, aes(x = iter, y = logloss, color = set, group = interaction(Fold, set))) +
      geom_line(alpha = 0.2) +
      scale_color_manual(values = c("train_logloss_mean" = "red",
                                    "test_logloss_mean" = "black"))
  }
  
  p + labs(title = "Logloss curves", x = "Boosting round", y = "Logloss")
}
plot_logloss_curve_final(logloss_final_log, type = "test", show_mean = TRUE)
plot_logloss_curve_final(logloss_final_log, type = "train", show_mean = FALSE)
plot_logloss_curve_final(logloss_final_log, type = "both", show_mean = FALSE)

# prepare logloss gap
# uses logloss_df, creates logloss_gap_df
# calculates the generalization gap (test loss - train loss) over boosting rounds
prepare_logloss_gap_final <- function(logloss_df) {
  logloss_df %>%
    group_by(iter) %>%
    summarise(
      mean_train = mean(train_logloss_mean, na.rm = TRUE),
      mean_test = mean(test_logloss_mean, na.rm = TRUE),
      gap = mean(test_logloss_mean - train_logloss_mean, na.rm = TRUE),
      .groups = "drop"
    )
}
logloss_gap_final_log <- prepare_logloss_gap_final(logloss_final_log)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
plot_logloss_gap_final <- function(logloss_gap_df) {
  ggplot(logloss_gap_df, aes(x = iter, y = gap)) +
    geom_line(color = "darkorange", linewidth = 1) +
    theme_minimal() +
    labs(
      title = "Mean logloss gap (test - train) across boosting rounds",
      x = "Boosting round",
      y = "Logloss gap"
    )
}

plot_logloss_gap_final(logloss_gap_final_log)


########################################################################################################
########   XGBOOST MODEL - EVALUATION OF MODEL WITH BEST HYPERPARAMETERS USING LOGLOSS - BO AUC  #######
########################################################################################################

# list of best hyperparameters 
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = bestparams_bo_auc$eta,
                    scale_pos_weight = bestparams_bo_auc$scale_pos_weight,
                    max_depth = bestparams_bo_auc$max_depth,
                    min_child_weight = bestparams_bo_auc$min_child_weight,
                    subsample = bestparams_bo_auc$subsample,
                    colsample_bytree = bestparams_bo_auc$colsample_bytree,
                    colsample_bynode = bestparams_bo_auc$colsample_bynode,
                    lambda = bestparams_bo_auc$lambda,
                    alpha = bestparams_bo_auc$alpha,
                    gamma = bestparams_bo_auc$alpha,
                    max_delta_step = bestparams_bo_auc$max_delta_step)

# metagen = data.frame of all features and labels (samples = rows, features + labels = columns)
metagen <- metagen

# predictor/feature columns (e.g., setdiff(colnames(metagen), c("condition", "condition_numeric")))
all_feat_cols <- all_feat_cols

# the label/variable as a factor, need to specify the order) (e.g., metagen$condition <- factor(metagen$condition, levels = c("healthy", "disease")))
target_var <- "condition"

# the label/variable as a numeric (e.g., metagen$condition_numeric <- as.numeric(metagen$condition) - 1)
target_var_numeric <- "condition_numeric"

### evaluate model using hyperparameter values determined using Bayesian optimization (AUC)
final_results_auc <- final_xgb_evaluation(best_params = best_params,
                                          metagen = metagen,
                                          all_feat_cols = all_feat_cols,
                                          target_var = "condition",
                                          target_var_numeric = "condition_numeric",
                                          n_repeats = 50,
                                          n_folds = 5)

### functions to analyze the output of final_xgb_evaluation

### summarize performance metrics
# uses ouput of final_xgb_evaluation
# mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds
summarize_performance(final_results_auc)


### feature importance (importance of features across parameter settings)
# uses output of final_xgb_evaluation, creates the all_importances_df
feature_freq_final_auc <- final_feature_importance(final_results_auc)
feature_freq_final_auc

# plot feature frequency/selection by importance(mean gain or mean cover)
# uses the all_importances_df
plot_feature_stability_final(feature_freq_final_auc, x = "freq_selected", y = "mean_gain")
plot_feature_stability_final(feature_freq_final_auc, x = "freq_selected", y = "mean_cover")

# number of features selected in more than threshold_frac folds
# uses all_importances_df
get_feature_stability_table_final(feature_freq_final_auc, threshold_frac = 0.4, n_repeats = 50)

# mean frequency of feature selection per parameter value
# uses all_importances_df
get_mean_feature_frequency_final(feature_freq_final_auc)


### logloss and overfitting analysis
# extract evaluation_logs (training/testing loss per boosting round)
# uses output of final_xgb_evaluation, creates logloss_df
logloss_final_auc <- extract_logloss_df_final(final_results_auc)
logloss_final_auc

# plot logloss
# uses logloss_df
# plots change of logloss over boosting rounds (type = "train", "test", or "both")
plot_logloss_curve_final(logloss_final_auc, type = "test", show_mean = TRUE)
plot_logloss_curve_final(logloss_final_auc, type = "train", show_mean = FALSE)
plot_logloss_curve_final(logloss_final_auc, type = "both", show_mean = FALSE)

# prepare logloss gap
# uses logloss_df, creates logloss_gap_df
# calculates the generalization gap (test loss - train loss) over boosting rounds
logloss_gap_final_auc <- prepare_logloss_gap_final(logloss_final_auc)

# plot generalization gap (visually assess overfitting)
# uses logloss_gap_df
plot_logloss_gap_final(logloss_gap_final_auc)


#######################################################################################
########   XGBOOST LOGLOSS MODEL - TRAIN FINAL MODEL WITH BEST HYPERPARAMETERS  #######
#######################################################################################

# best hyperparameter values determined by Bayesian optimization
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = bestparams_bo_auc$eta,
                    scale_pos_weight = bestparams_bo_auc$scale_pos_weight,
                    max_depth = bestparams_bo_auc$max_depth,
                    min_child_weight = bestparams_bo_auc$min_child_weight,
                    subsample = bestparams_bo_auc$subsample,
                    colsample_bytree = bestparams_bo_auc$colsample_bytree,
                    colsample_bynode = bestparams_bo_auc$colsample_bynode,
                    lambda = bestparams_bo_auc$lambda,
                    alpha = bestparams_bo_auc$alpha,
                    gamma = bestparams_bo_auc$alpha,
                    max_delta_step = bestparams_bo_auc$max_delta_step)

# optimum number of nrounds chosen during final evaluation
best_nrounds <- round(mean(unlist(final_results$final_evaluation$best_nrounds)))

# train the model on the full dataset
dtrain_full <- xgboost::xgb.DMatrix(data = as.matrix(metagen[, all_feat_cols]),
                                    label = metagen$condition_numeric)

final_model <- xgboost::xgb.train(params = best_params,
                                  data = dtrain_full,
                                  nrounds = best_nrounds,
                                  verbose = 1)


##########################################################################################
########   XGBOOST LOGLOSS MODEL - SHAP VALUES - DEPENDENCE AND INTERACTION PLOTS  #######
##########################################################################################

# compute tree SHAP values
shap_values <- predict(final_model,
                       newdata = dtrain_full,
                       predcontrib = TRUE)

dim(shap_values) # one row per sample, one column per feature + bias term (initial logit)
shap_summary <- as.data.frame(shap_values)
shap_summary$BIAS <- NULL  # remove bias term

# mean absolute SHAP value per feature
shap_mean_abs <- sort(colMeans(abs(shap_summary)), decreasing = TRUE)
shap_mean_abs[1:20]

shap_mean_abs <- as.data.frame(shap_mean_abs) # covert to data.frame
shap_mean_abs$feature <- rownames(shap_mean_abs)


# SHAP plots with SHAPforxgboost
library(SHAPforxgboost)
shap_long <- shap.prep(xgb_model = final_model, X_train = as.matrix(metagen[, all_feat_cols]))
shap.plot.summary(shap_long)

# shap.plot.summary with top features (SHAP value > 0.005)
top_shap <- shap_mean_abs %>%
  group_by(feature) %>%
  filter(shap_mean_abs > 0.005)

shap_long_filtered <- shap_long %>%
  filter(variable %in% top_shap$feature) %>%
  mutate(variable = factor(variable, levels = unique(variable)))

shap.plot.summary(shap_long_filtered)


### SHAP dependence plots
# how SHAP values for a given feature vary as the input values for the feature vary (CLR transformed relative abundance)

# plot dependence plot for feature of interest
shap.plot.dependence(data_long = shap_long, x = "Lachnoclostridium_sp._YL32", y = NULL)
shap.plot.dependence(data_long = shap_long, x = "Petrimonas_mucosa", y = NULL)


# filter for the feature of interest
feature_name <- "Lachnoclostridium_sp._YL32"
shap_dep <- shap_long %>%
  filter(variable == feature_name)

# plot dependence plot 
ggplot(shap_dep, aes(x = rfvalue, y = value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = "Feature value (CLR-transformed relative abundances)",
       y = "SHAP value") +
  theme_minimal()


### SHAP interaction values
# how pairs of features interact in affecting the prediction

# compute SHAP interaction values (3D array)
interaction_values <- predict(final_model,
                              newdata = as.matrix(metagen[, all_feat_cols]),
                              predinteraction = TRUE)


### top feature pairs based on interaction strength (pairwise)
# aggregate across samples (interaction matrix)
mean_interactions <- apply(abs(interaction_values), c(2, 3), mean)

feature_names <- colnames(mean_interactions) # get feature names

# convert to long format
interaction_df <- as.data.frame(mean_interactions)
interaction_df$Feature1 <- feature_names
interaction_long <- pivot_longer(interaction_df, 
                                 cols = -Feature1, 
                                 names_to = "Feature2", 
                                 values_to = "InteractionStrength")

# remove self-interactions (diagonal)
interaction_long <- interaction_long %>%
  filter(Feature1 != Feature2)

# remove duplicate symmetric pairs
interaction_long <- interaction_long %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Feature1, Feature2)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  ungroup()

# top interactions
top_interactions <- interaction_long %>%
  arrange(desc(InteractionStrength)) %>%
  head(20)
top_interactions

# plot top interactions (heatmap)
ggplot(top_interactions, aes(x = Feature1, y = Feature2, fill = InteractionStrength)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top SHAP feature interactions", x = "Feature 1", y = "Feature 2")


### show strongest-interacting features (feature-wise sum)
# compute sum of interactions for each feature, excluding self-interactions (diagonals)
interaction_strength_per_feature <- rowSums(mean_interactions) - diag(mean_interactions)

# rank features by interaction strength
top_features <- sort(interaction_strength_per_feature, decreasing = TRUE)[1:20]
top_feature_names <- names(top_features)
top_interaction_matrix <- mean_interactions[top_feature_names, top_feature_names] # subset matrix

# plot strongest interactions (heatmap)
heatmap(top_interaction_matrix,
        main = "SHAP interaction heatmap (Top features only)",
        col = viridis::viridis(100), 
        margins = c(8, 8))


### dependence plots for top interacting pairs
shap_f1 <- shap_long %>% filter(variable == "Acutalibacter_muris") # get SHAP values for Acutalibacter_muris
shap_f1$interactor_value <- metagen[["Enterocloster_clostridioformis"]] # add abundance values for Enterocloster_clostridioformis

# plot SHAP values versus abundance for Acutalibacter_muris and color by abundance of Enterocloster_clostridioformis
ggplot(shap_f1, aes(x = value, y = rfvalue, color = interactor_value)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(name = "Enterocloster_clostridioformis") +
  labs(title = paste("SHAP dependence plot:", "Acutalibacter_muris", "vs", "Enterocloster_clostridioformis"),
       x = paste("Acutalibacter_muris", "(abundance)"),
       y = paste("SHAP value for", "Acutalibacter_muris")) +
  theme_minimal()


# loop through top 20 interactions
for (i in 1:nrow(top_interactions)) {
  f1 <- top_interactions$Feature1[i]
  f2 <- top_interactions$Feature2[i]
  
  # get SHAP values for f1
  shap_f1 <- shap_long %>% filter(variable == f1)
  
  # add the interacting feature's abundance from the metagen dataset
  shap_f1$interactor_value <- metagen[[f2]]
  
  # plot one feature's SHAP values versus its abundance and colored by the interacting feature's abundance
  p <- ggplot(shap_f1, aes(x = value, y = rfvalue, color = interactor_value)) +
    geom_point(alpha = 0.6) +
    scale_color_viridis_c(name = f2) +
    labs(title = paste("SHAP dependence plot:", f1, "vs", f2),
         x = paste(f1, "(abundance)"),
         y = paste("SHAP value for", f1)) +
    theme_minimal()
  
  print(p)
}


### partial dependence plots (PDP)
# shows the average effect of the feature on predicted outcome
library(pdp)

# prepare data matrix without label
X <- as.matrix(metagen[, all_feat_cols])

# partial dependence for Lachnoclostridium_sp._YL32
pdp_lachno <- partial(object = final_model,
                      pred.var = "Lachnoclostridium_sp._YL32",
                      train = metagen[, all_feat_cols],
                      grid.resolution = 20, # number of grid points across the feature's range
                      type = "regression", # output of the model is probabilities
                      prob = TRUE, # return probabilities (not logits)
                      plot = FALSE)

ggplot(pdp_lachno, aes(x = Lachnoclostridium_sp._YL32, y = yhat)) +
  geom_line(color = "blue") + theme_minimal() +
  labs(title = "Partial dependence plot",
       x = "Lachnoclostridium_sp._YL32 (Abundance)", y = "Predicted probability of disease")


# partial dependence for Petrimonas_mucosa
pdp_petrim <- partial(object = final_model,
                      pred.var = "Petrimonas_mucosa",
                      train = metagen[, all_feat_cols],
                      grid.resolution = 20, # number of grid points across the feature's range
                      type = "regression", # output of the model is probabilities
                      prob = TRUE, # return probabilities (not logits)
                      plot = FALSE)

ggplot(pdp_petrim, aes(x = Petrimonas_mucosa, y = yhat)) +
  geom_line(color = "blue") + theme_minimal() +
  labs(title = "Partial dependence plot",
       x = "Petrimonas_mucosa (Abundance)", y = "Predicted probability of disease")


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
#   [1] pdp_0.8.2                     SHAPforxgboost_0.1.3          DiceKriging_1.6.0            
# [4] future.apply_1.20.0           future_1.58.0                 doParallel_1.0.17            
# [7] iterators_1.0.14              foreach_1.5.2                 ParBayesianOptimization_1.2.6
# [10] Matrix_1.7-3                  MLmetrics_1.1.3               pROC_1.18.5                  
# [13] caret_7.0-1                   lattice_0.22-7                xgboost_1.7.11.1             
# [16] compositions_2.0-8            lubridate_1.9.4               forcats_1.0.0                
# [19] stringr_1.5.1                 dplyr_1.1.4                   purrr_1.0.4                  
# [22] readr_2.1.5                   tidyr_1.3.1                   tibble_3.3.0                 
# [25] tidyverse_2.0.0               ggplot2_3.5.2                
# 
# loaded via a namespace (and not attached):
#   [1] gridExtra_2.3        rlang_1.1.6          magrittr_2.0.3       e1071_1.7-16        
# [5] compiler_4.5.0       mgcv_1.9-3           vctrs_0.6.5          reshape2_1.4.4      
# [9] lhs_1.2.0            pkgconfig_2.0.3      crayon_1.5.3         backports_1.5.0     
# [13] labeling_0.4.3       utf8_1.2.6           prodlim_2025.04.28   tzdb_0.5.0          
# [17] jsonlite_2.0.0       recipes_1.3.1        tweenr_2.0.3         broom_1.0.8         
# [21] R6_2.6.1             stringi_1.8.7        RColorBrewer_1.1-3   parallelly_1.45.0   
# [25] car_3.1-3            rpart_4.1.24         Rcpp_1.0.14          splines_4.5.0       
# [29] nnet_7.3-20          timechange_0.3.0     tidyselect_1.2.1     viridis_0.6.5       
# [33] rstudioapi_0.17.1    dichromat_2.0-0.1    abind_1.4-8          timeDate_4041.110   
# [37] codetools_0.2-20     listenv_0.9.1        plyr_1.8.9           withr_3.0.2         
# [41] survival_3.8-3       bayesm_3.1-6         proxy_0.4-27         polyclip_1.10-7     
# [45] pillar_1.10.2        ggpubr_0.6.0         carData_3.0-5        tensorA_0.36.2.1    
# [49] checkmate_2.3.2      stats4_4.5.0         generics_0.1.4       dbscan_1.2.2        
# [53] hms_1.1.3            scales_1.4.0         globals_0.18.0       class_7.3-23        
# [57] glue_1.8.0           tools_4.5.0          robustbase_0.99-4-1  data.table_1.17.4   
# [61] ModelMetrics_1.2.2.2 gower_1.0.2          ggsignif_0.6.4       grid_4.5.0          
# [65] ipred_0.9-15         nlme_3.1-168         BBmisc_1.13          ggforce_0.4.2       
# [69] Formula_1.2-5        cli_3.6.5            viridisLite_0.4.2    lava_1.8.1          
# [73] gtable_0.3.6         DEoptimR_1.1-3-1     rstatix_0.7.2        digest_0.6.37       
# [77] farver_2.1.2         lifecycle_1.0.4      hardhat_1.4.1        MASS_7.3-65

