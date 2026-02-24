# XGBoost Workflow - Metagenomics Data

# load libraries
library(ggplot2)
library(ggrepel)
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
library(purrr)
library(SHAPforxgboost)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name
meta$condition <- as.factor(meta$condition) # convert condition column into a factor
table(meta$condition)
meta <- meta[, -3] # remove unnecessary meta column

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
feat_clr <- as.data.frame(feat_clr)

# rownames of metadata need to match the rownames of the feature table
all(rownames(meta) == rownames(feat_clr))
feat_clr$sample_name <- rownames(feat_clr) # add column sample_name


### merge metadata and feature table
metagen <- merge(meta, feat_clr, by = "sample_name", all.x = TRUE)
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


##########################################################################################
###   BASELINE XGBOOST MODEL + EARLY STOPPING - 5-FOLD CROSS-VALIDATION + 50 REPEATS   ###
##########################################################################################

# data to be used in the model
str(metagen)

# set seed
set.seed(1234)

# 50 repeats of stratified 5-fold cross-validation
folds <- createMultiFolds(metagen$condition, k = 5, times = 50)

# create list to store performance metrics
xgb_metrics <- list() # list to store performance metrics
xgb_importances <- list() # list to store feature importances
eval_logs <- list() # list to store evaluation logs for overfitting analysis
best_nrounds_list <- list() # list to store best nrounds

# loop through the folds
for (r in names(folds)) {
  
  # splits the dataset into training and testing sets
  train_idx <- folds[[r]] # train indices
  train_data <- metagen[train_idx, ] # training data
  test_data <- metagen[-train_idx, ] # testing data
  
  # convert to DMatrix
  dtrain <- xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), label = train_data$condition_numeric)
  dtest <- xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), label = test_data$condition_numeric)
  
  # set params
  params <- list(objective = "binary:logistic", eval_metric = "logloss")
  
  # cross-validation with early stopping
  xgb_cv <- xgb.cv(data = dtrain,
                   params = params,
                   nrounds = 5000,
                   nfold = 5,
                   early_stopping_rounds = 50,
                   maximize = FALSE,
                   stratified = TRUE,
                   showsd = FALSE,
                   verbose = 0)
  
  # best number of boosting rounds
  best_nrounds <- xgb_cv$early_stop$best_iteration
  
  # store evaluation_log from xgb.cv()
  eval_log <- xgb_cv$evaluation_log
  eval_log$fold_id <- r 
  eval_logs[[r]] <- eval_log
  
  # train XGBoost model with best number of boosting rounds
  xgb_model <- xgb.train(data = dtrain,
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
  
  
  ### store metrics with repeat index
  repeat_id = as.numeric(sub(".*Rep", "", r)) # extract repeat id
  xgb_metrics[[r]] <- list(cm = cm, auc = auc_val, logloss = logloss, 
                           repeat_id = repeat_id) # store confusion matrices
  imp <- xgb.importance(feature_names = all_feat_cols, model = xgb_model)
  xgb_importances[[r]] <- imp # store feature importances
  best_nrounds_list[[r]] <- best_nrounds
}


### calculate performance statistics
# create data frame to store fold-level results
results_df <- data.frame(repeat_id = numeric(),
                         balanced_accuracy = numeric(),
                         f1_score = numeric(),
                         precision = numeric(),
                         sensitivity = numeric(),
                         specificity = numeric(), 
                         auc_values = numeric(),
                         logloss_values = numeric())

# extract fold-level metrics
for (r in names(xgb_metrics)) {
  
  cm <- xgb_metrics[[r]]$cm
  auc_val <- xgb_metrics[[r]]$auc
  logloss_val <- xgb_metrics[[r]]$logloss
  repeat_id <- xgb_metrics[[r]]$repeat_id
  nrounds_values = best_nrounds_list[[r]]
  
  results_df <- rbind(results_df, data.frame(repeat_id = repeat_id,
                                             balanced_accuracy = cm$byClass["Balanced Accuracy"],
                                             f1_score = cm$byClass["F1"],
                                             precision = cm$byClass["Precision"],
                                             sensitivity = cm$byClass["Sensitivity"],
                                             specificity = cm$byClass["Specificity"],
                                             auc_values = auc_val,
                                             logloss_values = logloss_val,
                                             nrounds_values = nrounds_values))
}

# average metrics within each repeat
repeat_summary <- results_df %>%
  group_by(repeat_id) %>%
  summarise(balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
            f1_score = mean(f1_score, na.rm = TRUE),
            precision = mean(precision, na.rm =TRUE),
            sensitivity = mean(sensitivity, na.rm = TRUE),
            specificity = mean(specificity, na.rm = TRUE),
            auc_values = mean(auc_values, na.rm = TRUE),
            logloss_values = mean(logloss_values, na.rm = TRUE),
            nrounds_values = mean(nrounds_values, na.rm = TRUE))

# final summary across repeats
xgb_summary <- data.frame(mean_bal_acc = mean(repeat_summary$balanced_accuracy),
                          sd_bal_acc = sd(repeat_summary$balanced_accuracy),
                          mean_f1 = mean(repeat_summary$f1_score),
                          sd_f1 = sd(repeat_summary$f1_score),
                          mean_precision = mean(repeat_summary$precision),
                          sd_precision = sd(repeat_summary$precision),
                          mean_sens = mean(repeat_summary$sensitivity),
                          sd_sens = sd(repeat_summary$sensitivity),
                          mean_spec = mean(repeat_summary$specificity),
                          sd_spec = sd(repeat_summary$specificity),
                          mean_auc = mean(repeat_summary$auc_values),
                          sd_auc = sd(repeat_summary$auc_values),
                          mean_logloss = mean(repeat_summary$logloss_values),
                          sd_logloss = sd(repeat_summary$logloss_values),
                          mean_nrounds = mean(repeat_summary$nrounds_values, na.rm = TRUE),
                          sd_nrounds = sd(repeat_summary$nrounds_values, na.rm = TRUE))

xgb_summary


### feature importances metrics (gain, cover, frequency)
all_importances <- bind_rows(xgb_importances)
mean_importance <- all_importances %>%
  group_by(Feature = Feature) %>%
  summarise(mean_gain = mean(Gain, na.rm = TRUE),
            mean_cover = mean(Cover, na.rm = TRUE),
            mean_freq = mean(Frequency, na.rm = TRUE),
            freq_selected = n()) %>%
  arrange(desc(freq_selected))
head(mean_importance, 20)

# plot feature importance by stability
top_mean_importance = mean_importance %>% 
  slice(1:20)

# feature stability by mean gain
ggplot(top_mean_importance, aes(x = freq_selected, y = mean_gain, label = Feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + 
  theme_minimal() + labs(title = "Feature importance versus stability", 
                         x = "Frequency selected", y = "Mean gain")

# feature stability by mean cover
ggplot(top_mean_importance, aes(x = freq_selected, y = mean_cover, label = Feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + 
  theme_minimal() + labs(title = "Feature importance versus stability", 
                         x = "Frequency selected", y = "Mean cover")


### overfitting analysis - logloss over boosting rounds
logloss_all <- bind_rows(eval_logs, .id = "Fold")

# compute the mean test logloss per iteration
mean_logloss <- logloss_all %>%
  group_by(iter) %>%
  summarise(mean_test_logloss = mean(test_logloss_mean, na.rm = TRUE), .groups = "drop")
logloss_all <- logloss_all %>% left_join(mean_logloss, by = "iter")

# plot test logloss curve
logloss_all_filt <- logloss_all %>% filter(iter <= 100) # number of boosting rounds shown
ggplot(logloss_all_filt, aes(x = iter, y = test_logloss_mean, group = Fold)) +
  geom_line(alpha = 0.2, color = "black") +
  geom_line(aes(y = mean_test_logloss), color = "blue", linewidth = 1) + # mean curve
  labs(title = "Test logloss over boosting rounds",
       x = "Boosting round", y = "Test logloss")

# plot train logloss curve
logloss_all_filt <- logloss_all %>% filter(iter <= 100) # number of boosting rounds shown
ggplot(logloss_all_filt, aes(x = iter, y = train_logloss_mean, group = Fold)) +
  geom_line(alpha = 0.2, color = "red") +
  labs(title = "Train logloss over boosting rounds",
       x = "Boosting round", y = "Train logloss")

# pivot to long format for logloss
logloss_long <- logloss_all %>%
  select(iter, Fold, train_logloss_mean, test_logloss_mean) %>%
  pivot_longer(cols = c(train_logloss_mean, test_logloss_mean), names_to = "set", values_to = "logloss") %>%
  mutate(set = ifelse(set == "train_logloss_mean", "Train", "Test"))

# plot test and train logloss curve
logloss_long_filt <- logloss_long %>% filter(iter <= 100) # number of boosting rounds shown
ggplot(logloss_long_filt, aes(x = iter, y = logloss, group = interaction(Fold, set), color = set)) +
  geom_line(alpha = 0.2) +
  scale_color_manual(values = c("Train" = "red", "Test" = "black")) +
  labs(title = "Train and test logloss over boosting rounds",
       x = "Boosting round", y = "Logloss", color = "Data")

# quantify generalization error over time per eta
logloss_gap <- logloss_long %>%
  pivot_wider(names_from = set, values_from = logloss) %>%
  mutate(gap = Test - Train) %>%
  group_by(iter) %>%
  summarise(mean_gap = mean(gap, na.rm = TRUE), .groups = "drop")

# plot logloss gap between test and train curves
logloss_gap_filt <- logloss_gap %>% filter(iter <= 100) # number of boosting rounds shown
ggplot(logloss_gap_filt, aes(x = iter, y = mean_gap)) +
  geom_line() +
  labs(title = "Mean logloss gap (test logloss - train logloss) across boosting rounds",
       y = "Gap", x = "Boosting round")


### best_nrounds is very low -> model overfitting very quickly (default hyperparameters too aggressive)


################################################################################################################
###   XGBOOST MODEL (EARLY STOPPING + 5-FOLD CV + 50 REPEATS) - FUNCTION TO GRID TUNE HYPERPARAMETERS - ETA  ###
################################################################################################################

### function to grid tune hyperparameters with xgboost
### early stopping + 5 fold cross-validation + 50 repeats
tune_xgb_param <- function(param_grid_name,
                           param_grid_values,
                           base_params = NULL,
                           dataframe,
                           all_feat_cols,
                           target_var,
                           target_var_numeric,
                           n_repeats = 50) {
  
  # default hyperparameter values
  default_params <- list(objective = "binary:logistic",
                         eval_metric = "logloss",
                         eta = 0.005, # actual default value is horrible for model
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
  
  
  # if no base_params provided, use values in default_params
  if (is.null(base_params)) {
    base_params <- default_params
  } else {
    # if hyperparameter values supplied, use those values
    base_params <- modifyList(default_params, base_params)
  }

  # remove the parameter being tuned
  base_params <- base_params[setdiff(names(base_params), param_grid_name)]
  
  # set seed
  set.seed(1234)
  
  # n_repeats (default 50) repeats of stratified 5-fold cross-validation
  folds <- caret::createMultiFolds(dataframe[[target_var]], k = 5, times = n_repeats)
  
  # list to store all results (metrics, feature importance, best nrounds, logloss) for each param grid value
  all_results <- list()
  
  # loop over val in param_grid_values
  for (val in param_grid_values) {
    cat("Running with", param_grid_name, "=", val, "\n")
    
    # update parameters with current val
    params <- base_params
    params[[param_grid_name]] <- val
    
    # create lists to store results 
    xgb_metrics <- list() # list to store performance metrics
    xgb_importances <- list() # list to store feature importances
    best_nrounds_list <- list() # list to store best nrounds
    eval_logs <- list() # list to store evaluation logs for overfitting analysis
    
    # loop over each of the cross-validation folds
    for (r in names(folds)) {
      
      # splits the dataset into training and testing sets
      train_idx <- folds[[r]] # train indices 
      train_data <- dataframe[train_idx, ] # training data
      test_data  <- dataframe[-train_idx, ] # testing data
      
      # convert to DMatrix
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), label = train_data[[target_var_numeric]])
      dtest <- xgboost::xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), label = test_data[[target_var_numeric]])
      
      # run internal 5-fold cross-validation to determine best nrounds
      xgb_cv <- xgboost::xgb.cv(data = dtrain,
                                params = params,
                                nrounds = 5000,
                                nfold = 5,
                                early_stopping_rounds = 50,
                                maximize = FALSE,
                                stratified = TRUE,
                                showsd = FALSE,
                                verbose = 0)
      
      # best number of boosting rounds
      best_nrounds <- xgb_cv$early_stop$best_iteration
      
      # store evaluation_log from xgb.cv()
      eval_log <- xgb_cv$evaluation_log
      eval_log$fold_id <- r
      eval_logs[[r]] <- eval_log
      
      # train final xgboost model with best nrounds
      xgb_model <- xgboost::xgb.train(data = dtrain,
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
      
      # calculate logloss
      eps <- 1e-15 # avoid log(0)
      prob_clipped <- pmin(pmax(preds_prob, eps), 1 - eps)
      logloss <- -mean(test_data[[target_var_numeric]] * log(prob_clipped) +
                         (1 - test_data[[target_var_numeric]]) * log(1 - prob_clipped))
      
      
      ### store performance metrics, importances, best nrounds and logloss
      repeat_id = as.numeric(sub(".*Rep", "", r)) # extract repeat id
      xgb_metrics[[r]] <- list(cm = cm, auc = auc_val, logloss = logloss,
                               repeat_id = repeat_id)
      imp <- xgboost::xgb.importance(feature_names = all_feat_cols, model = xgb_model)
      xgb_importances[[r]] <- imp
      best_nrounds_list[[r]] <- best_nrounds
    }
    
    
    ### calculate performance statistics
    # create data frame to store fold-level results
    results_df <- data.frame(repeat_id = numeric(),
                             balanced_accuracy = numeric(),
                             f1_score = numeric(),
                             precision = numeric(),
                             sensitivity = numeric(),
                             specificity = numeric(), 
                             auc_values = numeric(),
                             logloss_values = numeric())
    
    
    # loop through stored lists to extract metrics
    for (r in names(xgb_metrics)) {
      cm <- xgb_metrics[[r]]$cm
      auc_val <- xgb_metrics[[r]]$auc
      logloss_val <- xgb_metrics[[r]]$logloss
      repeat_id <- xgb_metrics[[r]]$repeat_id
      nrounds_values = best_nrounds_list[[r]]
      
      results_df <- rbind(results_df, data.frame(repeat_id = repeat_id,
                                                 balanced_accuracy = cm$byClass["Balanced Accuracy"],
                                                 f1_score = cm$byClass["F1"],
                                                 precision = cm$byClass["Precision"],
                                                 sensitivity = cm$byClass["Sensitivity"],
                                                 specificity = cm$byClass["Specificity"],
                                                 auc_values = auc_val,
                                                 logloss_values = logloss_val,
                                                 nrounds_values = nrounds_values))
    }
    
    # average metrics within each repeat
    repeat_summary <- results_df %>%
      group_by(repeat_id) %>%
      summarise(balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
                f1_score = mean(f1_score, na.rm = TRUE),
                precision = mean(precision, na.rm =TRUE),
                sensitivity = mean(sensitivity, na.rm = TRUE),
                specificity = mean(specificity, na.rm = TRUE),
                auc_values = mean(auc_values, na.rm = TRUE),
                logloss_values = mean(logloss_values, na.rm = TRUE),
                nrounds_values = mean(nrounds_values, na.rm = TRUE))
    
    # final summary across repeats
    summary_df <- data.frame(param = val,
                             mean_bal_acc = mean(repeat_summary$balanced_accuracy, na.rm = TRUE),
                             sd_bal_acc = sd(repeat_summary$balanced_accuracy, na.rm = TRUE),
                             mean_f1 = mean(repeat_summary$f1_score, na.rm = TRUE),
                             sd_f1 = sd(repeat_summary$f1_score, na.rm = TRUE),
                             mean_precision = mean(repeat_summary$precision, na.rm = TRUE),
                             sd_precision = sd(repeat_summary$precision, na.rm = TRUE),
                             mean_sens = mean(repeat_summary$sensitivity, na.rm = TRUE),
                             sd_sens = sd(repeat_summary$sensitivity, na.rm = TRUE),
                             mean_spec = mean(repeat_summary$specificity, na.rm = TRUE),
                             sd_spec = sd(repeat_summary$specificity, na.rm = TRUE),
                             mean_auc = mean(repeat_summary$auc_values, na.rm = TRUE),
                             sd_auc = sd(repeat_summary$auc_values, na.rm = TRUE),
                             mean_logloss = mean(repeat_summary$logloss_values, na.rm = TRUE),
                             sd_logloss = sd(repeat_summary$logloss_values, na.rm = TRUE),
                             mean_nrounds = mean(repeat_summary$nrounds_values, na.rm = TRUE),
                             sd_nrounds = sd(repeat_summary$nrounds_values, na.rm = TRUE))
    

    # aggregate feature importances
    all_importances <- dplyr::bind_rows(xgb_importances)
    mean_importance <- all_importances %>%
      dplyr::group_by(Feature = Feature) %>%
      dplyr::summarise(mean_gain = mean(Gain, na.rm = TRUE),
                       mean_cover = mean(Cover, na.rm = TRUE),
                       mean_freq = mean(Frequency, na.rm = TRUE),
                       freq_selected = dplyr::n()) %>%
      dplyr::arrange(desc(freq_selected))
    
    # store all values under current param value name
    all_results[[paste0(param_grid_name, "_", val)]] <- list(summary = summary_df,
                                                             feature_importance = mean_importance,
                                                             raw_metrics = xgb_metrics,
                                                             best_nrounds = best_nrounds_list,
                                                             eval_logs = eval_logs)
  }
  
  return(all_results)
}


### arguments to include
current_param <- "eta" # set current parameter
param_grid_values <- c(0.005, 0.01, 0.05, 0.1, 0.3) # grid of param values to test


### run model evaluation with eta parameter grid search
xgboost_model_eta <- tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                                    param_grid_values = param_grid_values, # grid of param values to test
                                    dataframe = metagen, # data to use
                                    all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                                    target_var = "condition",
                                    target_var_numeric = "condition_numeric")


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
tune_summarize_performance <- function(results_list) {
  dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
}
tune_summarize_performance(xgboost_model_eta)

 
### feature importance (importance of features across parameter settings)
tune_feature_importance <- function(results_list, param_name = current_param) {
  purrr::imap_dfr(results_list, function(res, name) {
    df <- res$feature_importance
    df[[param_name]] <- as.numeric(gsub(paste0(param_name, "_"), "", name))
    df
  })
}
feature_freq_eta <- tune_feature_importance(xgboost_model_eta, param_name = current_param)
feature_freq_eta


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency <- function(all_importances_df, feature, param_name = current_param) {
  filtered_df <- dplyr::filter(all_importances_df, Feature == feature)
  ggplot(filtered_df, aes(x = .data[[param_name]], y = freq_selected)) +
    geom_line() + geom_point() +
    labs(title = paste("Frequency of selection vs", param_name, "for", feature),
      x = param_name, y = "Frequency of selection")
}
tune_plot_feature_selection_frequency(feature_freq_eta, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability <- function(all_importances_df, x = "freq_selected", y = "mean_gain", color_by = current_param) {
  ggplot(all_importances_df, aes(x = .data[[x]], y = .data[[y]], color = as.factor(.data[[color_by]]))) +
    geom_point(alpha = 0.6) + theme_minimal() +
    labs(title = "Feature importance stability", x = x, y = y, color = color_by)
}
tune_plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_eta, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single <- function(all_importances_df, x = "freq_selected", y = "mean_gain",
                                               param_name, param_value, top_n = 20) {
  df_filtered <- all_importances_df %>%
    dplyr::filter(.data[[param_name]] == param_value) %>%
    dplyr::slice(1:top_n) # all_importances_df is arranged by freq_selected
  
  ggplot(df_filtered, aes(x = .data[[x]], y = .data[[y]], label = Feature)) +
    geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) +
    theme_minimal() + labs(title = paste0("Feature importance versus stability for ", param_name, " = ", param_value), 
                           x = x, y = y)
}
tune_plot_feature_stability_single(feature_freq_eta, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 0.005)
  

### logloss and overfitting analysis
tune_extract_logloss_df <- function(results_list, param_name = current_param) {
  purrr::imap_dfr(results_list, function(res, name) {
    logs <- res$eval_logs
    df <- dplyr::bind_rows(logs, .id = "Fold")
    df[[param_name]] <- name
    df
  })
}
logloss_eta <- tune_extract_logloss_df(xgboost_model_eta, param_name = current_param)


### plot logloss over boosting rounds (type = "train", "test", or "both")
tune_plot_logloss_curve <- function(logloss_df, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = NULL) {
  df <- logloss_df %>%
    mutate(!!param_name := as.numeric(gsub(paste0(param_name, "_"), "", .data[[param_name]])))
  
  # limit number of boosting rounds displayed
  if (!is.null(max_rounds)) {
    df <- df %>% dplyr::filter(iter <= max_rounds)
  }
  
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
    
    # limit df_long if max_rounds is set (for type = "both")
    if (!is.null(max_rounds)) {
      df_long <- df_long %>% dplyr::filter(iter <= max_rounds)
    }
    
    p <- ggplot(df_long, aes(x = iter, y = logloss, color = set, group = interaction(Fold, set))) +
      geom_line(alpha = 0.2) +
      scale_color_manual(values = c("train_logloss_mean" = "red", "test_logloss_mean" = "black")) +
      facet_wrap(as.formula(paste("~", param_name)), scales = "free_x")
  }
  
  p + labs(title = "Logloss curves", x = "Boosting round", y = "Logloss")
}
tune_plot_logloss_curve(logloss_eta, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 500)
tune_plot_logloss_curve(logloss_eta, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 500)
tune_plot_logloss_curve(logloss_eta, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
tune_logloss_gap <- function(logloss_df, param_name = current_param) {
  logloss_df %>%
    mutate(!!param_name := as.numeric(gsub(paste0(param_name, "_"), "", .data[[param_name]]))) %>%
    group_by_at(c(param_name, "iter")) %>%
    summarise(mean_train = mean(train_logloss_mean, na.rm = TRUE),
              mean_test = mean(test_logloss_mean, na.rm = TRUE),
              gap = mean(test_logloss_mean - train_logloss_mean, na.rm = TRUE),
              .groups = "drop")
}
logloss_gap_eta <- tune_logloss_gap(logloss_eta, param_name = current_param)


### plot logloss gap
tune_plot_logloss_gap <- function(logloss_gap_df, param_name = current_param, max_rounds = NULL) {
  
  # limit number of boosting rounds displayed
  if (!is.null(max_rounds)) {
    logloss_gap_df <- logloss_gap_df %>% dplyr::filter(iter <= max_rounds)
  }
  
  ggplot(logloss_gap_df, aes(x = iter, y = gap)) +
    geom_line() + theme_minimal() +
    facet_wrap(as.formula(paste("~", param_name)), scales = "free_x") +
    labs(title = paste("Mean logloss gap (test - train) across boosting rounds"),
      x = "Boosting round", y = "Gap")
}
tune_plot_logloss_gap(logloss_gap_eta, param_name = current_param, max_rounds = 500)


###################################################################################################
###   XGBOOST MODEL - HYPERPARAMETER TUNING + EARLY STOPPING + SHAP SUBSET - SCALE POS WEIGHT   ###
###################################################################################################

### arguments to include
current_param <- "scale_pos_weight" # set current parameter
param_grid_values <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 1.0)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


############################################################################################
###   XGBOOST MODEL - HYPERPARAMETER TUNING + EARLY STOPPING + SHAP SUBSET - MAX DEPTH   ###
############################################################################################

### arguments to include
current_param <- "max_depth" # set current parameter
param_grid_values <- c(2, 3, 4, 5, 6, 7) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 5)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


#####################################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - MIN CHILD WEIGHT   ###
#####################################################################################################

### arguments to include
current_param <- "min_child_weight" # set current parameter
param_grid_values <- c(1, 2, 3, 4, 5, 6) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 4)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


##############################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - SUBSAMPLE   ###
##############################################################################################

### arguments to include
current_param <- "subsample" # set current parameter
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 0.7)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


#####################################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - COLSAMPLE_BYTREE   ###
#####################################################################################################

### arguments to include
current_param <- "colsample_bytree" # set current parameter
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 0.8)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


#####################################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - COLSAMPLE_BYNODE   ###
#####################################################################################################

### arguments to include
current_param <- "colsample_bynode" # set current parameter
param_grid_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 0.6)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


###########################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - LAMBDA   ###
###########################################################################################

### arguments to include
current_param <- "lambda" # set current parameter
param_grid_values <- c(0.1, 0.5, 1, 1.5, 5, 10) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 1)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


##########################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - ALPHA   ###
##########################################################################################

### arguments to include
current_param <- "alpha" # set current parameter
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 0.5)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


##########################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING + SHAP SUBSET - GAMMA   ###
##########################################################################################

### arguments to include
current_param <- "gamma" # set current parameter
param_grid_values <- c(0, 0.1, 0.5, 1, 1.5, 5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = current_param, # name of parameter to tune
                      param_grid_values = param_grid_values,
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 1)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


#####################################################################################
###   XGBOOST MODEL - DEFAULT HYPERPARAMETERS + EARLY STOPPING - MAX_DELTA_STEP   ###
#####################################################################################

### arguments to include
current_param <- "max_delta_step" # set current parameter
param_grid_values <- c(0, 1, 2, 3, 4, 5) # grid of param values to test


### run model evaluation with current parameter grid search
assign(paste0("xgboost_model_", current_param),
       tune_xgb_param(param_grid_name = "current_param", # name of parameter to tune
                      param_grid_values = param_grid_values, # grid of param values to test
                      dataframe = metagen, # data to use
                      all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                      target_var = "condition",
                      target_var_numeric = "condition_numeric"))


### summarize performance metrics (mean and SD of balanced accuracy, f1, precision, sensitivity, specificity, auc, logloss, best nrounds)
shap_subset_object <- get(paste0("xgboost_model_", current_param)) # get model
tune_summarize_performance(shap_subset_object)


### feature importance (importance of features across parameter settings)
assign(paste0("feature_freq_", current_param),
       tune_feature_importance(shap_subset_object, param_name = current_param))
feature_freq_object <- get(paste0("feature_freq_", current_param))
feature_freq_object


### plot feature selection by frequency/selection across parameter values
tune_plot_feature_selection_frequency(feature_freq_object, feature = "Acutalibacter_muris", param_name = current_param)


### plot feature frequency/selection by importance (mean gain or mean cover)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_gain", color_by = current_param)
tune_plot_feature_stability(feature_freq_object, x = "freq_selected", y = "mean_cover", color_by = current_param)


### plot feature frequency/selection by importance for a single param value with feature names (mean gain or mean cover)
tune_plot_feature_stability_single(feature_freq_object, x = "freq_selected", y = "mean_gain", param_name = current_param, param_value = 3)


### logloss and overfitting analysis
assign(paste0("logloss_", current_param),
       tune_extract_logloss_df(shap_subset_object, param_name = current_param))
logloss_object <- get(paste0("logloss_", current_param))


### plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
tune_plot_logloss_curve(logloss_object, type = "test", show_mean = TRUE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "train", show_mean = FALSE, param_name = current_param, max_rounds = 1500)
tune_plot_logloss_curve(logloss_object, type = "both", show_mean = FALSE, param_name = current_param, max_rounds = 1500)


### calculate logloss gap (test loss - train loss) over boosting rounds and parameter values
assign(paste0("logloss_gap_", current_param),
       tune_logloss_gap(logloss_object, param_name = current_param))
logloss_gap_object <- get(paste0("logloss_gap_", current_param))

### plot logloss gap
tune_plot_logloss_gap(logloss_gap_object, param_name = current_param, max_rounds = 1500)


#################################################################################################
###   XGBOOST MODEL - NESTED CROSS VALIDATION WITH BAYESIAN OPTIMIZATION OF HYPERPARAMETERS   ###
#################################################################################################

# function to perform repeated nested cv with bayesian hyperparameter optimization
nested_cv_bayes <- function(dataframe,
                            all_feat_cols,
                            target_var = "condition", 
                            target_var_numeric = "condition_numeric",
                            bounds = bounds,
                            n_repeats = 10,
                            n_folds = 3,
                            bayes_n_repeats = 10,
                            bayes_n_folds = 3,
                            bayes_initPoints = 22,
                            bayes_iters.n = 10) {
  
  
  # create lists to store results
  xgb_metrics <- list() # list to store performance metrics
  xgb_importances <- list() # list to store feature importances
  outer_best_nrounds <- list() # list to store best nrounds
  outer_opt_params <- list() # list to store optimal hyperparameter values
  eval_logs <- list() # list to store evaluation logs for overfitting analysis
  
  
  ### outer loop
  # n_repeats of n_folds cross-validation
  for (r in 1:n_repeats) {
    
    # different seed for each repeat
    set.seed(1234 + r*10000)
    folds <- caret::createFolds(dataframe[[target_var]], k = n_folds, list = TRUE, returnTrain = FALSE)
    
    # loop over each cross-validation fold
    for (f in 1:n_folds) {
      cat("Fold", f, "of repeat", r, "\n")
      
      # split data into training and testing folds
      test_idx <- folds[[f]] # test indices
      test_data <- dataframe[test_idx, ] # training data
      train_data <- dataframe[-test_idx, ] # testing data 
      
      
      ### inner loop - Bayesian hyperparameter optimization
      
      # scoring function
      scoring_function <- function(eta, scale_pos_weight, max_depth, min_child_weight, 
                                   subsample, colsample_bytree, colsample_bynode, 
                                   lambda, alpha, gamma, max_delta_step) {
        
        # Make hyperparameter-specific seed
        seed_val <- as.integer(abs(sum(c(eta, scale_pos_weight, max_depth))) * 1e6)
        set.seed(seed_val)
        
        # convert to DMatrix
        dtrain <- xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), 
                              label = train_data[[target_var_numeric]])
        
        params <- list(objective = "binary:logistic",
                       eval_metric = "logloss",
                       nthread = 1, # to avoid oversubscription with BayesOpt
                       eta = eta,
                       scale_pos_weight = scale_pos_weight,
                       max_depth = as.integer(max_depth),
                       min_child_weight = as.integer(min_child_weight),
                       subsample = subsample,
                       colsample_bytree = colsample_bytree,
                       colsample_bynode = colsample_bynode,
                       lambda = lambda,
                       alpha = alpha,
                       gamma = gamma,
                       max_delta_step = as.integer(max_delta_step))
        
        # create numeric vector for logloss
        repeat_logloss <- numeric(bayes_n_repeats)
        
        # loop over each repeat
        for (br in 1:bayes_n_repeats) {
          
          # create bayes_n_folds-folds for cross-validation
          bayes_folds <- caret::createFolds(train_data[[target_var]], k = bayes_n_folds, list = TRUE, returnTrain = FALSE)
          
          xgb_cv <- tryCatch({
            xgboost::xgb.cv(params = params,
                            data = dtrain,
                            nrounds = 5000,
                            folds = bayes_folds,
                            stratified = TRUE,
                            showsd = FALSE,
                            early_stopping_rounds = 50,
                            verbose = 0)
          }, error = function(e) return(NULL))
          
          if (is.null(xgb_cv)) {
            repeat_logloss[br] <- Inf
          } else {
            
            # take the best iteration of logloss for the repeat
            repeat_logloss[br] <- xgb_cv$evaluation_log$test_logloss_mean[xgb_cv$early_stop$best_iteration]
          }
        }
        
        # average across repeats
        mean_logloss <- mean(repeat_logloss)
        
        # negative because Bayesian optimization maximizes the Score
        return(list(Score = -mean_logloss))
      }
      
      # register backend
      doParallel::registerDoParallel(parallel::detectCores() - 1)
      
      optObj <- bayesOpt(FUN = scoring_function,
                         bounds = bounds,
                         initPoints = bayes_initPoints,
                         iters.n = bayes_iters.n,
                         acq = "ei",
                         parallel = TRUE,
                         verbose = 1)
      
      # unregister the backend
      registerDoSEQ()
      
      # get optimal hyperparameter values 
      opt_param_values <- getBestPars(optObj)
      
      ### evaluation of model performance with Bayesian optimized hyperparameters
      # list of optimal hyperparameters determined by Bayesian optimization 
      opt_params <- c(list(objective = "binary:logistic", eval_metric = "logloss"), opt_param_values)
      
      # convert to DMAtrix
      dtrain_outer <- xgboost::xgb.DMatrix(data = as.matrix(train_data[, all_feat_cols]), 
                                           label = train_data[[target_var_numeric]])
      dtest_outer <- xgboost::xgb.DMatrix(data = as.matrix(test_data[, all_feat_cols]), 
                                          label = test_data[[target_var_numeric]])
      
      # run 5-fold cross-validation to determine best nrounds
      xgb_cv_outer <- xgboost::xgb.cv(data = dtrain_outer,
                                      params = opt_params,
                                      nrounds = 5000,
                                      nfold = 5,
                                      early_stopping_rounds = 50,
                                      maximize = FALSE,
                                      stratified = TRUE,
                                      showsd = FALSE,
                                      verbose = 0)
      
      # best number of boosting rounds
      best_nrounds_outer <- xgb_cv_outer$early_stop$best_iteration
      
      # store evaluation log
      eval_log <- xgb_cv_outer$evaluation_log
      eval_logs[[paste0("R", r, "_F", f)]] <- eval_log
      
      # train final model using best nrounds and optimal hyperparameters
      xgb_final_outer <- xgboost::xgb.train(data = dtrain_outer,
                                            params = opt_params,
                                            nrounds = best_nrounds_outer,
                                            verbose = 0)
      
      # evaluate on test set
      preds_prob <- predict(xgb_final_outer, dtest_outer)
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
      eps <- 1e-15
      prob_clipped <- pmin(pmax(preds_prob, eps), 1 - eps)
      logloss <- -mean(test_data[[target_var_numeric]] * log(prob_clipped) +
                         (1 - test_data[[target_var_numeric]]) * log(1 - prob_clipped))
      
      # store performance metrics (auc, cm and logloss)
      xgb_metrics[[paste0("R", r, "_F", f)]] <- list(cm = cm, auc = auc_val, logloss = logloss,
                                                     opt_params = opt_param_values,
                                                     best_nrounds = best_nrounds_outer)
      
      # feature importances 
      imp <- xgboost::xgb.importance(feature_names = all_feat_cols, model = xgb_final_outer)
      imp$Repeat_Fold <- paste0("R", r, "_F", f)
      xgb_importances[[paste0("R", r, "_F", f)]] <- imp
      
    }
  }
  
  ### summarize performance metrics
  per_fold_results <- list()
  
  # loop through stored lists to extract metrics
  for (key in names(xgb_metrics)) {
    cm <- xgb_metrics[[key]]$cm
    auc_val <- as.numeric(xgb_metrics[[key]]$auc)
    logloss_val <- xgb_metrics[[key]]$logloss
    best_nrounds <- xgb_metrics[[key]]$best_nrounds
    opt_params <- xgb_metrics[[key]]$opt_params
    
    # combine metrics + params into a single row
    per_fold_results[[key]] <- data.frame(repeat_fold = key,
                                          balanced_accuracy = cm$byClass["Balanced Accuracy"],
                                          f1_score = cm$byClass["F1"],
                                          precision = cm$byClass["Precision"],
                                          sensitivity = cm$byClass["Sensitivity"],
                                          specificity = cm$byClass["Specificity"],
                                          auc_values = auc_val,
                                          logloss_values = logloss_val,
                                          nrounds_values = best_nrounds,
                                          as.list(opt_params),
                                          check.names = FALSE,
                                          stringsAsFactors = FALSE,
                                          row.names = NULL)
  }
  
  # bind all folds together
  per_fold_df <- dplyr::bind_rows(per_fold_results) %>%
    arrange(desc(auc_values))
  
  # aggregate stored metrics into a data.frame
  summary_df <- per_fold_df %>%
    dplyr::summarise(mean_bal_acc = mean(balanced_accuracy, na.rm = TRUE),
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
    dplyr::group_by(feature = Feature) %>%
    dplyr::summarise(mean_gain = mean(Gain, na.rm = TRUE),
                     mean_cover = mean(Cover, na.rm = TRUE),
                     mean_freq = mean(Frequency, na.rm = TRUE),
                     freq_selected = dplyr::n()) %>%
    dplyr::arrange(desc(freq_selected))
  
  # return all metrics
  return(list(summary = summary_df,
              per_fold = per_fold_df,
              feature_importance = mean_importance,
              eval_logs = eval_logs))
  
}

# set bounds for bayesian hyperparameter optimization
bounds <- list(eta = c(0.001, 0.01),
               scale_pos_weight = c(0.5, 1.5),
               max_depth = c(2L, 7L),
               min_child_weight = c(1L, 3L),
               subsample = c(0.5, 1.0),
               colsample_bytree = c(0.5, 1.0),
               colsample_bynode = c(0.5, 1.0),
               lambda = c(0.1, 5),
               alpha = c(0, 1.5),
               gamma = c(0, 1.0),
               max_delta_step = c(0L, 5L))

bayes_obj <- nested_cv_bayes(dataframe = metagen,
                             all_feat_cols = setdiff(colnames(metagen), c("condition", "condition_numeric")),
                             target_var = "condition", 
                             target_var_numeric = "condition_numeric",
                             bounds = bounds,
                             n_repeats = 10,
                             n_folds = 3,
                             bayes_n_repeats = 10,
                             bayes_n_folds = 3,
                             bayes_initPoints = 22,
                             bayes_iters.n = 10)


### overall performance of model 
bayes_obj$summary

### per fold performance 
bayes_obj$per_fold

### feature importance (mean_gain, mean_cover, mean_freq)
feat_imp <- bayes_obj$feature_importance
feat_imp

### plot feature importance
tune_plot_feature_stability <- function(all_importances_df, x = "freq_selected", y = "mean_gain") {
  ggplot(all_importances_df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.6, color = "steelblue") +
    theme_minimal() +
    labs(title = "Feature importance stability", x = x, y = y)
}

tune_plot_feature_stability(bayes_obj$feature_importance, x = "freq_selected", y = "mean_gain")
tune_plot_feature_stability(bayes_obj$feature_importance, x = "freq_selected", y = "mean_cover")


### plot split frequency versus feature importance
# top 20 features by split frequency
top_features <- bayes_obj$feature_importance%>% 
  slice(1:15) # top 20 features by split frequency
top_features$proportion_selected <- top_features$freq_selected / 30

# plot mean cover versus proportion of models
ggplot(top_features, aes(x = proportion_selected, y = mean_gain, label = feature)) +
  geom_point(color = "steelblue", size = 3) + geom_text_repel(size = 3, max.overlaps = 8) + theme_minimal() +
  labs(x = "Proportion of models", y = "Mean gain",
       title = "Proportion of models versus mean gain")


### plot feature importance
plot_feature_importance <- function(feat_imp, 
                                    importance = c("gain", "cover", "frequency"), 
                                    tot_num_folds = 30, 
                                    importance_threshold = 0.05,
                                    freq_selected_threshold = 0.5) {
  
  # match argument
  importance <- match.arg(importance)
  
  # create data.frame
  all_feat_imp <- feat_imp %>%
    mutate(freq_selected = freq_selected / tot_num_folds) %>%
    mutate(in_atleast_X = ifelse(freq_selected >= freq_selected_threshold, 
                                 paste0("in_atleast_", freq_selected_threshold * 100, "%"), "other"))
  
  # choose column based on importance type
  y_col <- switch(importance,
                  gain = "mean_gain",
                  cover = "mean_cover",
                  frequency = "mean_freq")
  
  y_label <- switch(importance,
                    gain = "Mean gain",
                    cover = "Mean cover",
                    frequency = "Mean frequency")
  
  # filter by importance threshold
  plot_data <- all_feat_imp %>%
    filter(.data[[y_col]] > importance_threshold)
  
  # dynamic color mapping using setNames
  fill_colors <- setNames(c("indianred3", "steelblue"), 
                          c(paste0("in_atleast_", freq_selected_threshold * 100, "%"), "other"))
  
  # plot importance
  ggplot(plot_data,
         aes(x = reorder(.data$feature, .data[[y_col]]), y = .data[[y_col]], fill = .data$in_atleast_X)) +
    scale_fill_manual(values = fill_colors) +
    geom_col() + coord_flip() + theme_minimal(base_size = 12) +
    labs(title = y_label, x = "Feature", y = y_label, fill = "Feature selection frequency")
}

plot_feature_importance(feat_imp, importance = "gain", tot_num_folds = 50, importance_threshold = 0.02, freq_selected_threshold = 0.8)
plot_feature_importance(feat_imp, importance = "cover", tot_num_folds = 50, importance_threshold = 0.02, freq_selected_threshold = 0.8)
plot_feature_importance(feat_imp, importance = "frequency", tot_num_folds = 50, importance_threshold = 0.017, freq_selected_threshold = 0.8)


### calculate and plot logloss (plots change of logloss over boosting rounds (type = "train", "test", or "both"))
### compute logloss gap (test loss - train loss) over boosting rounds and plot
plot_logloss <- function(results_obj, type = "test", show_mean = TRUE, max_rounds = NULL) {
  
  # calculate logloss
  logloss_obj <- purrr::imap_dfr(results_obj$eval_logs, function(logs, name) {
    df <- dplyr::bind_rows(logs, .id = "Fold")
    df$repeat_fold <- name
    as.data.frame(df)
  })
  
  # plot logloss gap
  if (type == "gap") {
    
    # calculate logloss gap
    logloss_gap_df <- logloss_obj %>%
      dplyr::group_by(iter) %>%
      dplyr::summarise(mean_train = mean(train_logloss_mean, na.rm = TRUE),
                       mean_test = mean(test_logloss_mean, na.rm = TRUE),
                       gap = mean(test_logloss_mean - train_logloss_mean, na.rm = TRUE),
                       .groups = "drop")
    
    # limit number of boosting rounds displayed
    if (!is.null(max_rounds)) {
      logloss_gap_df <- logloss_gap_df %>% dplyr::filter(iter <= max_rounds)
    }
    
    # plot logloss gap
    p <- ggplot(logloss_gap_df, aes(x = iter, y = gap)) +
      geom_line(color = "black") + theme_minimal() +
      labs(title = "Mean logloss gap (test - train) across boosting rounds",
           x = "Boosting round",y = "Gap")
    
    return(p)
  }
  
  # otherwise, plot logloss 
  df <- logloss_obj
  
  # limit number of boosting rounds displayed
  if (!is.null(max_rounds)) {
    df <- df %>% dplyr::filter(iter <= max_rounds)
  }
  
  if (type == "test") {
    p <- ggplot(df, aes(x = iter, y = test_logloss_mean, group = repeat_fold)) +
      geom_line(alpha = 0.2, color = "black")
    
    if (show_mean) {
      mean_df <- df %>% dplyr::group_by(iter) %>%
        dplyr::summarise(mean_logloss = mean(test_logloss_mean, na.rm = TRUE), .groups = "drop")
      
      p <- p + geom_line(data = mean_df, aes(x = iter, y = mean_logloss), inherit.aes = FALSE,
                         color = "blue", linewidth = 1)
    }
    
  } else if (type == "train") {
    p <- ggplot(df, aes(x = iter, y = train_logloss_mean, group = repeat_fold)) +
      geom_line(alpha = 0.2, color = "red")
    
  } else if (type == "both") {
    df_long <- df %>%
      dplyr::select(iter, repeat_fold, train_logloss_mean, test_logloss_mean) %>%
      tidyr::pivot_longer(cols = c(train_logloss_mean, test_logloss_mean), names_to = "set", values_to = "logloss")
    
    # limit df_long if max_rounds is set (for type = "both")
    if (!is.null(max_rounds)) {
      df_long <- df_long %>% dplyr::filter(iter <= max_rounds)
    }
    
    p <- ggplot(df_long, aes(x = iter, y = logloss, color = set, group = interaction(repeat_fold, set))) + geom_line(alpha = 0.2) +
      scale_color_manual(values = c("train_logloss_mean" = "red", "test_logloss_mean" = "blue"))
  }
  
  p + labs(title = "Logloss curves", x = "Boosting round", y = "Logloss") + theme_minimal()
}

plot_logloss(bayes_obj, type = "test", show_mean = TRUE,  max_rounds = 1000)
plot_logloss(bayes_obj, type = "train", show_mean = FALSE, max_rounds = 1000)
plot_logloss(bayes_obj, type = "both", show_mean = FALSE,  max_rounds = 1000)
plot_logloss(bayes_obj, type = "gap",  max_rounds = 1000)


##############################################################################
###   XGBOOST LOGLOSS MODEL - TRAIN FINAL MODEL WITH BEST HYPERPARAMETERS  ###
##############################################################################

# data to be used in the model
str(metagen)

# optimum number of nrounds chosen during final evaluation
best_nrounds <- bayes_obj$per_fold$nrounds_values[1]

# best hyperparameter values determined by Bayesian optimization
best_params <- list(objective = "binary:logistic",
                    eval_metric = "logloss",
                    eta = bayes_obj$per_fold$eta[1],
                    scale_pos_weight = bayes_obj$per_fold$scale_pos_weight[1],
                    max_depth = bayes_obj$per_fold$max_depth[1],
                    min_child_weight = bayes_obj$per_fold$min_child_weight[1],
                    subsample = bayes_obj$per_fold$subsample[1],
                    colsample_bytree = bayes_obj$per_fold$colsample_bytree[1],
                    colsample_bynode = bayes_obj$per_fold$colsample_bynode[1],
                    lambda = bayes_obj$per_fold$lambda[1],
                    alpha = bayes_obj$per_fold$alpha[1],
                    gamma = bayes_obj$per_fold$gamma[1],
                    max_delta_step = bayes_obj$per_fold$max_delta_step[1])

# train the model on the full dataset
dtrain_full <- xgboost::xgb.DMatrix(data = as.matrix(metagen[, all_feat_cols]),
                                    label = metagen$condition_numeric)

final_model <- xgboost::xgb.train(params = best_params,
                                  data = dtrain_full,
                                  nrounds = best_nrounds,
                                  verbose = 1)

#################################################################################
###   XGBOOST LOGLOSS MODEL - SHAP VALUES - DEPENDENCE AND INTERACTION PLOTS  ###
#################################################################################

# compute tree SHAP values
shap_values <- predict(final_model, newdata = dtrain_full, predcontrib = TRUE)
shap_df <- as.data.frame(shap_values)
shap_df$'(Intercept)' <- NULL  # remove bias term

# prepare shap values into logn format for plotting
shap_long <- shap.prep(xgb_model = final_model, X_train = as.matrix(metagen[, all_feat_cols]))

# calculate mean absolute SHAP value per feature
shap_mean_abs <- sort(colMeans(abs(shap_df)), decreasing = TRUE)
shap_mean_abs <- as.data.frame(shap_mean_abs) # covert to data.frame
shap_mean_abs$feature <- rownames(shap_mean_abs)
shap_mean_abs <- shap_mean_abs %>% 
  arrange(desc(shap_mean_abs))

# subset shap_long to top 25 features by mean shap value
top25_feat <- shap_mean_abs$feature[1:25]
shap_long_top25 <- shap_long[variable %in% top25_feat]
shap_long_top25[, variable := droplevels(variable)] # remove extra levels

# plot summary of shap values 
shap.plot.summary(shap_long) # full dataset
shap.plot.summary(shap_long_top25) # reduced dataset


# recreate shap.plot.summary from treeshap (beeswarm-style plot) in ggpplot using reduced dataset
feature_order <- shap_mean_abs$feature
shap_plot <- shap_long_top25 %>%
  mutate(variable = factor(variable, levels = rev(feature_order)))

ggplot(shap_plot, aes(x = value, y = variable, color = rfvalue)) +
  geom_jitter(height = 0.2, alpha = 0.7, size = 1.2) +
  scale_color_viridis_c(option = "plasma", direction = -1) + theme_minimal() +
  labs(title = "SHAP summary plot", x = "SHAP value (impact on model output)", 
       color = "Feature value")


### plot mean absolute SHAP value per feature for top 25 features
shap_feat_plot_df <- head(shap_mean_abs, n = 25)
ggplot(shap_feat_plot_df, aes(x = reorder(feature, shap_mean_abs), y = shap_mean_abs)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = "Mean absolute SHAP value", 
       x = "Feature", y = "Mean absolute SHAP value")


### SHAP dependence plots (how SHAP values for a given feature vary as the input values for the feature vary)
# plot dependence plot for feature of interest (shap.plot.dependence and ggplot)
shap.plot.dependence(data_long = shap_long, x = "Lachnoclostridium_sp._YL32", y = NULL)
shap.plot.dependence(data_long = shap_long, x = "Petrimonas_mucosa", y = NULL)


### SHAP dependence plots plus interaction feature
# wide table of CLR-transformed relative abundance
rfvalue_wide <- shap_long %>%
  select(ID, variable, rfvalue) %>%
  pivot_wider(names_from = variable, values_from = rfvalue)

feature_name <- "Lachnoclostridium_sp._YL32"
shap_dep <- shap_long %>% filter(variable == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "Petrimonas_mucosa"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- CLR abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- CLR abun"))

feature_name <- "Petrimonas_mucosa"
shap_dep <- shap_long %>% filter(variable == feature_name)
plot_df <- shap_dep %>% left_join(rfvalue_wide, by = "ID")
interaction_feature <- "Acutalibacter_muris"

ggplot(plot_df, aes(x = rfvalue, y = value, color = .data[[interaction_feature]])) + theme_minimal() +
  geom_point(alpha = 0.8) + geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_color_distiller(palette = "RdBu", direction = 1) +
  labs(title = paste("SHAP dependence plot for", feature_name),
       x = paste(feature_name, "- CLR abun"), 
       y = paste(feature_name, "- SHAP value"),
       color = paste(interaction_feature, "- CLR abun"))


### SHAP interaction values (how pairs of features interact in affecting the prediction)
interaction_values <- predict(final_model,
                              newdata = as.matrix(metagen[, all_feat_cols]),
                              predinteraction = TRUE)

interaction_values <- interaction_values[, -ncol(interaction_values), -ncol(interaction_values)] # remove BIAS term
mean_interactions <- apply(abs(interaction_values), c(2, 3), mean) # average absolute interaction strengths

# set row and column names
feature_names <- colnames(metagen[, all_feat_cols])
rownames(mean_interactions) <- feature_names
colnames(mean_interactions) <- feature_names

# convert to long format for plotting
interaction_long <- as.data.frame(mean_interactions) %>%
  rownames_to_column("Feature1") %>%
  pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "InteractionStrength") 

# remove self-interactions (diagonal) - to see other interactions more clearly
interaction_long <- interaction_long %>%
  filter(Feature1 != Feature2)

# identify top 25 strongest interaction pairs
top25_pairs <- interaction_long %>%
  dplyr::arrange(dplyr::desc(InteractionStrength)) %>%
  head(n = 25)

# get all unique features involved in those pairs
top_feats_from_pairs <- unique(c(top25_pairs$Feature1, top25_pairs$Feature2))

# subset the original full interaction matrix to those features
mean_interactions_sub <- mean_interactions[top_feats_from_pairs, top_feats_from_pairs]

# zero out self-interactions for clearer color scale
diag(mean_interactions_sub) <- NA 

# convert subset mean_interactions to long format for plotting
interaction_sub_long <- as.data.frame(mean_interactions_sub) %>%
  tibble::rownames_to_column("Feature1") %>%
  tidyr::pivot_longer(cols = -Feature1, names_to = "Feature2", values_to = "InteractionStrength")

# plot interaction heatmap
ggplot(interaction_sub_long, aes(x = Feature1, y = Feature2, fill = InteractionStrength)) +
  geom_tile() + scale_fill_viridis_c() + theme_minimal() + coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "SHAP feature interactions", fill = "Mean absolute interaction")


### main (additive) effects versus interaction (synergy/modulation) effects
# calculate main SHAP effect
main_effects <- sapply(1:length(feature_names), function(i) {
  mean(abs(interaction_values[, i, i]))
})

# calculate strength of pairwise interactions
interaction_strength <- sapply(1:length(feature_names), function(i) {
  mean(rowSums(abs(interaction_values[, i, -i])))
})

effect_df <- data.frame(feature = feature_names,
                        main = main_effects,
                        interaction = interaction_strength)

effect_df$interaction_ratio <- effect_df$interaction / (effect_df$main + effect_df$interaction)

# plot how much the taxon influences prediction through interactions vs how strongly the taxon independently influences prediction
ggplot(effect_df %>% dplyr::filter(main > 0.025), aes(x = main, y = interaction)) +
  geom_point(color = "steelblue", size = 3) + theme_minimal() +
  geom_text_repel(aes(label = feature), size = 3, max.overlaps = 20) +
  labs(x = "Mean |main effect|", y = "Total mean |interaction effect|")


### feature-specific interactions
# target_feat <- "Acutalibacter_muris" # feature of interest
# target_feat <- "Ruminococcus_gauvreauii" # feature of interest
# target_feat <- "Enterocloster_clostridioformis" # feature of interest
# target_feat <- "Anaerostipes_caccae" # feature of interest
target_feat <- "Petrimonas_mucosa" # feature of interest

# data.frame with average feature specific interaction for target_feat
target_interact <- data.frame(feature = feature_names,
                              interaction = mean_interactions[which(feature_names == target_feat), ])
target_interact <- target_interact[target_interact$feature != target_feat, ] # remove self-interaction

# remove zero interactions
target_interaction_sub <- target_interact %>% 
  filter(interaction > 0) %>%
  arrange(desc(interaction))

# plot feature-specific interactions
ggplot(target_interaction_sub[1:10, ], aes(x = reorder(feature, interaction), y = interaction)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_minimal() +
  labs(title = paste("SHAP interactions with", target_feat),
       x = "Interacting feature", y = "Mean absolute interaction")


### mean absolute SHAP value per class per feature
shap_df$condition <- metagen$condition
abs_shap_by_class <- shap_df %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop")

# plot mean absolute SHAP value per class per feature
abs_shap_by_class_top25 <- abs_shap_by_class %>%
  filter(name %in% top25_feat) # filter features with top 25 abs shap values

ggplot(abs_shap_by_class_top25, aes(x = reorder(name, mean_abs_shap), y = mean_abs_shap, fill = condition)) +
  geom_col(position = "dodge") + coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean absolute SHAP values",
       x = "Feature", y = "Mean abs SHAP", fill = "Condition")


### mean SHAP value per class per feature
shap_df$condition <- metagen$condition
shap_by_class <- shap_df %>%
  pivot_longer(cols = -condition) %>%
  group_by(condition, name) %>%
  summarise(mean_shap = mean(value), .groups = "drop")

# plot mean SHAP value per class per feature
shap_by_class_top25 <- shap_by_class %>%
  filter(name %in% top25_feat) # filter features with top 25 abs shap values

ggplot(shap_by_class_top25, aes(x = reorder(name, mean_shap), y = mean_shap, fill = condition)) +
  geom_col(position = "dodge") + coord_flip() + theme_minimal() +
  labs(title = "Class-specific mean SHAP values",
       x = "Feature", y = "Mean SHAP", fill = "Condition")


### how the impact of a given feature on model prediction varies between healthy and disease samples
# histogram of SHAP values
ggplot(shap_df, aes(x = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - Lachnoclostridium_sp._YL32",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_df, aes(x = condition, y = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for Lachnoclostridium_sp._YL32 by condition",
       y = "SHAP value", x = "Condition")

# box plot of CLR-abundance
ggplot(metagen, aes(x = condition, y = Lachnoclostridium_sp._YL32, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "CLR abundance of Lachnoclostridium_sp._YL32 by condition",
       y = "CLR Abundance", x = "Condition")


# histogram of SHAP values
ggplot(shap_df, aes(x = Petrimonas_mucosa, fill = condition)) +
  geom_density(alpha = 0.6) + theme_minimal() +
  labs(title = "SHAP value distribution - Petrimonas_mucosa",
       x = "SHAP value", y = "Density", fill = "Condition")

# boxplot of SHAP values
ggplot(shap_df, aes(x = condition, y = Petrimonas_mucosa, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "SHAP values for Petrimonas_mucosa by condition",
       y = "SHAP value", x = "Condition")

# box plot of CLR-abundance
ggplot(metagen, aes(x = condition, y = Petrimonas_mucosa, fill = condition)) +
  geom_boxplot() + theme_minimal() +
  labs(title = "CLR abundance of Petrimonas_mucosa by condition",
       y = "CLR abundance", x = "Condition")


### SHAP values versus predicted probabilities
shap_df$ID <- rownames(shap_df)
shap_df$pred_prob <- predict(final_model, newdata = dtrain_full)

# Petrimonas_mucosa
plot_df <- shap_df %>%
  select(ID, Petrimonas_mucosa, condition, pred_prob)

ggplot(plot_df, aes(x = Petrimonas_mucosa, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Petrimonas_mucosa SHAP value versus prediction probability",
       x = "SHAP value for Petrimonas_mucosa", y = "Predicted probability (disease)")

# Lachnoclostridium_sp._YL32
plot_df <- shap_df %>%
  select(ID, Lachnoclostridium_sp._YL32, condition, pred_prob)

ggplot(plot_df, aes(x = Lachnoclostridium_sp._YL32, y = pred_prob, color = condition)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  labs(title = "Lachnoclostridium_sp._YL32 SHAP value versus prediction probability",
       x = "SHAP value for Lachnoclostridium_sp._YL32", y = "Predicted probability (disease)")


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

