# Logistic Regression Workflow - Hyperparameter Tuning

# load libraries
library(ggplot2)
library(tidyverse)
library(compositions)
library(glmnet)
library(caret)
library(pROC)
library(MLmetrics)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# metadata
meta <- read.table("metadata.txt", header = TRUE)
rownames(meta) <- meta$sample_name # add sample_names as rownames
meta$condition <- as.factor(meta$condition) # convert condition column into a factor
table(meta$condition)


# feature table
feat <- read.table("feature_table.txt", header = TRUE)

# filter species present in less than 10% of samples
dim(feat) # 2302   70
presence_threshold <- 0.10 * ncol(feat) # needs to be in at least 7 samples
features_to_keep <- rowSums(feat != 0) >= presence_threshold # determine if the feature is present >=10% of samples
feat_filt <- feat[features_to_keep, ] # subset the feature table to only include features present in at least 10% of samples
dim(feat_filt) # 935  70

# convert feature table to relative abundances
feat_rel <- sweep(feat_filt, 2, colSums(feat_filt), FUN = "/")
feat_rel <- as.data.frame(feat_rel)

# CLR transformation
feat_rel_pc <- feat_rel + 1e-6 # add pseudocount
feat_clr <- clr(t(feat_rel_pc)) # transpose data and perform CLR transformation
feat_final <- as.data.frame(feat_clr)

# rownames for metadata need to match rownames for feature table
all(rownames(meta) == rownames(feat_final))

# make sure feature names are syntactically valid 
colnames(feat_final) <- make.names(colnames(feat_final))


############################################################################
###   LOGISTIC REGRESSION + CROSS-VALIDATION - ALPHA AND LAMBDA TUNING   ###
############################################################################

# set.seed
set.seed(1234)

# prepare metadata and feature table for use with glmnet
feat <- as.matrix(feat_final)
label <- meta$condition

# set lambda and alpha values to test
alpha_values <- seq(0, 1, by = 0.2)
lambda_values <- 10^seq(-4, 1, length.out = 20) # logarithmically spaced

# create list to store performance metrics
performance_metrics <- list()

# loop through alpha values
for (a in alpha_values){
  
  # loop through lambda values
  for (l in lambda_values){
    cat("Alpha =", a, " Lambda =", l, "\n")
    
    # repeat cross-validation 50 times
    for (r in 1:50){
      
      # create stratified folds
      folds <- createFolds(label, k = 5, returnTrain = TRUE)
      
      # loop through folds
      for (f in 1:5) {
        
        # split the dataset into training and testing sets for the current fold
        train_idx <- folds[[f]] # train indices for the f-th fold
        test_idx <- setdiff(seq_along(label), train_idx)
        
        feat_train <- feat[train_idx, ] # training data
        label_train <- label[train_idx] # training labels
        feat_test <- feat[test_idx, ] # testing data
        label_test <- label[test_idx] # testing labels
        
        # fit logistic regression with alpha and lambda (binary classification)
        model <- glmnet(feat_train, label_train, family = "binomial",
                        alpha = a, lambda = l, standardize = TRUE)
        
        # evaluate on test set
        probabilities <- as.vector(predict(model, feat_test, type = "response"))
        predictions <- ifelse(probabilities > 0.5, "disease", "healthy")
        predictions <- factor(predictions, levels = levels(label))
        
        # generate confusion matrix
        cm <- confusionMatrix(predictions, label_test, positive = "disease")
        
        # calculate AUC
        label_test <- factor(label_test, levels = c("healthy", "disease"))
        roc_obj <- roc(response = label_test, predictor = probabilities, quiet = TRUE)
        auc_value <- auc(roc_obj)
        
        # store with alpha (a), lambda (l), repeat (r) and fold (f) index
        key <- paste0("alpha_", a, "_lambda_", l, "_Repeat_", r, "_Fold_", f)
        performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
      }
    }
  }
}


### calculate performance statistics
# create vectors to store metrics
alpha_lambda_repeat <- character()
alpha_lambda <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
precision <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract alpha, lambda and repeat name
  alr <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+")
  al <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+")
  
  alpha_lambda_repeat <- c(alpha_lambda_repeat, alr)
  alpha_lambda <- c(alpha_lambda, al)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  precision <- c(precision, cm$byClass["Precision"])
  auc_vals <- c(auc_vals, auc)
}


# combine metrics in a summary table
results_df <- data.frame(alpha_lambda_repeat = alpha_lambda_repeat,
                         alpha_lambda = alpha_lambda,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         prec = precision,
                         auc = auc_vals)

# summary of performance metrics (mean values for each repeat)
metric_summary <- results_df %>%
  group_by(alpha_lambda_repeat) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            mean_prec = mean(prec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(alpha_lambda_repeat)
metric_summary <- as.data.frame(metric_summary)


# summary of performance metrics (mean and standard deviation values for each alpha + lambda combo)
overall_metric_summary <- results_df %>%
  group_by(alpha_lambda) %>%
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
overall_metric_summary <- as.data.frame(overall_metric_summary)


############################################################################################
###   LOGISTIC REGRESSION + CROSS-VALIDATION - ALPHA AND LAMBDA TUNING + CLASS WEIGHTS   ###
############################################################################################

# set.seed
set.seed(1234)

# prepare metadata and feature table for use with glmnet
feat <- as.matrix(feat_final)
label <- meta$condition # labels moderately imbalanced

# set lambda and alpha values to test
alpha_values <- seq(0, 1, by = 0.2)
lambda_values <- 10^seq(-4, 1, length.out = 20) # logarithmically spaced

# create vector of sample weights (each observation gets a weight inversely proportional to its class frequency)
class_weights <- table(label)
weight_vector <- ifelse(label == "disease", 
                        1/class_weights["disease"], 1/class_weights["healthy"])
# normalize class weights to sum of number of observations
# keeps total loss scale equivalent to unweighted training, maintains comparability of penalization and regularization
weight_vector <- weight_vector * length(label)/sum(weight_vector) 

# create list to store performance metrics
performance_metrics <- list()

# loop through alpha values
for (a in alpha_values){
  
  # loop through lambda values
  for (l in lambda_values){
    cat("Alpha =", a, " Lambda =", l, "\n")
    
    # repeat cross-validation 50 times
    for (r in 1:50){
      
      # create stratified folds
      folds <- createFolds(label, k = 5, returnTrain = TRUE)
      
      # loop through folds
      for (f in 1:5) {
        
        # split the dataset into training and testing sets for the current fold
        train_idx <- folds[[f]] # train indices for the f-th fold
        test_idx <- setdiff(seq_along(label), train_idx)
        
        feat_train <- feat[train_idx, ] # training data
        label_train <- label[train_idx] # training labels
        feat_test <- feat[test_idx, ] # testing data
        label_test <- label[test_idx] # testing labels
        
        # extract weights for the training samples only (re-index for each training split)
        train_weights <- weight_vector[train_idx]
        
        # fit logistic regression with alpha and lambda and class weighting
        model <- glmnet(feat_train, label_train, family = "binomial",
                        alpha = a, lambda = l, standardize = TRUE,
                        weights = train_weights)
        
        # evaluate on test set
        probabilities <- as.vector(predict(model, feat_test, type = "response"))
        predictions <- ifelse(probabilities > 0.5, "disease", "healthy")
        predictions <- factor(predictions, levels = levels(label))
        
        # generate confusion matrix
        cm <- confusionMatrix(predictions, label_test, positive = "disease")
        
        # calculate AUC
        label_test <- factor(label_test, levels = c("healthy", "disease"))
        roc_obj <- roc(response = label_test, predictor = probabilities, quiet = TRUE)
        auc_value <- auc(roc_obj)
        
        # store with alpha (a), lambda (l), repeat (r) and fold (f) index
        key <- paste0("alpha_", a, "_lambda_", l, "_Repeat_", r, "_Fold_", f)
        performance_metrics[[key]] <- list(cm = cm, auc = auc_value) # store performance metrics
      }
    }
  }
}


### calculate performance statistics
# create vectors to store metrics
alpha_lambda_repeat <- character()
alpha_lambda <- character()
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()
precision <- numeric()
auc_vals <- numeric()

# loop through each stored confusion matrix + name
for (key in names(performance_metrics)) {
  result <- performance_metrics[[key]]
  cm <- result$cm
  auc <- result$auc
  
  # extract alpha, lambda and repeat name
  alr <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_[^_]+")
  al <- str_extract(key, "^[^_]+_[^_]+_[^_]+_[^_]+")
  
  alpha_lambda_repeat <- c(alpha_lambda_repeat, alr)
  alpha_lambda <- c(alpha_lambda, al)
  balanced_accuracy <- c(balanced_accuracy, cm$byClass["Balanced Accuracy"])
  f1_score <- c(f1_score, cm$byClass["F1"])
  sensitivity <- c(sensitivity, cm$byClass["Sensitivity"])
  specificity <- c(specificity, cm$byClass["Specificity"])
  precision <- c(precision, cm$byClass["Precision"])
  auc_vals <- c(auc_vals, auc)
}


# combine metrics in a summary table
results_df <- data.frame(alpha_lambda_repeat = alpha_lambda_repeat,
                         alpha_lambda = alpha_lambda,
                         bal_acc = balanced_accuracy,
                         f1 = f1_score,
                         sens = sensitivity,
                         spec = specificity,
                         prec = precision,
                         auc = auc_vals)

# summary of performance metrics (mean values for each repeat)
metric_summary <- results_df %>%
  group_by(alpha_lambda_repeat) %>%
  summarise(mean_bal_acc = mean(bal_acc, na.rm = TRUE),
            mean_f1 = mean(f1, na.rm = TRUE),
            mean_sens = mean(sens, na.rm = TRUE),
            mean_spec = mean(spec, na.rm = TRUE),
            mean_prec = mean(prec, na.rm = TRUE),
            mean_auc = mean(auc, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(alpha_lambda_repeat)
metric_summary <- as.data.frame(metric_summary)


# summary of performance metrics (mean and standard deviation values for each alpha + lambda combo)
overall_metric_summary <- results_df %>%
  group_by(alpha_lambda) %>%
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
overall_metric_summary <- as.data.frame(overall_metric_summary)


### determine the alpha and lambda values from the best performing model
best_params <- overall_metric_summary %>%
  slice(1)

best_alpha <- as.numeric(str_extract(best_params$alpha_lambda, "(?<=alpha_)\\d+\\.?\\d*"))
best_lambda <- as.numeric(str_extract(best_params$alpha_lambda, "(?<=lambda_)\\d+\\.?\\d*"))

### refit glmnet on the full dataset using the best parameters
final_model <- glmnet(x = feat, y = label, family = "binomial",
                      alpha = best_alpha, lambda = best_lambda,
                      weights = weight_vector,
                      standardize = TRUE) 

# extract coefficients (feature importance)
coefs <- coef(final_model)
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "coefficient"

# remove the intercept and zero coefficients
# rank by effect size
nonzero_coefs <- coefs_df %>%
  filter(feature != "(Intercept)" & coefficient != 0) %>%
  arrange(desc(abs(coefficient)))


###########################################################
###   LOGISTIC REGRESSION - CHECK WITH SIMPLER MODELS   ###
###########################################################

# set.seed
set.seed(1234)

# prepare metadata and feature table
feat <- as.matrix(feat_final)
label <- meta$condition 


### one stratified train/test split (no CV)
# stratified sampling: 80% train and 20% test
train_idx <- createDataPartition(label, p = 0.8, list = FALSE)
feat_train <- feat[train_idx, ]
feat_test  <- feat[-train_idx, ]
label_train <- label[train_idx]
label_test  <- label[-train_idx]

# fit basic logistic regression model (without regularization) with glm
df_train <- as.data.frame(feat_train) # convert to data.frame for glm
df_train$label <- label_train
glm_model <- glm(label ~ ., data = df_train, family = "binomial")

# evaluate on test set
probs_glm <- as.vector(predict(glm_model, newdata = as.data.frame(feat_test), type = "response"))
preds_glm <- ifelse(probs_glm > 0.5, "disease", "healthy")
preds_glm <- factor(preds_glm, levels = levels(label_test))

# generate confusion matrix
cm_glm <- confusionMatrix(preds_glm, label_test, positive = "disease")

# calculate AUC
label_test <- factor(label_test, levels = c("healthy", "disease"))
roc_glm <- roc(label_test, probs_glm, quiet = TRUE)
auc_glm <- auc(roc_glm)

# performance metrics
cm_glm # confusion matrx
auc_glm # AUC


### fit regularized model with cv.glmnet
cv_model <- cv.glmnet(feat_train, label_train, family = "binomial", alpha = 0.5)

# evaluate on test set using best lambda
probs_glmnet <- as.vector(predict(cv_model, feat_test, type = "response", s = "lambda.min"))
preds_glmnet <- ifelse(probs_glmnet > 0.5, "disease", "healthy")
preds_glmnet <- factor(preds_glmnet, levels = levels(label_test))

# generate confusion matrix
cm_glmnet <- confusionMatrix(preds_glmnet, label_test, positive = "disease")

# calculate AUC
label_test <- factor(label_test, levels = c("healthy", "disease"))
roc_glmnet <- roc(label_test, probs_glmnet, quiet = TRUE)
auc_glmnet <- auc(roc_glmnet)

# performance metrics 
cm_glmnet # confusion matrix
auc_glmnet # AUC


### data is not linearly separable


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
#   [1] MLmetrics_1.1.3    pROC_1.18.5        caret_7.0-1        lattice_0.22-7     glmnet_4.1-9      
# [6] Matrix_1.7-3       compositions_2.0-8 lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [11] dplyr_1.1.4        purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0      
# [16] tidyverse_2.0.0    ggplot2_3.5.2     
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6         shape_1.4.6.1        tensorA_0.36.2.1     recipes_1.3.1       
# [5] tzdb_0.5.0           vctrs_0.6.5          tools_4.5.0          generics_0.1.4      
# [9] stats4_4.5.0         parallel_4.5.0       proxy_0.4-27         DEoptimR_1.1-3-1    
# [13] ModelMetrics_1.2.2.2 pkgconfig_2.0.3      data.table_1.17.4    RColorBrewer_1.1-3  
# [17] lifecycle_1.0.4      compiler_4.5.0       farver_2.1.2         codetools_0.2-20    
# [21] class_7.3-23         prodlim_2025.04.28   pillar_1.10.2        MASS_7.3-65         
# [25] gower_1.0.2          iterators_1.0.14     rpart_4.1.24         foreach_1.5.2       
# [29] parallelly_1.45.0    nlme_3.1-168         lava_1.8.1           robustbase_0.99-4-1 
# [33] tidyselect_1.2.1     digest_0.6.37        stringi_1.8.7        future_1.58.0       
# [37] reshape2_1.4.4       listenv_0.9.1        splines_4.5.0        grid_4.5.0          
# [41] cli_3.6.5            magrittr_2.0.3       dichromat_2.0-0.1    survival_3.8-3      
# [45] e1071_1.7-16         future.apply_1.20.0  withr_3.0.2          scales_1.4.0        
# [49] timechange_0.3.0     globals_0.18.0       nnet_7.3-20          timeDate_4041.110   
# [53] hms_1.1.3            hardhat_1.4.1        rlang_1.1.6          Rcpp_1.0.14         
# [57] glue_1.8.0           bayesm_3.1-6         ipred_0.9-15         rstudioapi_0.17.1   
# [61] R6_2.6.1             plyr_1.8.9 

