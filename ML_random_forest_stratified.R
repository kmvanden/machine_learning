# Random Forest Workflow - Stratification by Metadata Variables

# load libraries
library(curatedMetagenomicData)
library(ggplot2)
library(tidyverse)
library(randomForest)
library(caret)
library(vegan)
library(mice)
library(e1071)
library(VennDiagram)


# set.seed
set.seed(1234)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data from curatedMetagenomicData
ibdmdb <- curatedMetagenomicData("2021-10-14.HMP_2019_ibdmdb.relative_abundance", dryrun = FALSE, counts = TRUE, rownames = "short")
ibd <- ibdmdb[[1]]
str(ibd)


### METADATA
# extract metadata
ibd_meta <- colData(ibd)
ibd_meta <- as.data.frame(ibd_meta)
ibd_meta$sample_id <- rownames(ibd_meta)

# keep only one sample from the first visit for each subject
table(ibd_meta$visit_number)
length(unique(ibd_meta$subject_id))

ibd_meta_first <- ibd_meta %>%
  filter(visit_number == 1) # only keep samples from the first visit

ibd_meta_filt <- ibd_meta_first %>%
  group_by(subject_id) %>%
  filter(n() == 1 | grepl("_P$", sample_id)) %>%
  ungroup() # only keep one sample per subject from the first visit 
ibd_meta_filt <- as.data.frame(ibd_meta_filt)
rownames(ibd_meta_filt) <- ibd_meta_filt$sample_id
head(ibd_meta_filt)

# remove metadata columns that only have one value
ibd_meta_filt <- ibd_meta_filt[, sapply(ibd_meta_filt, function(col) length(unique(col)) > 1)]

# remove metadata categories that won't be used in model fitting
# study_condition == disease
ibd_meta_filt <- ibd_meta_filt[, !(colnames(ibd_meta_filt) %in% c("study_condition", "PMID", "number_reads", 
                                                                  "number_bases", "minimum_read_length", 
                                                                  "median_read_length", "sample_id", 
                                                                  "age_category", "subject_id"))]
# change NAs to "healthy" in disease_subtype
ibd_meta_filt$disease_subtype[is.na(ibd_meta_filt$disease_subtype)] <- "healthy"

# shorten column name
colnames(ibd_meta_filt)[which(colnames(ibd_meta_filt) == "antibiotics_current_use")] <- "antibiotics"

# format variables correctly for modeling (as factor or numeric)
head(ibd_meta_filt)
ibd_meta_filt <- ibd_meta_filt %>%
  mutate(across(c(antibiotics, disease, gender, location, disease_subtype), as.factor)) # change variables to factors
ibd_meta_filt$age <- as.numeric(ibd_meta_filt$age) # change variable to numeric
str(ibd_meta_filt)

# add sample name column to metadata
ibd_meta_filt$sample_name <- rownames(ibd_meta_filt)


### FEATURE TABLE
# extract feature table
ibd_feat <- assay(ibd, "relative_abundance")
ibd_feat <- as.data.frame(ibd_feat)

# format rownames
rownames(ibd_feat)[1:10]
rownames(ibd_feat) <- rownames(ibd_feat) %>%
  sub("^species:", "", .) %>%
  gsub("\\[|\\]", "", .) %>%
  gsub(" sp\\. ", "_", .) %>%
  gsub(" ", "_", .) %>%
  gsub(":", "_", .)
rownames(ibd_feat)[1:10]

# subset feature table to samples that were retained in the metadata (one sample per subject)
dim(ibd_feat) # 585 1627
all(rownames(ibd_meta_filt) %in% colnames(ibd_feat))
ibd_feat_filt <- ibd_feat[, colnames(ibd_feat) %in% rownames(ibd_meta_filt)]
all(colnames(ibd_feat_filt) == rownames(ibd_meta_filt)) # column names of feature data should exactly match rownames of metadata
dim(ibd_feat_filt) # 585 130 (from 1627 samples to 130 --> one sample for each subject)

# convert feature table to relative abundances
ibd_feat_rel <- sweep(ibd_feat_filt, 2, colSums(ibd_feat_filt), FUN = "/")
ibd_feat_rel <- as.data.frame(ibd_feat_rel)

### add pseudocount and perform log-transform
log_n0 <- 1e-6 # pseudocount
ibd_feat_log <- t(ibd_feat_rel) # transpose feature table
ibd_feat_log <- log(ibd_feat_log + log_n0)
ibd_feat_log <- as.data.frame(ibd_feat_log)

# perform row-wise L2 normalization (SIAMCAT)
n_p <- 2 # L2 norm
ibd_row_norms <- sqrt(rowSums(ibd_feat_log^n_p))
ibd_feat_lognorm <- sweep(ibd_feat_log, 1, ibd_row_norms, FUN = "/")
ibd_feat <- as.data.frame(ibd_feat_lognorm)


# rownames of metadata need to match the rownames of the feature table
all(rownames(ibd_meta_filt) == rownames(ibd_feat))
ibd_feat$sample_name <- rownames(ibd_feat) # add sample_name column to feature table

### merge metadata and feature table
ibd_merge <- merge(ibd_meta_filt, ibd_feat, by = "sample_name", all.x = TRUE)
ibd_merge <- ibd_merge[,-1] # remove sample_name
ibd_merge[1:10, 1:10]


############################################
###   HANDLING MISSING METADATA VALUES   ###
############################################

### MISSING VALUES: BMI
# over 1/3 of samples do not have a BMI measurement
sum(is.na(ibd_merge$BMI))/length(ibd_merge$BMI) # 0.3692308

# create column indicating whether BMI measurement is present
ibd_merge$BMI_avail <- !is.na(ibd_merge$BMI)


### BMI_avail significantly associated with other metadata variables
# antibiotics
table(ibd_merge$antibiotics, ibd_merge$BMI_avail)
chisq.test(table(ibd_merge$antibiotics, ibd_merge$BMI_avail)) # p-value = 1

# disease
table(ibd_merge$disease, ibd_merge$BMI_avail)
chisq.test(table(ibd_merge$disease, ibd_merge$BMI_avail)) # p-value = 0.2686

# age
t.test(age ~ BMI_avail, data = ibd_merge) # p-value = 0.3145

# gender
table(ibd_merge$gender, ibd_merge$BMI_avail)
chisq.test(table(ibd_merge$gender, ibd_merge$BMI_avail)) # p-value = 0.9621

# location
table_location <- table(ibd_merge$location, ibd_merge$BMI_avail)
fisher.test(table_location) # p-value = 0.0001104

# disease subtype
table_disease_subtype <- table(ibd_merge$disease_subtype, ibd_merge$BMI_avail)
fisher.test(table_disease_subtype) # p-value = 0.3981

### BMI_avail is significantly associated with location (but not other metadata variables)


### BMI_avail significantly associated with microbiome structure 
# feature table should only include microbiome features
ibd_feat_adonis <- ibd_feat # use transposed feature table used in merge
ibd_feat_adonis <- ibd_feat_adonis[,-586] # removed sample_name

# metadata should only contain metadata
ibd_meta_adonis <- ibd_meta_filt[, -8] # remove sample name
ibd_meta_adonis$BMI_avail <- ibd_merge$BMI_avail # add BMI_avail
  
# rownames need to match 
all(rownames(ibd_feat_adonis) == rownames(ibd_meta_adonis))

# PERMANOVA using euclidean method (to handle the log unit normalization)
adonis2(ibd_feat_adonis ~ BMI_avail, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.00766 | p-value: 0.435

### BMI_avail is not significantly associated with the microbiome structure

### remove BMI as a metadata variable
ibd_merge <- ibd_merge[, -c(6, 593)] # remove BMI and BMI_avail from merged data.frame (to be used in modeling)
ibd_meta_adonis <- ibd_meta_adonis[, -c(6, 8)] # remove BMI and BMI_avail from ibd_meta_adonis (used for PERMANOVA)


### MISSING VALUES: AGE
# missing metadata values for age (<5% of samples don't have a value for age)
sum(is.na(ibd_merge$age))/length(ibd_merge$age) # 0.04615385

### impute missing values using MICE
# subset metadata columns for imputation (age + variables that could help impute age values)
meta_impute <- ibd_merge[, c("age", "antibiotics", "disease", "gender", "location", "disease_subtype")]

# impute age values using MICE and pmm method (predictive mean matching for numeric variables)
imputed_data <- mice(meta_impute, m = 50, method = "pmm", maxit = 10, seed = 123) 

# run linear regression on each imputed dataset and pool results
fit <- with(imputed_data, lm(age ~ antibiotics + disease + gender + location + disease_subtype))
pooled_results <- pool(fit)
pooled_results # check pooled results

# extract completed datasets into a list
imputed_list <- lapply(1:50, function(i) mice::complete(imputed_data, i))

# average imputed ages across all imputations
age_matrix <- sapply(imputed_list, function(df) df$age)
mean_imputed_age <- rowMeans(age_matrix)
mean_imputed_age <- round(mean_imputed_age)

# replace age with age_imputed (NAs in age column replaced by average imputed age values)
ibd_merge$age <- mean_imputed_age
colnames(ibd_merge)[colnames(ibd_merge) == "age"] <- "age_imputed" # change name to age_imputed

ibd_meta_adonis$age <- mean_imputed_age
colnames(ibd_meta_adonis)[colnames(ibd_meta_adonis) == "age"] <- "age_imputed" # change name to age_imputed


##########################################
###   CHECKING METADATA ASSOCIATIONS   ###
##########################################

# identify and account for confounders

### are metadata variables associated with the microbiome structure
# ibd_feat_adonis (feature data (numeric))
# ibd_meta_adonis (metadata variables (antibiotics, disease, age, gender, location, disease_subtype and age_imputed)
all(rownames(ibd_feat_adonis) == rownames(ibd_meta_adonis)) # rownames (sample_ids need to match)

# antibiotics
adonis2(ibd_feat_adonis ~ antibiotics, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.00949 | p-value = 0.158

# imputed age
adonis2(ibd_feat_adonis ~ age_imputed, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.01484 | p-value = 0.019

# gender
adonis2(ibd_feat_adonis ~ gender, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.01642 | 0.006

# location
adonis2(ibd_feat_adonis ~ location, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.02813 | p-value = 0.094

# disease
adonis2(ibd_feat_adonis ~ disease, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.01153 | p-value = 0.054

# disease_subtype
adonis2(ibd_feat_adonis ~ disease_subtype, data = ibd_meta_adonis, method = "euclidean") # R2 = 0.02124 | p-value = 0.049

### use disease as label for modeling
# perform subgroup analysis with disease_subtype later


### are metadata variables significantly different between healthy and IBD patients
# ibd_meta_adonis: contains only the metadata values (not the features)

# imputed age
kruskal.test(age_imputed ~ disease, data = ibd_meta_adonis) # p-value = 0.9519

# antibiotics
fisher.test(table(ibd_meta_adonis$antibiotics, ibd_meta_adonis$disease)) # p-value = 0.1215

# gender 
# >5 values for each category --> chisq.test
chisq.test(table(ibd_meta_adonis$gender, ibd_meta_adonis$disease)) # p-value = 0.7319

# location
fisher.test(table(ibd_meta_adonis$location, ibd_meta_adonis$disease)) # p-value = 0.005405


### retain disease, disease_subtype, age_imputed, gender and location in the metadata
# location is significantly different between healthy and IBD patients and slightly associated with the microbiome structure (possible confounder)
# age_imputed and gender significantly associated with the microbiome structure (predictive signal)
   # disease and disease_subtype are slightly associated with microbiome structure
# antibiotics is not significantly associated with disease_subtype or with microbiome structure

ibd_merge <- ibd_merge[, -1] # remove antibiotics from merged data set

### check visually for batch effects using NMDS 
# calculate dissimilarity matrix (Euclidean)
dist_matrix <- dist(ibd_feat_adonis, method = "euclidean")

# perform NMDS
nmds <- metaMDS(dist_matrix, k = 2, trymax = 100)

# prepare NMDS data for plotting
nmds_df <- as.data.frame(nmds$points)
nmds_df$location <- ibd_merge$location
nmds_df$disease <- ibd_merge$disease

# plot NMDS
ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = location, shape = disease)) +
  geom_point(size = 3) + theme_minimal() +
  labs(title = "NMDS of microbiome structure (Euclidean)")


##################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - DISEASE   ###
##################################################################################

# set seed
set.seed(1234)

# column names for microbiome features
micro_feat_cols <- setdiff(colnames(ibd_merge), c("disease", "disease_subtype", "age_imputed", "gender", "location"))

# column names for microbiome features + relevant metadata (full predictor set)
all_feat_cols <- c(micro_feat_cols, "age_imputed", "gender", "location")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on disease)
  folds <- createFolds(ibd_merge$disease, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- ibd_merge[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- ibd_merge[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by all_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, all_feat_cols], 
                             y = as.factor(train_data$disease), 
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
    cm <- confusionMatrix(predictions, as.factor(test_data$disease), positive = "IBD")
    
    # store with repeat (r) and fold (f) index
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature freqeuncies
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
# mean_MeanDecreaseAccuracy: overall importance of feature on model accuracy
# mean_MeanDecreaseGini: frequency and usefulness in splitting (how much a feature reduces impurity when used to split the decision trees)
   # Gini is sensitive to splits, NOT predictive value
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_IBD = mean(IBD, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseGini))
head(mean_importance, 10)


# same data from full model for comparison
full_model_feature_split_summary <- feature_split_summary
full_model_metric_summary <- metric_summary
full_model_mean_importance <- mean_importance


### age_imputed and gender both have negative mean_MeanDecreaseAccuracy values (but age_imputed has a high mean_MeanDecreaseGini value)
# age_imputed is useful for splitting the microbiome data, but not useful for predicting for disease classification
# the association of age_imputed and gender with the microbiome is likely orthogonal to disease


##########################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - DISEASE_SUBTYPE   ###
##########################################################################################

# multi-class classification

# set seed
set.seed(1234)

# column names for microbiome features
micro_feat_cols <- setdiff(colnames(ibd_merge), c("disease", "disease_subtype", "age_imputed", "gender", "location"))

# column names for microbiome features + relevant metadata (full predictor set)
all_feat_cols <- c(micro_feat_cols, "age_imputed", "gender", "location")

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on disease_subtype)
  folds <- createFolds(ibd_merge$disease_subtype, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- ibd_merge[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- ibd_merge[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by all_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, all_feat_cols], 
                             y = as.factor(train_data$disease_subtype), 
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
    cm <- confusionMatrix(predictions, as.factor(test_data$disease_subtype))
    
    # store with repeat (r) and fold (f) index
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature freqeuncies
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


### calculate performance statistics (multi-class)
# each confusion matrix now returns a matrix with cm$byClass
# get macro-averaged metrics across the 3 classes (healthy, CD, UC)
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()

# extract metrics from the stored confusion matrices (50 repeats x 5 folds = 250 values)
for (cm in performance_metrics) {
  if (is.matrix(cm$byClass)) {
  balanced_accuracy <- c(balanced_accuracy, mean(cm$byClass[,"Balanced Accuracy"], na.rm = TRUE))
  f1_score <- c(f1_score, mean(cm$byClass[,"F1"], na.rm = TRUE))
  sensitivity <- c(sensitivity, mean(cm$byClass[,"Sensitivity"], na.rm = TRUE))
  specificity <- c(specificity, mean(cm$byClass[,"Specificity"], na.rm = TRUE))
  }
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
# mean_MeanDecreaseAccuracy: overall importance of feature on model accuracy
# mean_MeanDecreaseGini: frequency and usefulness in splitting (how much a feature reduces impurity when used to split the decision trees)
# Gini is sensitive to splits, NOT predictive value
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_CD = mean(CD, na.rm = TRUE),
            mean_UC = mean(UC, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseGini))
head(mean_importance, 10)


# same data from full model for comparison
disease_sub_feature_split_summary <- feature_split_summary
disease_sub_metric_summary <- metric_summary
disease_sub_mean_importance <- mean_importance


### comparison of disease versus disease_subtype as model label
full_model_feature_split_summary # disease
disease_sub_feature_split_summary # disease_subtype

full_model_metric_summary # disease
disease_sub_metric_summary # disease_subtype

full_model_mean_importance # disease
disease_sub_mean_importance # disease_subtype


### disease has very high sensitivity, but very low specificity
### disease_subtype has lower overall performance, but performance across classes is better balanced (doesn't over predict disease as much)
### age_imputed has a very high Gini (appears in every tree), but decreases model accuracy
### location has a modes impact on accuracy
### gender has a very low Gini and decreases model accuracy


############################################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CROSS-VALIDATION + 50 REPEATS - DISEASE_SUBTYPE - MICROBIOME ONLY   ###
############################################################################################################

# using only microbiome features to train the model

# set seed
set.seed(1234)

# column names for microbiome features
micro_feat_cols <- setdiff(colnames(ibd_merge), c("disease", "disease_subtype", "age_imputed", "gender", "location"))

# column names for microbiome features + relevant metadata (full predictor set)
all_feat_cols <- micro_feat_cols  ### only include microbiome features in this model

# create lists to store metrics
feature_importances <- list() # list to store feature importances
performance_metrics <- list() # list to store performance metrics
feature_frequencies <- list() # list to store feature selection frequencies

# repeat cross-validation 50 times
for (r in 1:50) {
  cat("Repeat:", r, "\n")
  
  # create 5-folds for cross-validation (stratified on disease_subtype)
  folds <- createFolds(ibd_merge$disease_subtype, k = 5, list = TRUE)
  
  # loop through the folds
  for (f in 1:5) {
    
    # splits the dataset into training and testing sets for the current fold
    test_idx <- folds[[f]] # test indices for the f-th fold
    train_data <- ibd_merge[-test_idx, ] # training data (all rows not in fold f)
    test_data  <- ibd_merge[test_idx, ] # testing data (fold f)
    
    # train random forest model
    # x = all data in data.frame subset by all_feat_cols (predictor values)
    # y = target variable as factor
    rf_model <- randomForest(x = train_data[, all_feat_cols], 
                             y = as.factor(train_data$disease_subtype), 
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
    cm <- confusionMatrix(predictions, as.factor(test_data$disease_subtype))
    
    # store with repeat (r) and fold (f) index
    # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
    key <- paste0("Repeat_", r, "_Fold_", f)
    feature_frequencies[[key]] <- as.data.frame(split_counts) # store feature freqeuncies
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


### calculate performance statistics (multi-class)
# each confusion matrix now returns a matrix with cm$byClass
# get macro-averaged metrics across the 3 classes (healthy, CD, UC)
# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()

# extract metrics from the stored confusion matrices (50 repeats x 5 folds = 250 values)
for (cm in performance_metrics) {
  if (is.matrix(cm$byClass)) {
    balanced_accuracy <- c(balanced_accuracy, mean(cm$byClass[,"Balanced Accuracy"], na.rm = TRUE))
    f1_score <- c(f1_score, mean(cm$byClass[,"F1"], na.rm = TRUE))
    sensitivity <- c(sensitivity, mean(cm$byClass[,"Sensitivity"], na.rm = TRUE))
    specificity <- c(specificity, mean(cm$byClass[,"Specificity"], na.rm = TRUE))
  }
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
# mean_MeanDecreaseAccuracy: overall importance of feature on model accuracy
# mean_MeanDecreaseGini: frequency and usefulness in splitting (how much a feature reduces impurity when used to split the decision trees)
# Gini is sensitive to splits, NOT predictive value
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_CD = mean(CD, na.rm = TRUE),
            mean_UC = mean(UC, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseGini))
head(mean_importance, 10)


# same data from full model for comparison
microbiome_only_feature_split_summary <- feature_split_summary
microbiome_only_metric_summary <- metric_summary
microbiome_only_mean_importance <- mean_importance


### comparison of disease_subtype full model versus microbiome only model
disease_sub_feature_split_summary # full model
microbiome_only_feature_split_summary # microbiome only

disease_sub_metric_summary # full model
microbiome_only_metric_summary # microbiome only

disease_sub_mean_importance # full model
microbiome_only_mean_importance # microbiome only


### model performance metrics are essentially identical
   # microbiome features alone are capturing most of the discriminative signal
### top microbiome features are stable across both models


#########################################################################################################
###   OVERALL RANDOM FOREST - 5-FOLD CV + 50 REPEATS - DISEASE_SUBTYPE - MICROBIOME ONLY - STRATIFIED ###
#########################################################################################################

# stratified performance analysis (sex and age)
# evaluate whether model performance and feature importance varies across subgroups

# split age into two groups by the median
median_age <- median(ibd_merge$age_imputed, na.rm = TRUE)
ibd_merge$age_group <- ifelse(ibd_merge$age_imputed <= median_age, "younger", "older")

# set seed
set.seed(1234)

# column names for microbiome features
micro_feat_cols <- setdiff(colnames(ibd_merge), c("disease", "disease_subtype", "age_imputed", "gender", "location"))

# column names for microbiome features + relevant metadata (full predictor set)
all_feat_cols <- micro_feat_cols  ### only include microbiome features in this model

### set stratification variable
strat_var <- "gender" # "gender" or "age_group"
# get the levels to compare
subgroups <- unique(ibd_merge[[strat_var]])

# create lists to store stratified metrics
stratified_metrics <- list()
stratified_importance <- list()

# loop over subgroups
for (group in subgroups){
  cat("Analyzing subgroup:", group, "\n")
  
  # subset data by the stratificiation variable
  data_sub <- ibd_merge[ibd_merge[[strat_var]] == group, ]
  
  # create lists to store metrics
  feature_importances <- list() # list to store feature importances
  performance_metrics <- list() # list to store performance metrics
  
  # repeat cross-validation 50 times
  for (r in 1:50) {
    cat("Repeat:", r, "\n")
    
    # create 5-folds for cross-validation (stratified on disease_subtype) using data subset on the stratification variable
    folds <- createFolds(data_sub$disease_subtype, k = 5, list = TRUE)
    
    # loop through the folds
    for (f in 1:5) {
      
      # splits the dataset into training and testing sets for the current fold (data already subset on stratification variable)
      test_idx <- folds[[f]] # test indices for the f-th fold
      train_data <- data_sub[-test_idx, ] # training data (all rows not in fold f)
      test_data  <- data_sub[test_idx, ] # testing data (fold f)
      
      # train random forest model
      # x = all data in data.frame subset by all_feat_cols (predictor values)
      # y = target variable as factor
      rf_model <- randomForest(x = train_data[, all_feat_cols], 
                               y = as.factor(train_data$disease_subtype), 
                               ntree = 500, importance = TRUE) 
      
      # evaluate on test set
      predictions <- predict(rf_model, newdata = test_data[, all_feat_cols])
      
      # generate confusion matrix
      cm <- confusionMatrix(predictions, as.factor(test_data$disease_subtype))
      
      # store with repeat (r) and fold (f) index
      # performance_metrics and feature_importances will be lists of 250 elements (50 repeats x 5 folds)
      key <- paste0("Repeat_", r, "_Fold_", f)
      performance_metrics[[key]] <- cm # store performance metrics
      feature_importances[[key]] <- importance(rf_model)  # store feature importances
    }
  }


### calculate performance statistics

# create vectors to store metrics
balanced_accuracy <- numeric()
f1_score <- numeric()
sensitivity <- numeric()
specificity <- numeric()

# extract metrics from the stored confusion matrices (50 repeats x 5 folds = 250 values)
for (cm in performance_metrics) {
  if (is.matrix(cm$byClass)) {
    balanced_accuracy <- c(balanced_accuracy, mean(cm$byClass[,"Balanced Accuracy"], na.rm = TRUE))
    f1_score <- c(f1_score, mean(cm$byClass[,"F1"], na.rm = TRUE))
    sensitivity <- c(sensitivity, mean(cm$byClass[,"Sensitivity"], na.rm = TRUE))
    specificity <- c(specificity, mean(cm$byClass[,"Specificity"], na.rm = TRUE))
  }
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
# add group name to summary metrics info
metric_summary$Group <- group

# store stratified performance metrics (each group has a list of two elements: summary and full metrics)
stratified_metrics[[group]] <- list(summary = metric_summary,
                                    full_metrics = performance_metrics)


### calculate feature importances
# combine all feature_importances data.frames into one data.frame
all_features_importances <- do.call(rbind, lapply(names(feature_importances), function(name) {
  df <- as.data.frame(feature_importances[[name]])
  df$Feature <- rownames(df)
  df$Repeat_Fold <- name
  return(df)
}))

# group importance metrics by feature and sort by overall importance 
# mean_MeanDecreaseAccuracy: overall importance of feature on model accuracy
# mean_MeanDecreaseGini: frequency and usefulness in splitting (how much a feature reduces impurity when used to split the decision trees)
# Gini is sensitive to splits, NOT predictive value
mean_importance <- all_features_importances %>%
  group_by(Feature) %>%
  summarise(mean_healthy = mean(healthy, na.rm = TRUE),
            mean_CD = mean(CD, na.rm = TRUE),
            mean_UC = mean(UC, na.rm = TRUE),
            mean_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            mean_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(mean_MeanDecreaseGini))
head(mean_importance, 10)

stratified_importance[[group]] <- mean_importance
}

### does the model perform statistically differently across the subgroups
# function to extract metrics from confusion matrices
extract_metric <- function(metric_name, cm_list) {
  sapply(cm_list, function(cm) {
    if (!is.null(cm) && is.matrix(cm$byClass)) {
      mean(cm$byClass[, metric_name], na.rm = TRUE)
    } else {
      NA
    }
  })
}

# function to extract and average metric across folds per repeat (per repeat averages)
get_repeat_averages <- function(group, metric) {
  vals <- extract_metric(metric, stratified_metrics[[group]]$full_metrics)
  rowMeans(matrix(vals, nrow = 50, byrow = TRUE), na.rm = TRUE)
}

# MALE VS FEMALE
# balanced accuracy
bal_acc_male_avg <- get_repeat_averages("male", "Balanced Accuracy")
bal_acc_female_avg <- get_repeat_averages("female", "Balanced Accuracy")
wilcox.test(bal_acc_male_avg, bal_acc_female_avg) # p-value = 1.47e-14

perf_df <- data.frame(BalancedAccuracy = c(bal_acc_male_avg, bal_acc_female_avg),
                      Group = rep(c("Male", "Female"), each = 50))
ggplot(perf_df, aes(x = Group, y = BalancedAccuracy, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink")) +
  labs(title = "Balanced accuracy by gender", y = "Balanced accuracy", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# f1 score
f1_male_avg <- get_repeat_averages("male", "F1")
f1_female_avg <- get_repeat_averages("female", "F1")
wilcox.test(f1_male_avg, f1_female_avg) # p-value = 0.759

perf_df <- data.frame(F1_score = c(f1_male_avg, f1_female_avg),
                      Group = rep(c("Male", "Female"), each = 50))
ggplot(perf_df, aes(x = Group, y = F1_score, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink")) +
  labs(title = "F1 score by gender", y = "F1 score", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# sensitivity
sens_male_avg <- get_repeat_averages("male", "Sensitivity")
sens_female_avg <- get_repeat_averages("female", "Sensitivity")
wilcox.test(sens_male_avg, sens_female_avg) # p-value = 2.509e-13

perf_df <- data.frame(Sens = c(sens_male_avg, sens_female_avg),
                      Group = rep(c("Male", "Female"), each = 50))
ggplot(perf_df, aes(x = Group, y = Sens, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink")) +
  labs(title = "Sensitivity by gender", y = "Sensitivity", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# specificity
spec_male_avg <- get_repeat_averages("male", "Specificity")
spec_female_avg <- get_repeat_averages("female", "Specificity")
wilcox.test(spec_male_avg, spec_female_avg) # p-value = 5.037e-16

perf_df <- data.frame(Spec = c(spec_male_avg, spec_female_avg),
                      Group = rep(c("Male", "Female"), each = 50))
ggplot(perf_df, aes(x = Group, y = Spec, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Male" = "skyblue", "Female" = "pink")) +
  labs(title = "Specificity by gender", y = "Specificity", x = NULL) +
  theme_minimal() + theme(legend.position = "none")


# YOUNGER VS OLDER
# balanced accuracy
bal_acc_younger_avg <- get_repeat_averages("younger", "Balanced Accuracy")
bal_acc_older_avg <- get_repeat_averages("older", "Balanced Accuracy")
wilcox.test(bal_acc_younger_avg, bal_acc_older_avg) # p-value = 0.003999

perf_df <- data.frame(BalancedAccuracy = c(bal_acc_younger_avg, bal_acc_older_avg),
                      Group = rep(c("Younger", "Older"), each = 50))
ggplot(perf_df, aes(x = Group, y = BalancedAccuracy, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Younger" = "skyblue", "Older" = "pink")) +
  labs(title = "Balanced accuracy by age", y = "Balanced accuracy", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# f1 score
f1_younger_avg <- get_repeat_averages("younger", "F1")
f1_older_avg <- get_repeat_averages("older", "F1")
wilcox.test(f1_younger_avg, f1_older_avg) # p-value < 2.2e-16

perf_df <- data.frame(F1 = c(f1_younger_avg, f1_older_avg),
                      Group = rep(c("Younger", "Older"), each = 50))
ggplot(perf_df, aes(x = Group, y = F1, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Younger" = "skyblue", "Older" = "pink")) +
  labs(title = "F1 score by age", y = "F1 score", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# sensitivity
sens_younger_avg <- get_repeat_averages("younger", "Sensitivity")
sens_older_avg <- get_repeat_averages("older", "Sensitivity")
wilcox.test(sens_younger_avg, sens_older_avg) # p-value = 0.08532

perf_df <- data.frame(Sens = c(sens_younger_avg, sens_older_avg),
                      Group = rep(c("Younger", "Older"), each = 50))
ggplot(perf_df, aes(x = Group, y = Sens, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Younger" = "skyblue", "Older" = "pink")) +
  labs(title = "Sensitivity by age", y = "Sensitivity score", x = NULL) +
  theme_minimal() + theme(legend.position = "none")

# specificity
spec_younger_avg <- get_repeat_averages("younger", "Specificity")
spec_older_avg <- get_repeat_averages("older", "Specificity")
wilcox.test(spec_younger_avg, spec_older_avg) # p-value = 1.866e-06

perf_df <- data.frame(Spec = c(spec_younger_avg, spec_older_avg),
                      Group = rep(c("Younger", "Older"), each = 50))
ggplot(perf_df, aes(x = Group, y = Spec, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("Younger" = "skyblue", "Older" = "pink")) +
  labs(title = "Specificity by age", y = "Specificity score", x = NULL) +
  theme_minimal() + theme(legend.position = "none")


### are different features predictive in different subgroups

# MALE VS FEMALE
### gender_stratified_importance <- bind_rows(stratified_importance, .id = "Group")
gender_stratified_importance <- gender_stratified_importance %>%
  arrange(desc(mean_MeanDecreaseAccuracy))

# extract top 50 features in males
top_features_male <- stratified_importance[["male"]] %>%
  arrange(desc(mean_MeanDecreaseAccuracy)) %>%
  slice_head(n = 50) %>%
  pull(Feature)

# extract top 50 features in females
top_features_female <- stratified_importance[["female"]] %>%
  arrange(desc(mean_MeanDecreaseAccuracy)) %>%
  slice_head(n = 50) %>%
  pull(Feature)

# shared and unique features
shared_features <- intersect(top_features_male, top_features_female)
unique_to_male <- setdiff(top_features_male, top_features_female)
unique_to_female <- setdiff(top_features_female, top_features_male)

# Venn diagram of shared and unique features
venn.plot <- venn.diagram(
  x = list(Male = top_features_male, Female = top_features_female),
  filename = NULL, fill = c("skyblue", "pink"), alpha = 0.5,
  cex = 1.5, cat.cex = 1.2, main = paste("Top 50 features for males and females"))
grid.draw(venn.plot)


# YOUNGER VS OLDER
### age_stratified_importance <- bind_rows(stratified_importance, .id = "Group")
age_stratified_importance <- age_stratified_importance %>%
  arrange(desc(mean_MeanDecreaseAccuracy))

# extract top 50 features in males
top_features_younger <- stratified_importance[["younger"]] %>%
  arrange(desc(mean_MeanDecreaseAccuracy)) %>%
  slice_head(n = 50) %>%
  pull(Feature)

# extract top 50 features in females
top_features_older <- stratified_importance[["older"]] %>%
  arrange(desc(mean_MeanDecreaseAccuracy)) %>%
  slice_head(n = 50) %>%
  pull(Feature)

# shared and unique features
shared_features <- intersect(top_features_younger, top_features_older)
unique_to_younger <- setdiff(top_features_younger, top_features_older)
unique_to_older <- setdiff(top_features_older, top_features_younger)

# Venn diagram of shared and unique features
venn.plot <- venn.diagram(
  x = list(Younger = top_features_younger, Older = top_features_older),
  filename = NULL, fill = c("skyblue", "pink"), alpha = 0.5,
  cex = 1.5, cat.cex = 1.2, main = paste("Top 50 features for younger and older"))
grid.draw(venn.plot)


# IN ALL
### top features found in all stratified models
in_age <- intersect(top_features_younger, top_features_older)
in_gender <- intersect(top_features_male, top_features_female)
in_all <- intersect(in_age, in_gender)
in_all

ggplot(ibd_merge, aes(x = disease_subtype, y = Bacteroides_stercoris, fill = disease_subtype)) +
  geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("healthy" = "skyblue", "CD" = "pink", "UC" = "plum")) +
  labs(title = "Relative abundance of Bacteroides stercoris", 
       y = "Relative abundance", x = NULL) + 
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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] VennDiagram_1.7.3               futile.logger_1.4.3             pROC_1.18.5                    
# [4] cowplot_1.1.3                   e1071_1.7-16                    caret_7.0-1                    
# [7] lattice_0.22-7                  mice_3.18.0                     vegan_2.7-1                    
# [10] permute_0.9-7                   randomForest_4.7-1.2            lubridate_1.9.4                
# [13] forcats_1.0.0                   stringr_1.5.1                   dplyr_1.1.4                    
# [16] purrr_1.0.4                     readr_2.1.5                     tidyr_1.3.1                    
# [19] tibble_3.3.0                    tidyverse_2.0.0                 ggplot2_3.5.2                  
# [22] curatedMetagenomicData_3.16.1   TreeSummarizedExperiment_2.16.1 Biostrings_2.76.0              
# [25] XVector_0.48.0                  SingleCellExperiment_1.30.1     SummarizedExperiment_1.38.1    
# [28] Biobase_2.68.0                  GenomicRanges_1.60.0            GenomeInfoDb_1.44.0            
# [31] IRanges_2.42.0                  S4Vectors_0.46.0                BiocGenerics_0.54.0            
# [34] generics_0.1.4                  MatrixGenerics_1.20.0           matrixStats_1.5.0              
# 
# loaded via a namespace (and not attached):
#   [1] ggtext_0.1.2                fs_1.6.6                    DirichletMultinomial_1.50.0
# [4] httr_1.4.7                  RColorBrewer_1.1-3          tools_4.5.0                
# [7] backports_1.5.0             utf8_1.2.6                  R6_2.6.1                   
# [10] lazyeval_0.2.2              mgcv_1.9-3                  jomo_2.7-6                 
# [13] withr_3.0.2                 gridExtra_2.3               cli_3.6.5                  
# [16] textshaping_1.0.1           formatR_1.14                sandwich_3.1-1             
# [19] labeling_0.4.3              slam_0.1-55                 mvtnorm_1.3-3              
# [22] proxy_0.4-27                systemfonts_1.2.3           yulab.utils_0.2.0          
# [25] dichromat_2.0-0.1           scater_1.35.0               decontam_1.28.0            
# [28] parallelly_1.45.0           readxl_1.4.5                fillpattern_1.0.2          
# [31] rstudioapi_0.17.1           RSQLite_2.4.1               shape_1.4.6.1              
# [34] rbiom_2.2.0                 Matrix_1.7-3                ggbeeswarm_0.7.2           
# [37] DECIPHER_3.4.0              abind_1.4-8                 lifecycle_1.0.4            
# [40] multcomp_1.4-28             yaml_2.3.10                 recipes_1.3.1              
# [43] SparseArray_1.8.0           BiocFileCache_2.16.0        blob_1.2.4                 
# [46] ExperimentHub_2.16.0        crayon_1.5.3                mitml_0.4-5                
# [49] beachmat_2.24.0             KEGGREST_1.48.0             pillar_1.10.2              
# [52] boot_1.3-31                 estimability_1.5.1          future.apply_1.20.0        
# [55] codetools_0.2-20            pan_1.9                     glue_1.8.0                 
# [58] data.table_1.17.4           MultiAssayExperiment_1.34.0 vctrs_0.6.5                
# [61] png_0.1-8                   treeio_1.32.0               Rdpack_2.6.4               
# [64] cellranger_1.1.0            gtable_0.3.6                cachem_1.1.0               
# [67] gower_1.0.2                 rbibutils_2.3               S4Arrays_1.8.1             
# [70] mime_0.13                   prodlim_2025.04.28          coda_0.19-4.1              
# [73] reformulas_0.4.1            survival_3.8-3              timeDate_4041.110          
# [76] iterators_1.0.14            hardhat_1.4.1               lava_1.8.1                 
# [79] bluster_1.18.0              TH.data_1.1-3               ipred_0.9-15               
# [82] nlme_3.1-168                bit64_4.6.0-1               filelock_1.0.3             
# [85] irlba_2.3.5.1               vipor_0.4.7                 rpart_4.1.24               
# [88] DBI_1.2.3                   nnet_7.3-20                 tidyselect_1.2.1           
# [91] emmeans_1.11.1              bit_4.6.0                   compiler_4.5.0             
# [94] curl_6.3.0                  glmnet_4.1-9                BiocNeighbors_2.2.0        
# [97] xml2_1.3.8                  DelayedArray_0.34.1         scales_1.4.0               
# [100] rappdirs_0.3.3              digest_0.6.37               minqa_1.2.8                
# [103] pkgconfig_2.0.3             lme4_1.1-37                 sparseMatrixStats_1.20.0   
# [106] dbplyr_2.5.0                fastmap_1.2.0               rlang_1.1.6                
# [109] UCSC.utils_1.4.0            DelayedMatrixStats_1.30.0   farver_2.1.2               
# [112] zoo_1.8-14                  jsonlite_2.0.0              BiocParallel_1.42.1        
# [115] ModelMetrics_1.2.2.2        BiocSingular_1.24.0         magrittr_2.0.3             
# [118] scuttle_1.18.0              GenomeInfoDbData_1.2.14     patchwork_1.3.0            
# [121] Rcpp_1.0.14                 ape_5.8-1                   ggnewscale_0.5.1           
# [124] viridis_0.6.5               stringi_1.8.7               MASS_7.3-65                
# [127] AnnotationHub_3.16.0        plyr_1.8.9                  parallel_4.5.0             
# [130] listenv_0.9.1               ggrepel_0.9.6               splines_4.5.0              
# [133] gridtext_0.1.5              hms_1.1.3                   igraph_2.1.4               
# [136] reshape2_1.4.4              ScaledMatrix_1.16.0         pkgload_1.4.0              
# [139] futile.options_1.0.1        BiocVersion_3.21.1          lambda.r_1.2.4             
# [142] BiocManager_1.30.26         nloptr_2.2.1                tzdb_0.5.0                 
# [145] foreach_1.5.2               future_1.58.0               BiocBaseUtils_1.10.0       
# [148] rsvd_1.0.5                  broom_1.0.8                 xtable_1.8-4               
# [151] tidytree_0.4.6              viridisLite_0.4.2           class_7.3-23               
# [154] ragg_1.4.0                  memoise_2.0.1               beeswarm_0.4.0             
# [157] AnnotationDbi_1.70.0        cluster_2.1.8.1             timechange_0.3.0           
# [160] globals_0.18.0              mia_1.16.0

