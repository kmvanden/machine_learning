# SIAMCAT: Overall Workflow, Holdout Testing and ML Pitfalls

# https://siamcat.embl.de/articles/SIAMCAT_vignette.html
# https://siamcat.embl.de/articles/SIAMCAT_holdout.html
# https://siamcat.embl.de/articles/SIAMCAT_ml_pitfalls.html

### SIAMCAT input
# features: matrix, data.frame or otu_table
# features in rows and samples in columns
### features need to be RELATIVE ABUNDANCES

# metadata: matrix or data.frame
# samples in rows and metadata in columns
### rownames in the metadata should MATCH the colnames in the feature matrix

# label: named vector or metadata column
### should match rownames in the metadata and colnames in the feature matrix


################################
### OVERALL SIAMCAT WORKFLOW ###
################################


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(SIAMCAT)
library(ggplot2)
library(tidyverse)

# load data
# microbiome data from colorectal cancer study
data("feat_crc_zeller", package = "SIAMCAT") # feature data
feat.crc.zeller[1:3, 1:3]
dim(feat.crc.zeller) # 1754  141
data("meta_crc_zeller", package = "SIAMCAT") # metadata
head(meta.crc.zeller)

# create label object from Group column in metadata
label.crc.zeller <- create.label(meta = meta.crc.zeller, label = "Group", case = "CRC")


# CREATE SIAMCAT OBJECT
sc.obj <- siamcat(feat = feat.crc.zeller, label = label.crc.zeller, meta = meta.crc.zeller)
show(sc.obj)
# siamcat-class object
# label()                Label object:         88 CTR and 53 CRC samples
# 
# contains phyloseq-class experiment-level object @phyloseq:
# phyloseq@otu_table()   OTU Table:            [ 1754 taxa and 141 samples ]
# phyloseq@sam_data()    Sample Data:          [ 141 samples by 6 sample variables ]


# perform unsupervised feature selection 
# remove features whose maximum abundance is never above 0.001
sc.obj <- filter.features(sc.obj, filter.method = "abundance", cutoff = 0.001)


### ASSOCIATION TESTING
# computes associations between microbial species  and the disease (condition)
# uses a non-parametric Wilcoxon test and different effect sizes for the association (AUC of fold change)
sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = "fc", panels = c("fc", "prevalence", "auroc"))


### CONFOUNDER TESTING
# check for other variables that influence the microbiome other than disease (condition)
# tests the associated metadata variables for potential confounding influence
  # do the distributions look similar between the case and the controls?
# conditional entropy: primarily serves to remove nonsensical variables from subsequent checks
# variance explained by [factor]: plots the variance explained by the disease (condition) against the variance explained by [factor] for each individual feature
# variables with many features in the upper left-hand corner might be confounding the disease associations

check.confounders(sc.obj, meta.in = NULL, feature.type = "filtered", fn.plot = "confounder_plots.pdf")


######################### FROM SIAMCAT SOURCE CODE #########################

######################### GET factorize.metadata()

factorize.metadata <- function(meta, verbose) {
  if ('BMI' %in% toupper(colnames(meta))) {
    idx <- match('BMI', toupper(colnames(meta)))
    meta[,idx] <- factorize.bmi(meta[,idx])
  }
  
  factorized <- as.data.frame(lapply(meta, FUN=function(x) {
    if (is.numeric(x) & (length(unique(x)) > 5)) {
      quart <- quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE)
      temp <- cut(x, unique(quart), include.lowest = TRUE)
      return(factor(temp, labels = seq_along(levels(temp))))
    } else {
      return(as.factor(x))
    }
  }))
  rownames(factorized) <- rownames(meta)
  
  n.levels <- vapply(colnames(factorized), FUN=function(x){
    length(levels(factorized[[x]]))
  }, FUN.VALUE = integer(1))
  
  if (any(n.levels > 0.9 * nrow(meta))) {
    remove.meta <- names(which(n.levels > 0.9 * nrow(meta)))
    if (verbose > 1){
      message("++ metadata variables:\n\t",
              paste(remove.meta, collapse = " & "),
              "\n++ have too many levels and have been removed")
    }
    factorized <- factorized[, !(colnames(factorized) %in% remove.meta)]
  }
  
  return(factorized)
}

##################################################

# get label, meta and feat
label <- label(sc.obj)
meta <- meta(sc.obj)
feat <- get.filt_feat.matrix(sc.obj)

# make sure sample names match
all(colnames(feat) == names(label$label)) 

### get variance explained by label
var.label <- vapply(rownames(feat), FUN = function(x){
  x_vals <- rank(feat[x,]) / length(feat[x,])
  ss.tot <- sum((x_vals - mean(x_vals))^2) / length(x_vals)
  ss.o.i <- sum(vapply(unique(label$label), function(s){
    sum((x_vals[label$label == s] - mean(x_vals[label$label == s]))^2)
  }, FUN.VALUE = double(1))) / length(x_vals)
  return(1 - ss.o.i / ss.tot)
}, FUN.VALUE = double(1))

### get variance explained by Age
Age <- meta$Age
Age <- factorize.metadata(data.frame(Age = Age), verbose = 0)$Age
names(Age) <- rownames(meta)

# ensure names stil match
all(names(Age) %in% colnames(feat))
all(colnames(feat) %in% names(Age))

var.age <- vapply(rownames(feat), FUN = function(x){
  x_vals <- rank(feat[x, names(Age)]) / length(Age)
  ss.tot <- sum((x_vals - mean(x_vals))^2) / length(x_vals)
  ss.o.i <- sum(vapply(levels(Age), function(s){
    sum((x_vals[Age == s] - mean(x_vals[Age == s]))^2)
  }, FUN.VALUE = double(1))) / length(x_vals)
  return(1 - ss.o.i / ss.tot)
}, FUN.VALUE = double(1))


# plot size as average log abundance values
r.mean   <- rowMeans(log10(feat + 1e-05))
mean.size <- (r.mean + 5) * 8/5 + 1

# define maximum values for axes
lim <- max(c(var.label, var.age), na.rm = TRUE) 

plot(var.label, var.age, type = "n",
     xlab = "Variance explained by label",
     ylab = "Variance expained by age",
     xlim = c(0, lim), ylim = c(0, lim))
symbols(x = var.label, y = var.age, circles = mean.size, inches = 1/9,
        bg = alpha("darkgrey", 0.4), fg = alpha("black", 0.7), add = TRUE)
abline(0,1, lty = 3, col = "black")


sort_var <- sort(var.label, decreasing = TRUE)
sort_var[1:5]


###################### ABOVE FROM SIAMCAT SOURCE CODE ######################

### MODEL BUILDING
# data normalization
# normalization methods available: "rank.unit", "rank.std", "log.std", "log.unit", "log.clr", "std", "pass"
sc.obj <- normalize.features(sc.obj, norm.method = "log.unit",
                             norm.param = list(log.n0 = 1e-06, n.p = 2,norm.margin = 1))

# split the data into cross-validation folds
# num.folds = number of folds | num.resample = number of resampling
# stratify = boolean, stratify the splits so that an equal proportion of classes are present in each fold
# inseparable = string, name of metadata variable to be inseparable (e.g., replicates or times series)
sc.obj <-  create.data.split(sc.obj, num.folds = 5, num.resample = 2)

# model training
# method  = "lasso", "enet", "ridge", "lasso_ll", "ridge_ll", "randomForest"
# param.set = other parameters to run depending on the model used
sc.obj <- train.model(sc.obj, method = "lasso")
model_type(sc.obj) # lasso"

# access models
models <- models(sc.obj)
models[[1]]$model

# predictions: apply the model to the test folds (left out data)
# predictions are saved in the pred_matrix slot of the SIAMCAT object
sc.obj <- make.predictions(sc.obj)
pred_matrix <- pred_matrix(sc.obj)
head(pred_matrix)

# model evaluation and interpretation
# calculate AUROC and the precision recall curve for each resampled CV run
sc.obj <-  evaluate.predictions(sc.obj)

# evaluation plots (AUROC and PR curves)
model.evaluation.plot(sc.obj) 

# for the top selected features 
    # barplot of model weights and how robust they are (in how many models in the repeated CV that feature was selected into the model)
    # heatmap of z-scores or fold changes (distribution of metadata is shown below heatmap)
    # box plot showing the proportion of the the model weight that is explained by the selected features
model.interpretation.plot(sc.obj, fn.plot = "interpretation.pdf",
                          consens.thres = 0.5, limits = c(-3, 3), 
                          heatmap.type = "zscore") # interpretation plot


##########################################
### SIAMCAT WORKFLOW - HOLDOUT TESTING ###
##########################################

# train a model on one dataset and then apply it to another similarly processed holdout dataset
# two datasets from two studies on colorectal cancer (France and China)
# both datasets were profiles with the same taxonomic profiling tool (mOTUs)

# load libraries
library(SIAMCAT)

# load data
data.loc <- "https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/"

# study from France
fn.meta.fr  <- paste0(data.loc, "meta_Zeller.tsv")
fn.feat.fr  <- paste0(data.loc, "specI_Zeller.tsv")

# study from China
fn.meta.cn  <- paste0(data.loc, "meta_Yu.tsv")
fn.feat.cn  <- paste0(data.loc, "specI_Yu.tsv")

### build a SIAMCAT object using the data from the study in France
# read in feature data
feat.fr  <- read.table(fn.feat.fr, sep = "\t", quote = "", check.names = FALSE, 
                       stringsAsFactors = FALSE)
head(feat.fr)

# change the data.frame from counts to relative abundances
feat.fr.rel <- prop.table(as.matrix(feat.fr), 2)
head(feat.fr.rel)

# read in metadata
meta.fr  <- read.table(fn.meta.fr, sep = "\t", quote = "",
                       check.names = FALSE, stringsAsFactors = FALSE)

# create a SIAMCAT object for the study in France
siamcat.fr <- siamcat(feat = feat.fr.rel, meta = meta.fr, label = "Group", case = "CRC")


### build a SIAMCAT object using the data from the study in China
# read in feature data
feat.cn  <- read.table(fn.feat.cn, sep = "\t", quote = "", check.names = FALSE, 
                       stringsAsFactors = FALSE)
head(feat.cn)

# change the data.frame from counts to relative abundances
feat.cn.rel <- prop.table(as.matrix(feat.cn), 2)
head(feat.cn.rel)

# read in metadata
meta.cn  <- read.table(fn.meta.cn, sep = "\t", quote = "",
                       check.names = FALSE, stringsAsFactors = FALSE)

# create SIAMCAT object for the study in China
siamcat.cn <- siamcat(feat = feat.cn.rel, meta = meta.cn, label = "Group", case = "CRC")


### build the model using the French dataset
# filter data by abundance 
siamcat.fr <- filter.features(siamcat.fr, filter.method = "abundance",
                              cutoff = 0.001, rm.unmapped = TRUE, verbose = 2)
# normalize data
siamcat.fr <- normalize.features(siamcat.fr, norm.method = "log.std",
                                 norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1), 
                                 verbose = 2)
# split the data
siamcat.fr <-  create.data.split(siamcat.fr, num.folds = 5, num.resample = 2)
# train the data using a LASSO model
siamcat.fr <- train.model(siamcat.fr, method = "lasso")

# make and evaluate predictions for each CV fold
siamcat.fr <- make.predictions(siamcat.fr)
siamcat.fr <-  evaluate.predictions(siamcat.fr)


### APPLY MODEL TO HOLDOUT DATASET
# # apply the model built for the French dataset to the Chinese dataset

# normalize the dataset from China using hte same parameters that were used for the dataset from France 
siamcat.cn <- normalize.features(siamcat.cn, norm.param = norm_params(siamcat.fr),
                                 feature.type = "original", verbose = 2)

# apply the trained model to the holdout dataset
siamcat.cn <- make.predictions(siamcat = siamcat.fr, siamcat.holdout = siamcat.cn,
                               normalize.holdout = FALSE)

# evaluate the predictions
siamcat.cn <- evaluate.predictions(siamcat.cn)

# compare the performance of the classifier on the original and on the holdout dataset
# supply the two SIAMCAT objects and the model evaluation is plotted on the same plot for both
model.evaluation.plot("FR-CRC" = siamcat.fr, "CN-CRC" = siamcat.cn, colours = c("dimgrey", "orange"))


####################################################
### SIAMCAT WORKFLOW - MACHINE LEARNING PITFALLS ###
####################################################

# two pitfalls in machine learning that can lead to overly optimistic performance estimates
  # 1. supervised feature selection (label information is taken into account before CV)
      # features are selected if they are associated with the label
      # CORRECT: nest the selection step into the CV procedure --> feature selection is performed for each training fold separately
  # 2. naive splitting of dependent data

# load libraries
library(tidyverse)
library(SIAMCAT)
library(curatedMetagenomicData)
library(dplyr)

# load data
# two datasets of colorectal cancer 

### THOMAS
# determine the current name of the dataset in curatedMetagenomicData
datasets <- curatedMetagenomicData(pattern = ".", dryrun = TRUE)
grep("ThomasAM_2018", datasets, value = TRUE)

feat.t <- curatedMetagenomicData("2021-03-31.ThomasAM_2018a.relative_abundance", dryrun = FALSE)
tse <- feat.t[[1]] # extract the single TreeSumarizedExoeriment from the list
feat.t <- assay(tse, "relative_abundance") # extract the abundance matrix

# clean up metaphlan profiles to contain only species-level abundances
feat.t <- feat.t[grep(x = rownames(feat.t), pattern = "s__"),]
feat.t <- feat.t[grep(x = rownames(feat.t),pattern = "t__", invert = TRUE),]
feat.t <- t(t(feat.t)/100)

# get metadata
meta.t <- as.data.frame(colData(tse))
# filter metadata
meta.t <- meta.t %>%
  filter(study_condition %in% c("control", "CRC"))
rownames(meta.t) <- meta.t$subject_id # should be subject_id not sampleID

# edit colnames of the feature data to match that of the rownames of the column data
colnames(feat.t) <- sub("^[^_]+_([^_]+)_[^_]+$", "\\1", colnames(feat.t))
rownames(meta.t) %in% colnames(feat.t) # check that they match


### ZELLER
datasets <- curatedMetagenomicData(pattern = ".", dryrun = TRUE)
grep("ZellerG_2014", datasets, value = TRUE)

feat.z <- curatedMetagenomicData("2021-03-31.ZellerG_2014.relative_abundance", dryrun = FALSE)
zse <- feat.z[[1]] # extract the single TreeSumarizedExoeriment from the list
feat.z <- assay(zse, "relative_abundance") # extract the abundance matrix

# clean up metaphlan profiles to contain only species-level abundances
feat.z <- feat.z[grep(x = rownames(feat.z), pattern = "s__"),]
feat.z <- feat.z[grep(x = rownames(feat.z),pattern = "t__", invert = TRUE),]
feat.z <- t(t(feat.z)/100) 

# get metadata
meta.z <- as.data.frame(colData(zse))
# filter metadata
meta.z <- meta.z %>%
  filter(study_condition %in% c("control", "CRC"))
rownames(meta.z) <- meta.z$subject_id # should be subject_id not sampleID
colnames(feat.z) # not the subject_id

names <- as.data.frame(colData(zse)) # get the subject_id
colnames(feat.z) <- names$subject_id # change colnames of feature table to subject_id
rownames(meta.z) %in% colnames(feat.z) # check that they match


# in order to be able to use the datasets as training and external test sets in SIAMCAT, the features in both datasets need to overlap completely
species.union <- union(rownames(feat.t), rownames(feat.z))
# add Zeller-only species to the Thomas matrix
add.species <- setdiff(species.union, rownames(feat.t))
feat.t <- rbind(feat.t, matrix(0, nrow = length(add.species), ncol = ncol(feat.t),
                       dimnames = list(add.species, colnames(feat.t))))
# add Thomas-only species to the Zeller matrix
add.species <- setdiff(species.union, rownames(feat.z))
feat.z <- rbind(feat.z, matrix(0, nrow = length(add.species), ncol = ncol(feat.z),
                       dimnames = list(add.species, colnames(feat.z))))


### MODEL TRAINING (LASSO) ####################################

# chose three different feature selection cutoff values
fs.cutoff <- c(20, 100, 250)
# create a tibble to hold the results
auroc.all <- tibble(cutoff = character(0), type = character(0), 
                    study.test = character(0), AUC = double(0))

### TRAIN MODEL WITHOUT FEATURE SELECTION
# train model without feature selection and add to results model twice (correct and incorrect) for easier plotting
sc.obj.t <- siamcat(feat = feat.t, meta = meta.t, label = "study_condition", case = "CRC")
sc.obj.t <- filter.features(sc.obj.t, filter.method = "prevalence", cutoff = 0.01)
sc.obj.t <- normalize.features(sc.obj.t, norm.method = "log.std", 
                               norm.param = list(log.n0 = 1e-05, sd.min.q = 0))
sc.obj.t <- create.data.split(sc.obj.t, num.folds = 10, num.resample = 10)
sc.obj.t <- train.model(sc.obj.t, method = "lasso")
sc.obj.t <- make.predictions(sc.obj.t)
sc.obj.t <- evaluate.predictions(sc.obj.t)

auroc.all <- auroc.all %>% 
  add_row(cutoff = "full", type = "correct", study.test = "Thomas_2018", 
          AUC = as.numeric(sc.obj.t@eval_data$auroc)) %>% 
  add_row(cutoff = "full", type = "incorrect", study.test = "Thomas_2018", 
          AUC = as.numeric(sc.obj.t@eval_data$auroc)) 


### apply the model to the external dataset and record the generalization
sc.obj.z <- siamcat(feat = feat.z, meta = meta.z, label = "study_condition", case = "CRC")
sc.obj.z <- make.predictions(sc.obj.t, sc.obj.z)
sc.obj.z <- evaluate.predictions(sc.obj.z)
auroc.all <- auroc.all %>% 
  add_row(cutoff = "full", type = "correct", 
          study.test = "Zeller_2014", 
          AUC = as.numeric(sc.obj.z@eval_data$auroc)) %>% 
  add_row(cutoff = "full", type = "incorrect", 
          study.test = "Zeller_2014", 
          AUC = as.numeric(sc.obj.z@eval_data$auroc)) 


### INCORRECT WORKFLOW - TRAIN WITH SUPERVISED FEATURE SECTION 
# test the features for differential abundance using the complete dataset and then chose the top associated features 

sc.obj.t <- check.associations(sc.obj.t)
# association.plot(sc.obj.t, plot.type = "quantile.rect", fn.plot = "./association_plot.pdf")
mat.assoc <- associations(sc.obj.t)
mat.assoc$species <- rownames(mat.assoc)
mat.assoc <- mat.assoc %>% as_tibble() %>% arrange(p.val) # sort by p-value


# based on the p-values from check.associations, x number of features are chosen to train the model
# Error in h(simpleError(msg, call)) : 
#   error in evaluating the argument 'x' in selecting a method for function 'slice': Rle of type 'list' is not supported
# tried converting mat.assoc to a data.frame before loop, but that didn't fix the error
# changed the script to use base R indexing rather than slice()

for (x in fs.cutoff){
  # select x number of features based on p-value ranking
  feat.train.red <- feat.t[mat.assoc[seq_len(x), ] %>%
                             pull(species), ]
  sc.obj.t.fs <- siamcat(feat = feat.train.red, meta = meta.t, 
                         label = "study_condition", case = "CRC")
  # normalize the features without filtering
  sc.obj.t.fs <- normalize.features(sc.obj.t.fs, norm.method = "log.std",
                                    norm.param = list(log.n0 = 1e-05, sd.min.q = 0), 
                                    feature.type = "original")
  # take the same cross validation split as before
  data_split(sc.obj.t.fs) <- data_split(sc.obj.t)
  # train
  sc.obj.t.fs <- train.model(sc.obj.t.fs, method = "lasso")
  # make predictions
  sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
  # evaluate predictions and record the result
  sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
  auroc.all <- auroc.all %>% 
    add_row(cutoff = as.character(x), type = "incorrect", 
            study.test = "Thomas_2018",
            AUC = as.numeric(sc.obj.t.fs@eval_data$auroc))
  # apply to the external dataset and record the result
  sc.obj.z <- siamcat(feat = feat.z, meta = meta.z,
                      label = "study_condition", case = "CRC")
  sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
  sc.obj.z <- evaluate.predictions(sc.obj.z)
  auroc.all <- auroc.all %>% 
    add_row(cutoff = as.character(x), type = "incorrect", 
            study.test = "Zeller_2014", 
            AUC = as.numeric(sc.obj.z@eval_data$auroc))
}


### CORRECT WORKFLOW - TRAIN MODEL WITH NESTED FEATURE SELECTION 
# feature selection needs to be nested within the cross-validation procedure
# in SIAMCAT, set perform.fs = TRUE in train.model()
# updated arguments in train.model to: no_features = x, method = "AUC"

for (x in fs.cutoff){
  # train using the original SIAMCAT object 
  # with correct version of feature selection
  sc.obj.t.fs <- train.model(sc.obj.t, method = "lasso", perform.fs = TRUE, 
                             param.fs = list(no_features = x, method = "AUC", direction = "absolute"))
  # make predictions
  sc.obj.t.fs <- make.predictions(sc.obj.t.fs)
  # evaluate predictions and record the result
  sc.obj.t.fs <- evaluate.predictions(sc.obj.t.fs)
  auroc.all <- auroc.all %>% 
    add_row(cutoff = as.character(x), type = "correct", 
            study.test = "Thomas_2018",
            AUC = as.numeric(sc.obj.t.fs@eval_data$auroc))
  # apply to the external dataset and record the result
  sc.obj.z <- siamcat(feat = feat.z, meta = meta.z,
                      label = "study_condition", case = "CRC")
  sc.obj.z <- make.predictions(sc.obj.t.fs, sc.obj.z)
  sc.obj.z <- evaluate.predictions(sc.obj.z)
  auroc.all <- auroc.all %>% 
    add_row(cutoff = as.character(x), type = "correct", 
            study.test = "Zeller_2014", 
            AUC = as.numeric(sc.obj.z@eval_data$auroc))
}


### PLOT THE RESULTS
# plot the resulting performance estimates for the CV and the external validation

auroc.all %>%
  # facetting for plotting
  mutate(split = case_when(study.test == "Thomas_2018"~
                         "Cross validation (Thomas 2018)",
                       TRUE ~ "External validation (Zeller 2014)")) %>%
  # convert to factor to enforce ordering
  mutate(cutoff = factor(cutoff, levels = c(fs.cutoff, "full"))) %>%
  ggplot(aes(x = cutoff, y = AUC, col = type)) +
  geom_point() + geom_line(aes(group = type)) +
  facet_grid(~split) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0)) +
  xlab("Features selected") +
  ylab("AUROC") +
  theme_bw() + 
  scale_colour_manual(values = c("correct" = "blue", "incorrect" = "red"),
                      name = "Feature selection procedure") + 
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

### incorrect feature selection results in inflated AUROC values in CV
### incorrect feature selection results in lower generalization to truly external datasets
### both of these errors are larger when fewer features are selected


### NAIVE SPLITTING OF DEPENDENT DATA
# if samples are not independent (same individual, different time points) they will usually be more similar to each other than samples from other individuals
# if the samples from the same individual are split between training and testing datasets 
  # it will train the model to generalize across time points for the same individual
  # rather than training the model to distinguish the label across individuals

### dependent measurements need to be blocked during CV to ensure that the samples from the same individual reamin in the same fold

# load data (Crohn's disease data - filtered and cleaned)
data.loc <- "https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/"

# get feature data
feat.motus <- read.table(paste0(data.loc, "feat_rel_filt_cd.tsv"),
                         sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
str(feat.motus)

# get metadata
meta.all <- read_tsv(paste0(data.loc, "meta_all_cd.tsv"))
head(meta.all)

x <- meta.all %>% 
  group_by(Study, Group) %>% 
  summarise(n.all = n(), .groups = "drop")
y <- meta.all %>% 
  select(Study, Group, Individual_ID) %>% 
  distinct() %>% 
  group_by(Study, Group) %>% 
  summarize(n.indi = n(),  .groups = "drop")
full_join(x,y) # several individuals (especially in the HMP2 study) were sampled multiple times

# the number of samples per individual in the HMP2 study varies a lot
meta.all %>% 
  filter(Study == "HMP2") %>% 
  group_by(Individual_ID) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pull(n) %>% hist(20)

# randomly select a set of five samples per individual for training dataset
meta.train <- meta.all %>% 
  filter(Study == "HMP2") %>% 
  group_by(Individual_ID) %>%
  sample_n(5, replace = TRUE) %>%
  distinct() %>%
  as.data.frame()
rownames(meta.train) <- meta.train$Sample_ID

# for evaluation, only include one sample per individual
# create a new matrix from meta.all by removing repeated samples (keep only the earliest time point)
meta.ind <- meta.all %>% 
  group_by(Individual_ID) %>% 
  filter(Timepoint == min(Timepoint)) %>% 
  ungroup()

# create a tibble to hold the resulting AUROC values
auroc.all <- tibble(type=character(0), study.test = character(0), AUC = double(0))


### TRAIN WITH NAIVE CROSS-VALIDATION
# not taking into account the dependency between samples
sc.obj <- siamcat(feat = feat.motus, meta = meta.train, label = "Group", case = "CD")
sc.obj <- normalize.features(sc.obj, norm.method = "log.std",
                             norm.param = list(log.n0 = 1e-05, sd.min.q = 1), 
                             feature.type = "original")
sc.obj.naive <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)
sc.obj.naive <- train.model(sc.obj.naive, method = "lasso")
sc.obj.naive <- make.predictions(sc.obj.naive)
sc.obj.naive <- evaluate.predictions(sc.obj.naive)
auroc.all <- auroc.all %>% 
  add_row(type = "naive", study.test = "HMP2",
          AUC = as.numeric(eval_data(sc.obj.naive)$auroc))



### TRAIN WITH BLOCKED CROSS-VALIDATION
# training the model after blocking the CV procedure by individuals
# samples from the same individuals will always end up in the same fold
# specify the inseparable parameter in create.data.split()
sc.obj.block <- create.data.split(sc.obj, num.folds = 10, num.resample = 10,
                                  inseparable = "Individual_ID")
sc.obj.block <- train.model(sc.obj.block, method = "lasso")
sc.obj.block <- make.predictions(sc.obj.block)
sc.obj.block <- evaluate.predictions(sc.obj.block)
auroc.all <- auroc.all %>% 
  add_row(type = "blocked", study.test = "HMP2", 
          AUC = as.numeric(eval_data(sc.obj.block)$auroc))


### APPLY TO EXTERNAL DATASETS
# apply the naive (incorrect) and the blocked (correct) models to external datasets
for (i in setdiff(unique(meta.all$Study), "HMP2")){
  meta.test <- meta.ind %>% 
    filter(Study == i) %>% 
    as.data.frame()
  rownames(meta.test) <- meta.test$Sample_ID
  # apply naive model
  sc.obj.test <- siamcat(feat = feat.motus, meta = meta.test, 
                         label = "Group", case = "CD")
  sc.obj.test <- make.predictions(sc.obj.naive, sc.obj.test)
  sc.obj.test <- evaluate.predictions(sc.obj.test)
  auroc.all <- auroc.all %>% 
    add_row(type = "naive", study.test = i,
            AUC = as.numeric(eval_data(sc.obj.test)$auroc))
  # apply blocked model
  sc.obj.test <- siamcat(feat = feat.motus, meta = meta.test, 
                         label = "Group", case = "CD")
  sc.obj.test <- make.predictions(sc.obj.block, sc.obj.test)
  sc.obj.test <- evaluate.predictions(sc.obj.test)
  auroc.all <- auroc.all %>% 
    add_row(type = "blocked", study.test = i,
            AUC = as.numeric(eval_data(sc.obj.test)$auroc))
}


### PLOT THE RESULTS
# compare the AUROC values between the two different approaches
auroc.all %>%
  # convert to factor to enforce ordering
  mutate(type = factor(type, levels = c("naive", "blocked"))) %>%
  # facetting for plotting
  mutate(CV = case_when(study.test == "HMP2" ~ "CV", 
                      TRUE ~ "External validation")) %>%
  ggplot(aes(x = study.test, y = AUC, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(), col = "black") +
  theme_bw() +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_fill_manual(values = c("red", "blue"), name = "") +
  facet_grid(~CV, space = "free", scales = "free") +
  xlab("") + ylab("AUROC") +
  theme(legend.position = c(0.8, 0.8))

## naive CV resulted in inflated performance estimation compared to the blocked CV
## naive CV resulted in worse generalization to external datasets


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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ExperimentHub_2.16.0            AnnotationHub_3.16.0           
# [3] BiocFileCache_2.16.0            dbplyr_2.5.0                   
# [5] cowplot_1.1.3                   ggembl_0.1.2                   
# [7] ggpubr_0.6.0                    curatedMetagenomicData_3.16.1  
# [9] TreeSummarizedExperiment_2.16.1 Biostrings_2.76.0              
# [11] XVector_0.48.0                  SingleCellExperiment_1.30.1    
# [13] SummarizedExperiment_1.38.1     Biobase_2.68.0                 
# [15] GenomicRanges_1.60.0            GenomeInfoDb_1.44.0            
# [17] IRanges_2.42.0                  S4Vectors_0.46.0               
# [19] BiocGenerics_0.54.0             generics_0.1.4                 
# [21] MatrixGenerics_1.20.0           matrixStats_1.5.0              
# [23] SIAMCAT_2.12.0                  phyloseq_1.52.0                
# [25] mlr3_1.0.0                      lubridate_1.9.4                
# [27] forcats_1.0.0                   stringr_1.5.1                  
# [29] dplyr_1.1.4                     purrr_1.0.4                    
# [31] readr_2.1.5                     tidyr_1.3.1                    
# [33] tibble_3.3.0                    ggplot2_3.5.2                  
# [35] tidyverse_2.0.0                
# 
# loaded via a namespace (and not attached):
#   [1] ggtext_0.1.2                fs_1.6.6                    DirichletMultinomial_1.50.0
# [4] httr_1.4.7                  RColorBrewer_1.1-3          numDeriv_2016.8-1.1        
# [7] tools_4.5.0                 mlr3learners_0.12.0         backports_1.5.0            
# [10] utf8_1.2.6                  R6_2.6.1                    vegan_2.7-1                
# [13] lazyeval_0.2.2              mgcv_1.9-3                  rhdf5filters_1.20.0        
# [16] permute_0.9-7               withr_3.0.2                 prettyunits_1.2.0          
# [19] gridExtra_2.3               cli_3.6.5                   textshaping_1.0.1          
# [22] sandwich_3.1-1              labeling_0.4.3              slam_0.1-55                
# [25] mvtnorm_1.3-3               mlr3tuning_1.4.0            systemfonts_1.2.3          
# [28] yulab.utils_0.2.0           paradox_1.0.1               dichromat_2.0-0.1          
# [31] scater_1.35.0               decontam_1.28.0             parallelly_1.45.0          
# [34] readxl_1.4.5                fillpattern_1.0.2           rstudioapi_0.17.1          
# [37] RSQLite_2.4.1               shape_1.4.6.1               vroom_1.6.5                
# [40] car_3.1-3                   rbiom_2.2.0                 Matrix_1.7-3               
# [43] biomformat_1.36.0           ggbeeswarm_0.7.2            DECIPHER_3.4.0             
# [46] abind_1.4-8                 infotheo_1.2.0.1            lifecycle_1.0.4            
# [49] multcomp_1.4-28             yaml_2.3.10                 carData_3.0-5              
# [52] rhdf5_2.52.1                SparseArray_1.8.0           grid_4.5.0                 
# [55] blob_1.2.4                  LiblineaR_2.10-24           crayon_1.5.3               
# [58] lattice_0.22-7              beachmat_2.24.0             KEGGREST_1.48.0            
# [61] pillar_1.10.2               beanplot_1.3.1              boot_1.3-31                
# [64] estimability_1.5.1          codetools_0.2-20            glue_1.8.0                 
# [67] data.table_1.17.4           MultiAssayExperiment_1.34.0 vctrs_0.6.5                
# [70] png_0.1-8                   treeio_1.32.0               Rdpack_2.6.4               
# [73] cellranger_1.1.0            gtable_0.3.6                cachem_1.1.0               
# [76] rbibutils_2.3               S4Arrays_1.8.1              mime_0.13                  
# [79] coda_0.19-4.1               reformulas_0.4.1            survival_3.8-3             
# [82] iterators_1.0.14            bluster_1.18.0              TH.data_1.1-3              
# [85] nlme_3.1-168                bit64_4.6.0-1               bbotk_1.6.0                
# [88] progress_1.2.3              PRROC_1.4                   filelock_1.0.3             
# [91] irlba_2.3.5.1               vipor_0.4.7                 mlr3measures_1.0.0         
# [94] colorspace_2.1-1            DBI_1.2.3                   ade4_1.7-23                
# [97] mlr3misc_0.18.0             tidyselect_1.2.1            emmeans_1.11.1             
# [100] bit_4.6.0                   compiler_4.5.0              curl_6.3.0                 
# [103] glmnet_4.1-9                BiocNeighbors_2.2.0         lgr_0.4.4                  
# [106] xml2_1.3.8                  DelayedArray_0.34.1         checkmate_2.3.2            
# [109] scales_1.4.0                rappdirs_0.3.3              palmerpenguins_0.1.1       
# [112] digest_0.6.37               minqa_1.2.8                 pkgconfig_2.0.3            
# [115] lme4_1.1-37                 sparseMatrixStats_1.20.0    fastmap_1.2.0              
# [118] rlang_1.1.6                 UCSC.utils_1.4.0            DelayedMatrixStats_1.30.0  
# [121] farver_2.1.2                zoo_1.8-14                  jsonlite_2.0.0             
# [124] BiocParallel_1.42.1         BiocSingular_1.24.0         magrittr_2.0.3             
# [127] Formula_1.2-5               scuttle_1.18.0              GenomeInfoDbData_1.2.14    
# [130] patchwork_1.3.0             Rhdf5lib_1.30.0             Rcpp_1.0.14                
# [133] ape_5.8-1                   ggnewscale_0.5.1            viridis_0.6.5              
# [136] stringi_1.8.7               pROC_1.18.5                 MASS_7.3-65                
# [139] plyr_1.8.9                  parallel_4.5.0              listenv_0.9.1              
# [142] ggrepel_0.9.6               splines_4.5.0               gridtext_0.1.5             
# [145] multtest_2.64.0             hms_1.1.3                   igraph_2.1.4               
# [148] uuid_1.2-1                  ggsignif_0.6.4              pkgload_1.4.0              
# [151] reshape2_1.4.4              ScaledMatrix_1.16.0         BiocVersion_3.21.1         
# [154] BiocManager_1.30.26         nloptr_2.2.1                tzdb_0.5.0                 
# [157] foreach_1.5.2               future_1.58.0               gridBase_0.4-7             
# [160] BiocBaseUtils_1.10.0        rsvd_1.0.5                  broom_1.0.8                
# [163] xtable_1.8-4                tidytree_0.4.6              rstatix_0.7.2              
# [166] viridisLite_0.4.2           ragg_1.4.0                  lmerTest_3.1-3             
# [169] memoise_2.0.1               beeswarm_0.4.0              AnnotationDbi_1.70.0       
# [172] cluster_2.1.8.1             corrplot_0.95               timechange_0.3.0           
# [175] globals_0.18.0              mia_1.16.0  

