############################################
### SIAMCAT - INFLAMMATORY BOWEL DISEASE ###
############################################

# Re-creating the figures in Microbiome meta-analysis and cross-disease comparison enabled by the SIAMCAT machine learning toolbox
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02306-1
# microbiome data necessitates data normalization for ML algorithms to work well
   # microbiome data is: compositional, zero-inflated, and non-Gaussian distribution
# SIAMCAT input: feature matrix (abundance of microbial taxa), group label, and meta-variables (e.g., demographics, clinical measurements)


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(SIAMCAT)
library(tidyverse)
library(ggpubr)
library(ggembl)
library(cowplot)
library(curatedMetagenomicData)

# load data
nielsen_data <- curatedMetagenomicData("NielsenHB_2014.marker_abundance", dryrun = FALSE)
nielsen_se <- mergeData(nielsen_data) # merge into a SummarizedExperiment object
meta.nielsen.full <- as.data.frame(colData(nielsen_se)) # extract sample metadata

head(meta.nielsen.full) # view the metadata
table(meta.nielsen.full$country)
table(meta.nielsen.full$disease, meta.nielsen.full$disease_subtype)

# there are more samples than individual subjects
length(unique(meta.nielsen.full$subject_id)) # 318 | subjectID changed to subject_id
nrow(meta.nielsen.full) # 396


# remove repeated samples for the same subject
# for subjects that were sampled multiple times, only take the sample from the first visit
colnames(meta.nielsen.full)
head(meta.nielsen.full)
table(meta.nielsen.full$days_from_first_collection) # sampleID changed to days_from_first_collection

meta.nielsen <- meta.nielsen.full %>%
  select(days_from_first_collection, subject_id, study_condition, disease_subtype,
         disease, age, country, number_reads, median_read_length, BMI) %>%
  mutate(visit = case_when(is.na(days_from_first_collection) ~ 0, TRUE ~ days_from_first_collection)) %>%
  group_by(subject_id) %>%
  filter(visit == min(visit)) %>%
  ungroup() %>%
  mutate(Group = case_when(disease == "healthy" ~ "CTR", TRUE ~ disease_subtype))

length(unique(meta.nielsen$subject_id)) # 318
View(meta.nielsen)


### TAXONOMIC PROFILES
# load taxonomic profiles generated with MetaPhlAn2
# x <- "2021-03-31.NielsenHB_2014.relative_abundance"
# no longer available on curatedMetagenomicData
# on ExperimentHub, NielsenHB_2014.marker_abundance doens't have visit info, so need to use the version on curatedMetagenomicData

# load ExperimentHub
library(ExperimentHub)
eh <- ExperimentHub()

query(eh, "NielsenHB_2014.metaphlan_bugs_list.stool") # EH301
eh_bugs <- "EH301" # fetch data
feat <- eh[[eh_bugs]]@assayData$exprs

# only retain the species level profiles
feat <- feat[grep(x = rownames(feat), pattern = "s__"),]
feat <- feat[grep(x = rownames(feat),pattern="t__", invert = TRUE),]
feat <- t(t(feat)/100) # convert to relative percentages
feat <- as.data.frame(feat)

# shorten the species names
rownames(feat) <- str_extract(rownames(feat), "s__.*$")


###########################
### SPANISH COHORT ONLY ###
###########################

# restrict dataset to only the Spanish samples 
# restrict dataset to only control (CTR) and ulcerative colitis (UC)
table(meta.nielsen$country)
table(meta.nielsen$Group)

meta.nielsen.uc <- meta.nielsen %>%
  filter(Group %in% c("CTR", "UC")) %>%
  filter(country == "ESP") %>%
  as.data.frame()
rownames(meta.nielsen.uc) <- meta.nielsen.uc$subject_id
meta.nielsen.uc$days_from_first_collection <- NULL

View(meta.nielsen.uc)


### CREATE THE SIAMCAT OBJECT
# stores the feature matrix, the meta-variables and the label (Group from metadata)
all(order(rownames(meta.nielsen.uc)) %in% order(colnames(feat)))

# reformat data (need to get the data from two different sources)
feat <- feat[, order(names(feat))] # order dataframe by colnames
meta.nielsen.uc1 <- meta.nielsen.uc[order(rownames(meta.nielsen.uc)), ] # order dataframe by rownames
colnames(feat) <- sub("^([^_]+_[^_]+)_[^_]+$", "\\1", colnames(feat)) # match nomenclature
all(rownames(meta.nielsen.uc) %in% colnames(feat))

sc.obj <- siamcat(feat = feat, meta = meta.nielsen.uc, label = "Group", case = "UC")


### FEATURE FILTERING
# filter features with low overall abundance and prevalence 
sc.obj <- filter.features(sc.obj, cutoff = 1e-04, filter.method = "abundance")
sc.obj <- filter.features(sc.obj, cutoff = 0.05, filter.method = "prevalence", 
                          feature.type = "filtered")


### ASSOCIATION PLOT (FIGURE 1B)
# univariate association testing results 
# distribution of microbial abundance data differing significantly between the groups 
# calculates the significance of enrichment after multiple testing correction
# calculates metrics of association: generalized fold change (non-parametric measure of effect size) and single-feature AUROC

sc.obj <- check.associations(sc.obj, alpha = 0.1)
association.plot(sc.obj, plot.type = "quantile.rect", fn.plot = "./association_plot.pdf")


### CONFOUNDER ANALYSIS (FIGURE 1C)
# check the metadata for potential confounding
# statistical testing and diagnostic visualizations to identify potential confounders 
# tests for associations between meta-variables and disease labels
# produces one plot for each meta-variable
check.confounders(sc.obj, fn.plot = "./confounders.pdf")


### MACHINE LEARNING WORKFLOW
# 1. feature normalization
sc.obj <- normalize.features(sc.obj, norm.method = "log.std", 
                             norm.param = list(log.n0 = 1e-06, sd.min.q = 0))

# 2. data splitting for cross-validation
sc.obj <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)

# 3. model training
sc.obj <- train.model(sc.obj, method = "lasso")

# 4.model predictions (from left-out data)
sc.obj <- make.predictions(sc.obj)

# 5. evaluate model predictions (using AUROC and AUPRC)
sc.obj <- evaluate.predictions(sc.obj)


# MODEL EVALUATION (FIGURE 1D)
# ROC curve + precision-recall curve
# cross-validation error as a receiver operating characteristic curve (ROC) | 95% CI shaded in gray
model.evaluation.plot(sc.obj, fn.plot = "./eval_plot.pdf", show.all = TRUE)
# show.all = FALSE to re-create figure in paper


# MODEL INTERPRETATION PLOT (FIGURE 1E)
# feature importance barplot
# feature robustness (in the repeated CV, percent indicates the fraction of models containing that feature)
# normalized feature abundances across all samples (heatmap)
   # classification result (test predictions) and user defined meta-variables
# proportion of the model weight that is explained by the selected features (boxplot)

model.interpretation.plot(sc.obj, consens.thres = 0.8, fn.plot = "./interpretation.pdf")


##################################################
### SPANISH COHORT ONLY - GGPLOT2 FOR FIGURE 1 ###
##################################################

# make data.frame with species names and results from association testing
temp <- associations(sc.obj)
temp$species <- rownames(temp)
temp <- temp %>%
  mutate(species_short = str_remove(species, ".*s__")) %>%
  filter(p.adj < 0.1) %>%
  arrange(fc) %>%
  mutate(species = factor(species, levels = species))

# plot pf generalized fold change as a non-parametric measure of effect size
g.fc <- temp %>%
  ggplot(aes(x = species, y = fc, fill = fc > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_publication() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("gFC") + xlab("") +
  scale_fill_manual(values = c("#C8A2C885", "#FFA07A85"), guide = FALSE)

# plot significance (after multiple testing correction) as horizontal bars
g.pval <- temp %>%
  ggplot(aes(x = species, y = -log10(p.adj))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_publication() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("-") + xlab("")  +
  geom_hline(yintercept = -log10(0.1), colour = "red")

# create feature matrix with subject ids
df.label <- enframe(sc.obj@label$label, name = "SampleID", value = "Label")
mat.ab <- get.filt_feat.matrix(sc.obj)[as.character(temp$species),] %>%
  as.data.frame()
mat.ab$species <- rownames(mat.ab)

# create data.frame with species, sampleid and abundance
df.plot <- mat.ab %>%
  mutate(species = str_remove(species, ".*s__")) %>%
  gather(key = "SampleID", value = "abundance", -species) %>%
  mutate(abundance = log10(abundance+1e-06)) %>%
  mutate(species = factor(species, levels = temp$species_short))

# plot of the distributions of microbial abundance data differing significantly between the groups
g.ab <- df.plot %>%
  full_join(df.label) %>%
  mutate(Label = as.factor(Label)) %>%
  mutate(species = factor(species, levels = temp$species_short)) %>%
  ggplot(aes(x = species, y = abundance, fill = Label)) +
  geom_boxplot(outlier.shape = NA, size = 0.1) +
  geom_jitter(aes(col = Label), position = position_jitterdodge(jitter.width = 0.08), size = 0.6, stroke = 0) +
  coord_flip() + theme_publication() + 
  scale_fill_manual(values = c("#C8A2C885", "#FFA07A85"), guide = FALSE) +
  scale_colour_manual(values = c("#C8A2C8", "#FFA07A"), guide = FALSE) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  ylab("log10(rel. ab.)") + xlab("")

### FIGURE 1B: univariate association testing results
g.all <- plot_grid(g.ab, g.pval, g.fc, nrow = 1, rel_widths = c(0.5, 0.25, 0.25))


### FIGURE 1C: comparison of BMI between study groups
g.bmi <- meta(sc.obj) %>%
  ggplot(aes(x = study_condition, y = BMI)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, stroke = 0, size = 1.3) +
  theme_publication() +
  xlab("") + ylab("BMI")

# create data.frame of feature weights
feat.weights <- feature_weights(sc.obj)
feat.weights$species <- rownames(feat.weights)
rownames(feat.weights) <- NULL
feat.weight.selected <- feat.weights %>%
  mutate(species = str_remove(species, ".*s__")) %>%
  filter(percentage > 0.8) %>%
  arrange(median.rel.weight)

# FIGURE 1D: CV error as a ROC curve (95% confidence interval in gray and AUROC below the curve
model.evaluation.plot(sc.obj, fn.plot = "./roc.pdf", show.all = TRUE)

### FIGURE 1E: barplot of feature importance for the features that are included in the majority of models fitted during CV 
# percentages indicate the respective fraction of models containing the feature
g.weights <- feat.weight.selected %>%
  mutate(species = factor(species, levels = rev(species))) %>%
  ggplot(aes(x = species, y = -median.rel.weight)) +
  geom_bar(stat = "identity") +
  theme_publication() + coord_flip() +
  xlab("") + ylab("median rel. feat weight") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_text(aes(label = paste0(percentage*100, "%"), y = 0.08))


##################################
### SPANISH AND DANISH COHORTS ###
##################################

# include both the Spanish and Danish samples 
# restrict dataset to only control (CTR) and ulcerative colitis (UC)
table(meta.nielsen$country)
table(meta.nielsen$Group)

meta.nielsen.uc.all <- meta.nielsen %>%
  filter(Group %in% c("CTR", "UC")) %>%
  as.data.frame()
rownames(meta.nielsen.uc.all) <- meta.nielsen.uc.all$subject_id
meta.nielsen.uc.all$days_from_first_collection <- NULL

head(meta.nielsen.uc.all)


### CREATE THE SIAMCAT OBJECT
sc.obj.all <- siamcat(feat = feat, meta = meta.nielsen.uc.all, label = "Group", case = "UC")


### FEATURE FILTERING
# filter features with low overall abundance and prevalence 
sc.obj.all <- filter.features(sc.obj.all, cutoff = 1e-04, filter.method = "abundance")
sc.obj.all <- filter.features(sc.obj.all, cutoff = 0.05, filter.method = "prevalence",
                              feature.type = "filtered")


### CONFOUNDER ANALYSIS
check.confounders(sc.obj.all, fn.plot = "./confounders_all.pdf")
# meta-variable country might cause problems


### ASSOCIATION TESTING
sc.obj.all <- check.associations(sc.obj.all, alpha = 0.1)
association.plot(sc.obj.all, plot.type = "quantile.rect", fn.plot = "./association_plot_all.pdf")


# confounders can lead to biases in association testing
# extract association metrics from the SIAMCAT object after checking for associations in only ESP and in ESP + DNK
assoc.esp <- associations(sc.obj)
assoc.esp$species <- rownames(assoc.esp)

assoc.all <- associations(sc.obj.all)
assoc.all$species <- rownames(assoc.all)

df.plot <- full_join(assoc.esp, assoc.all, by = "species")
df.plot %>%
  mutate(highlight = str_detect(species, "Dorea_formicigenerans")) %>%
  ggplot(aes(x = -log10(p.adj.x), y = -log10(p.adj.y), col = highlight)) +
  geom_point(alpha = 0.3) +
  xlab("Spanish samples only\n-log10(q)") +
  ylab("Both Spanish and Danish samples\n-log10(q)") +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "lightgrey"),
        aspect.ratio = 1.3) +
  scale_colour_manual(values = c("darkgrey", "#D41645"), guide = "none") +
  annotate("text", x = 0.7, y = 8, label = "Dorea formicigenerans")

# the plot demonstrates that several species are only significant if the Danish control samples are included
   # e.g., D. formicigenerans is not significant in the Spanish cohort, but is highly significant when the Danish samples are included


### EXTRACT INFORMATION FROM THE SIAMCAT OBJECT
# plot log10 of D. formicigenerans in DNK_CTR, ESP_CTR and ESP_UC
feat.all <- get.filt_feat.matrix(sc.obj.all)
label.all <- label(sc.obj.all)$label
country <- meta(sc.obj.all)$country
names(country) <- rownames(meta(sc.obj.all))

df.plot <- tibble(dorea = log10(feat.all["s__Dorea_formicigenerans", names(label.all)] + 1e-05),
                  label = label.all, country = country) %>%
  mutate(label = case_when(label == "-1" ~ "CTR", TRUE ~ "UC")) %>%
  mutate(x_value = paste0(country, "_", label))

df.plot %>%
  ggplot(aes(x = x_value, y = dorea)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, stroke = 0, alpha = 0.2) + theme_classic() + xlab("") +
  ylab("log10(Dorea_formicigenerans)") + theme(aspect.ratio = 1.3) +
  stat_compare_means(comparisons = list(c("DNK_CTR", "ESP_CTR"),
                                        c("DNK_CTR", "ESP_UC"),
                                        c("ESP_CTR", "ESP_UC")))


### MACHINE LEARNING
# the results from machine learning workflows can be biased by confounders (e.g., country)
   # leads to exaggerated performance estimates

# 1. feature normalization
sc.obj.all <- normalize.features(sc.obj.all, norm.method = "log.std",
                                 norm.param = list(log.n0 = 1e-06, sd.min.q = 0))

# 2. data splitting for cross-validation
sc.obj.all <- create.data.split(sc.obj.all, num.folds = 10, num.resample = 10)

# 3. model training
sc.obj.all <- train.model(sc.obj.all, method = "lasso")

# 4. model predictions (from left-out data)
sc.obj.all <- make.predictions(sc.obj.all)

# 5. evaluate model predictions (using AUROC and AUPRC)
sc.obj.all <- evaluate.predictions(sc.obj.all)


# MODEL EVALUATION 
# ROC curve + precision-recall curve
model.evaluation.plot("Spanish samples only" = sc.obj, 
                      "Danish and Spanish samples" = sc.obj.all,
                      fn.plot = "./eval_plot_all.pdf", show.all = TRUE)

# the model with the samples from both ESP and DNK seems to perform better
   # higher AUROC value
# previous analysis tells us that this performance estimate is biased and exaggerated because of differences between ESP and DNK samples


# TRAIN THE MODEL TO DISTINGUISH BETWEEN THE CONTROL SAMPLES BETWEEN ESP AND DNK

# include data from both countries, but only from the control group
meta.nielsen.country <- meta.nielsen %>%
  filter(Group %in% c("CTR")) %>% # only control samples
  as.data.frame()
rownames(meta.nielsen.country) <- meta.nielsen.country$subject_id
meta.nielsen.country$days_from_first_collection <- NULL

# create the SIAMCAT object
sc.obj.country <- siamcat(feat = feat, meta = meta.nielsen.country,
                          label = "country", case = "ESP")


# filter features 
sc.obj.country <- filter.features(sc.obj.country, cutoff = 1e-04, filter.method = "abundance")
sc.obj.country <- filter.features(sc.obj.country, cutoff = 0.05,
                                  filter.method = "prevalence", feature.type = "filtered")

# normalize features
sc.obj.country <- normalize.features(sc.obj.country, norm.method = "log.std",
                                     norm.param = list(log.n0 = 1e-06, sd.min.q = 0))

# split object for cross-validation
sc.obj.country <- create.data.split(sc.obj.country, num.folds = 10, num.resample = 10)

# train model
sc.obj.country <- train.model(sc.obj.country, method = "lasso")

# make and evaluate predictions
sc.obj.country <- make.predictions(sc.obj.country)
sc.obj.country <- evaluate.predictions(sc.obj.country)
model.evaluation.plot(sc.obj.country, fn.plot = "./eval_plot_country.pdf", show.all = TRUE)

# the model can distinguish between the two countries with almost perfect accuracy


############################################################
### SPANISH AND DANISH COMPARISON - GGPLOT2 FOR FIGURE 2 ###
############################################################

feat <- get.filt_feat.matrix(sc.obj.all) # feature table for ESP and DNK
label <- label(sc.obj.all)$label # subjectid vector for ESP and DNK
country <- as.factor(meta(sc.obj.all)$country) # country vector
names(country) <- rownames(meta(sc.obj.all)) # country and subjectid vector

# species abundance for each subject_id
var.label <- vapply(rownames(feat), FUN = function(x){
  x <- feat[x,]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(unique(label), function(s){
    sum((x[label == s] - mean(x[label == s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1-ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.label))){
  var.label[is.infinite(var.label)] <- NA
}

# batch effect measure for each bacteria (variance explained by country grouping)
var.batch <- vapply(rownames(feat), FUN = function(x){
  x <- feat[x, names(country)]
  x <- rank(x)/length(x)
  ss.tot <- sum((x - mean(x))^2)/length(x)
  ss.o.i <- sum(vapply(levels(country), function(s){
    sum((x[country == s] - mean(x[country == s]))^2)
  }, FUN.VALUE = double(1)))/length(x)
  return(1 - ss.o.i/ss.tot)
}, FUN.VALUE = double(1))
if (any(is.infinite(var.batch))){
  var.batch[is.infinite(var.batch)] <- NA
}

# FIGURE 2B: variance explained by country and by disease
df.plot <- tibble(label = var.label, batch = var.batch, species = names(var.label))
temp <- get.filt_feat.matrix(sc.obj.all)
df.plot$mean <- rowMeans(log10(temp + 1e-05))
df.plot <- df.plot %>%
  mutate(species = str_remove(species, "^.*s__"))
g.conf <- df.plot %>%
  ggplot(aes(x = label, y = batch, size = mean)) +
  geom_point(stroke = 0, alpha = 0.2) +
  theme_publication() + theme(aspect.ratio = 1) +
  xlim(0, 0.3) + ylim(0, 0.3) +
  scale_size_continuous(range = c(0.2, 4))

### FIGURE 2E and 2F: machine learning models is confounded by country
model.evaluation.plot(sc.obj.all, fn.plot = "./roc_all.pdf", show.all = TRUE)
model.evaluation.plot(sc.obj.country, fn.plot = "./roc_country.pdf", show.all = TRUE)


### FIGURE 2C: comparison of padj of SIAMCAT associations between ESP and DNK dataset and just the ESP dataset 
temp <- associations(sc.obj)
temp$species <- rownames(temp)
temp2 <- associations(sc.obj.all)
temp2$species <- rownames(temp2)
df.plot <- full_join(temp, temp2, by = "species")

g.assoc <- df.plot %>%
  ggplot(aes(x = -log10(p.adj.x), y = -log10(p.adj.y))) +
  geom_abline(intercept = 0, slope = 1, col = "darkgrey", linetype = 3) +
  geom_point(alpha = 0.2, stroke = 0) + theme_publication() +
  xlab("Spanish samples only\n-log10(q val)") +
  ylab("Spanish and Danish samples\n-log10(q val)")

# FIGURE 2D: relative abundance differences for D. formicigenerans between DNK_CTR, ESP_CTR and ESP_UC
temp <- get.filt_feat.matrix(sc.obj.all)
label <- label(sc.obj.all)$label
country <- meta(sc.obj.all)[["country"]]
x <- which(str_detect(rownames(temp), "Dorea_formicigenerans"))

df.plot <- tibble(bug = log10(temp[x,names(label)] + 1e-05), label = label, country = country) %>%
  mutate(x_value = paste0(country, "_", label))

g <- df.plot %>%
  ggplot(aes(x = x_value, y = bug)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.08, stroke = 0, alpha = 0.2) +
  theme_publication() + xlab("") + ylab("log10(Dorea_formicigenerans)") +
  stat_compare_means(comparisons = list(c("DNK_-1", "ESP_-1"),
                                        c("DNK_-1", "ESP_1"),
                                        c("ESP_-1", "ESP_1")))
wilcox.test(df.plot$bug ~ df.plot %>% 
              mutate(x_value = case_when(x_value == "ESP_-1" ~ "DNK_-1", TRUE~x_value)) %>% 
              pull(x_value)) # right Wilcoxon test


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
# [5] curatedMetagenomicData_3.16.1   TreeSummarizedExperiment_2.16.1
# [7] Biostrings_2.76.0               XVector_0.48.0                 
# [9] SingleCellExperiment_1.30.1     SummarizedExperiment_1.38.1    
# [11] Biobase_2.68.0                  GenomicRanges_1.60.0           
# [13] GenomeInfoDb_1.44.0             IRanges_2.42.0                 
# [15] S4Vectors_0.46.0                BiocGenerics_0.54.0            
# [17] generics_0.1.4                  MatrixGenerics_1.20.0          
# [19] matrixStats_1.5.0               cowplot_1.1.3                  
# [21] ggembl_0.1.2                    ggpubr_0.6.0                   
# [23] lubridate_1.9.4                 forcats_1.0.0                  
# [25] stringr_1.5.1                   dplyr_1.1.4                    
# [27] purrr_1.0.4                     readr_2.1.5                    
# [29] tidyr_1.3.1                     tibble_3.3.0                   
# [31] ggplot2_3.5.2                   tidyverse_2.0.0                
# [33] SIAMCAT_2.12.0                  phyloseq_1.52.0                
# [35] mlr3_1.0.0                     
# 
# loaded via a namespace (and not attached):
#   [1] ggtext_0.1.2                fs_1.6.6                    DirichletMultinomial_1.50.0
# [4] httr_1.4.7                  RColorBrewer_1.1-3          numDeriv_2016.8-1.1        
# [7] tools_4.5.0                 mlr3learners_0.12.0         backports_1.5.0            
# [10] utf8_1.2.6                  R6_2.6.1                    vegan_2.7-1                
# [13] lazyeval_0.2.2              mgcv_1.9-3                  rhdf5filters_1.20.0        
# [16] permute_0.9-7               withr_3.0.2                 prettyunits_1.2.0          
# [19] gridExtra_2.3               textshaping_1.0.1           cli_3.6.5                  
# [22] sandwich_3.1-1              labeling_0.4.3              slam_0.1-55                
# [25] mvtnorm_1.3-3               systemfonts_1.2.3           mlr3tuning_1.4.0           
# [28] yulab.utils_0.2.0           paradox_1.0.1               dichromat_2.0-0.1          
# [31] scater_1.35.0               decontam_1.28.0             parallelly_1.45.0          
# [34] readxl_1.4.5                fillpattern_1.0.2           rstudioapi_0.17.1          
# [37] RSQLite_2.4.1               shape_1.4.6.1               car_3.1-3                  
# [40] rbiom_2.2.0                 Matrix_1.7-3                biomformat_1.36.0          
# [43] ggbeeswarm_0.7.2            DECIPHER_3.4.0              abind_1.4-8                
# [46] infotheo_1.2.0.1            lifecycle_1.0.4             multcomp_1.4-28            
# [49] yaml_2.3.10                 carData_3.0-5               rhdf5_2.52.1               
# [52] SparseArray_1.8.0           grid_4.5.0                  blob_1.2.4                 
# [55] LiblineaR_2.10-24           crayon_1.5.3                lattice_0.22-7             
# [58] beachmat_2.24.0             KEGGREST_1.48.0             pillar_1.10.2              
# [61] beanplot_1.3.1              boot_1.3-31                 estimability_1.5.1         
# [64] codetools_0.2-20            glue_1.8.0                  data.table_1.17.4          
# [67] MultiAssayExperiment_1.34.0 vctrs_0.6.5                 png_0.1-8                  
# [70] treeio_1.32.0               Rdpack_2.6.4                cellranger_1.1.0           
# [73] gtable_0.3.6                cachem_1.1.0                mime_0.13                  
# [76] rbibutils_2.3               S4Arrays_1.8.1              coda_0.19-4.1              
# [79] reformulas_0.4.1            survival_3.8-3              iterators_1.0.14           
# [82] bluster_1.18.0              TH.data_1.1-3               nlme_3.1-168               
# [85] bit64_4.6.0-1               bbotk_1.6.0                 progress_1.2.3             
# [88] PRROC_1.4                   filelock_1.0.3              irlba_2.3.5.1              
# [91] vipor_0.4.7                 mlr3measures_1.0.0          colorspace_2.1-1           
# [94] DBI_1.2.3                   ade4_1.7-23                 mlr3misc_0.18.0            
# [97] tidyselect_1.2.1            emmeans_1.11.1              bit_4.6.0                  
# [100] compiler_4.5.0              curl_6.3.0                  glmnet_4.1-9               
# [103] BiocNeighbors_2.2.0         lgr_0.4.4                   xml2_1.3.8                 
# [106] DelayedArray_0.34.1         checkmate_2.3.2             scales_1.4.0               
# [109] rappdirs_0.3.3              palmerpenguins_0.1.1        digest_0.6.37              
# [112] minqa_1.2.8                 pkgconfig_2.0.3             lme4_1.1-37                
# [115] sparseMatrixStats_1.20.0    fastmap_1.2.0               rlang_1.1.6                
# [118] UCSC.utils_1.4.0            DelayedMatrixStats_1.30.0   farver_2.1.2               
# [121] zoo_1.8-14                  jsonlite_2.0.0              BiocParallel_1.42.1        
# [124] BiocSingular_1.24.0         magrittr_2.0.3              Formula_1.2-5              
# [127] scuttle_1.18.0              GenomeInfoDbData_1.2.14     patchwork_1.3.0            
# [130] Rhdf5lib_1.30.0             Rcpp_1.0.14                 ape_5.8-1                  
# [133] ggnewscale_0.5.1            viridis_0.6.5               stringi_1.8.7              
# [136] pROC_1.18.5                 MASS_7.3-65                 plyr_1.8.9                 
# [139] parallel_4.5.0              listenv_0.9.1               ggrepel_0.9.6              
# [142] splines_4.5.0               gridtext_0.1.5              multtest_2.64.0            
# [145] hms_1.1.3                   igraph_2.1.4                uuid_1.2-1                 
# [148] ggsignif_0.6.4              reshape2_1.4.4              ScaledMatrix_1.16.0        
# [151] BiocVersion_3.21.1          BiocManager_1.30.26         nloptr_2.2.1               
# [154] tzdb_0.5.0                  foreach_1.5.2               future_1.58.0              
# [157] gridBase_0.4-7              BiocBaseUtils_1.10.0        rsvd_1.0.5                 
# [160] broom_1.0.8                 xtable_1.8-4                tidytree_0.4.6             
# [163] rstatix_0.7.2               ragg_1.4.0                  viridisLite_0.4.2          
# [166] lmerTest_3.1-3              memoise_2.0.1               beeswarm_0.4.0             
# [169] AnnotationDbi_1.70.0        cluster_2.1.8.1             corrplot_0.95              
# [172] timechange_0.3.0            globals_0.18.0              mia_1.16.0

