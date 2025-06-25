# IBD Meta-analysis using SIAMCAT
# metagenomics datasets from five studies of Crohn's disease (CD)
# https://github.com/zellerlab/siamcat_paper/blob/master/ibd_meta_analysis/ibd_meta_analysis.R
# https://siamcat.embl.de/articles/SIAMCAT_meta.html


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(SIAMCAT)
library(tidyverse)

# load data 
# raw data was taxonomically classified using mOTUs and aggregated at genus level
data.loc <- "https://zenodo.org/api/files/d81e429c-870f-44e0-a44a-2a4aa541b6c1/"
datasets <- c("metaHIT", "Lewis_2015", "He_2017", "Franzosa_2019", "HMP2")

meta.all <- read_tsv(paste0(data.loc, "meta_all_cd.tsv"))
head(meta.all)
length(meta.all$Individual_ID) # 1597

feat <- read.table(paste0(data.loc, "feat_genus_cd.tsv"), check.names = FALSE, 
                   stringsAsFactors = FALSE, quote = "", sep = "\t")
feat <- as.matrix(feat)
str(feat) 

# check that column names in feature table are the same as the SampleIDs in the metadata
all(colnames(feat) == meta.all$Sample_ID)

# distribution of control and CD groups across the studies
table(meta.all$Study, meta.all$Group)

# create a new metadata table retaining only one sample for each subject (first sample when there are multiple timepoints) 
meta.ind <- meta.all %>% 
  group_by(Individual_ID) %>% 
  filter(Timepoint == min(Timepoint)) %>% 
  ungroup()
length(meta.ind$Individual_ID) # 504
table(meta.ind$Study, meta.ind$Group)


### COMPARE ASSOCIATIONS
### create individual SIAMCAT objects for each dataset and extract associations 

assoc.list <- list()
for (d in datasets){
  # filter metadata and convert to dataframe
  meta.train <- meta.ind %>% 
    filter(Study == d) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  
  # create SIAMCAT object
  sc.obj <- siamcat(feat = feat, meta = meta.train, label = "Group", case = "CD")
  # test for associations
  sc.obj <- check.associations(sc.obj, log.n0 = 1e-05, 
                               feature.type = "original")
  # extract the associations and save them in the assoc.list
  temp <- associations(sc.obj)
  temp$genus <- rownames(temp)
  assoc.list[[d]] <- temp %>% 
    select(genus, fc, auc, p.adj) %>% 
    mutate(Study = d)
}

# combine all associations 
df.assoc <- bind_rows(assoc.list)
length(df.assoc$genus) # 10830
df.assoc <- df.assoc %>% filter(genus != "unclassified")
head(df.assoc)


### extract features that are strongly associated with the disease (control, CD) in at least one study
# plot generalized fold change in a heatmap
genera.of.interest <- df.assoc %>% 
  group_by(genus) %>% 
  summarise(m = mean(auc), n.filt = any(auc < 0.25 | auc > 0.75), .groups = "keep") %>% 
  filter(n.filt) %>% 
  arrange(m)

df.assoc %>% 
  # take only genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # convert to factor to enforce an ordering by mean AUC
  mutate(genus = factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  # convert to factor to enforce ordering again
  mutate(Study = factor(Study, levels = datasets)) %>% 
  # annotate the cells in the heatmap with stars
  mutate(l = case_when(p.adj < 0.01 ~ "*", TRUE ~ "")) %>%  
  ggplot(aes(y = genus, x = Study, fill = fc)) + geom_tile() + 
  scale_fill_gradient2(low = "#3B6FB6", high = "#D41645", mid = "white", 
                       limits=c(-2.7, 2.7), name='Generalized\nfold change') + 
  theme_minimal() + geom_text(aes(label = l)) + theme(panel.grid = element_blank()) + 
  xlab("") + ylab("") + theme(axis.text = element_text(size = 6))


### EXAMINE STUDY AS A POTENTIAL CONFOUNDER
# examine how differences between the studies might influence the variance of specific genera

# create a single SIAMCAT object for all of the datasets
df.meta <- meta.ind %>% 
  as.data.frame()
rownames(df.meta) <- df.meta$Sample_ID
sc.obj <- siamcat(feat = feat, meta = df.meta, label = "Group", case = "CD")

# check for confounders
check.confounders(sc.obj, fn.plot = "./confounder_plot_cd_meta.pdf",
                  feature.type = "original")
### some genera are strongly impacted by differences between the studies 
### genera that show the most impact from disease do not shoe a lot of variance across studies


#####################################################
### MACHINE LEARNING META-ANALYSIS - LASSO MODELS ###
#####################################################


# train one model for each of the datasets, then apply it to the other datasets (holdout testing functionality of SIAMCAT)
# block the cross-validation of subjects with repeated samples to not bias the results

# create tibble to store all of the predictions
auroc.all.lasso <- tibble(study.train = character(0), 
                          study.test = character(0),
                          AUC = double(0))

# create a list to save the trained SIAMCAT objects
sc.list.lasso <- list()

for (i in datasets){
  # restrict to a single study
  meta.train <- meta.all %>% 
    filter(Study == i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  
  ## take into account repeated sampling by blocking CV by Individual_ID (for studies with repeated samples)
  block <- NULL
  if (i %in% c("metaHIT", "Lewis_2015", "HMP2")){
    block <- "Individual_ID"
    if (i == "HMP2"){ 
      # for HMP2 dataset, reduce number of repeated samples per subject (some have up to 20 samples) 
      meta.train <- meta.all %>% 
        filter(Study == "HMP2") %>% 
        group_by(Individual_ID) %>% 
        sample_n(5, replace = TRUE) %>% 
        distinct() %>% 
        as.data.frame()
      rownames(meta.train) <- meta.train$Sample_ID
    }
  }
  # create SIAMCAT object
  sc.obj.train <- siamcat(feat = feat, meta = meta.train, 
                          label = "Group", case = "CD")
  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = "log.std",
                                     norm.param = list(log.n0 = 1e-05, sd.min.q = 0),
                                     feature.type = "original")
  # split data
  sc.obj.train <- create.data.split(sc.obj.train,
                                    num.folds = 10, num.resample = 10, inseparable = block)
  # train LASSO model
  sc.obj.train <- train.model(sc.obj.train, method = "lasso")
  
  ## apply trained models to other datasets
  # loop through datasets again
  for (i2 in datasets){
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all.lasso <- auroc.all.lasso %>% 
        add_row(study.train = i, study.test = i,
                AUC = eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate cross-validation predictions (external datasets)
      # use meta.ind here, since we want only one sample per subject
      meta.test <- meta.ind %>% 
        filter(Study == i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Sample_ID
      sc.obj.test <- siamcat(feat = feat, meta = meta.test,
                             label = "Group", case = "CD")
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all.lasso <- auroc.all.lasso %>% 
        add_row(study.train = i, study.test = i2,
                AUC = eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list.lasso[[i]] <- sc.obj.train
}

# calculate the test average for each dataset
test.average.lasso <- auroc.all.lasso %>% 
  filter(study.train != study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC = mean(AUC), .groups = "drop") %>% 
  mutate(study.train = "Average")

# combine AUROC values with test averages
bind_rows(auroc.all.lasso, test.average.lasso) %>% 
  # highlight cross validation versus transfer results
  mutate(CV = study.train == study.test) %>%
  # for faceting
  mutate(split = case_when(study.train == "Average" ~ "Average", TRUE ~ "none")) %>% 
  mutate(split = factor(split, levels = c("none", "Average"))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train = factor(study.train, levels = c(datasets, "Average"))) %>% 
  mutate(study.test = factor(study.test, levels = c(rev(datasets), "Average"))) %>% 
  ggplot(aes(y = study.test, x = study.train, fill = AUC, size = CV, color = CV)) +
  geom_tile() + theme_minimal() +
  geom_text(aes_string(label = "format(AUC, digits = 2)"), col = "white", size = 2) +
  scale_fill_gradientn(colours = rev(c("darkgreen","forestgreen", "chartreuse3", "lawngreen", "yellow"))) +
  scale_x_discrete(position = "top") + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle = 45, hjust = .1), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank()) + 
  xlab("Training Set") + ylab("Test Set") + 
  scale_color_manual(values = c("#FFFFFF00", "grey"), guide = "none") + 
  scale_size_manual(values = c(0, 1), guide = "none") + 
  facet_grid(~split, scales = "free", space = "free")


### INVESTIGATE FEATURE WEIGHTS
# trained models are saved in sc.list
# extract the models weights and compare them to the associations calculated above

weight.list.lasso <- list()

for (d in datasets){
  sc.obj.train <- sc.list[[d]]
  # extract the feature weights out of the SIAMCAT object
  temp <- feature_weights(sc.obj.train)
  temp$genus <- rownames(temp)
  # save selected info in the weight.list.lasso
  weight.list.lasso[[d]] <- temp %>% 
    select(genus, median.rel.weight, mean.rel.weight, percentage) %>% 
    mutate(Study = d) %>% 
    mutate(r.med = rank(-abs(median.rel.weight)), 
           r.mean = rank(-abs(mean.rel.weight)))
}

# combine all feature weights into a single tibble
df.weights.lasso <- bind_rows(weight.list.lasso)
sum(df.weights.lasso$genus == "unclassified")
df.weights.lasso <- df.weights.lasso %>% filter(genus != "unclassified")


### plot a heatmap with the weights for the genera of interest (those that had their associations plotted above)
# compute absolute feature weights
abs.weights.lasso <- df.weights.lasso %>% 
  group_by(Study) %>% 
  summarise(sum.median = sum(abs(median.rel.weight)),
            sum.mean = sum(abs(mean.rel.weight)),
            .groups = "drop")

df.weights.lasso %>% 
  full_join(abs.weights.lasso) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight = median.rel.weight/sum.median) %>% 
  # only include genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med = case_when(r.med > 20~NA_real_, TRUE ~ r.med)) %>%
  # enforce the correct ordering by converting to factors 
  mutate(genus = factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  mutate(Study = factor(Study, levels = datasets)) %>% 
  ggplot(aes(y = genus, x = Study, fill = median.rel.weight)) + geom_tile() + 
  scale_fill_gradientn(colours=rev(c("#007A53", "#009F4D", "#6CC24A", "white",
                                     "#EFC06E", "#FFA300", "#BE5400")), 
    limits = c(-0.15, 0.15)) + theme_minimal() + 
  geom_text(aes(label = r.med), col = "black", size = 2) +
  theme(panel.grid = element_blank()) + xlab("") + ylab("") +
  theme(axis.text = element_text(size = 6))


###########################################################
### MACHINE LEARNING META-ANALYSIS - ELASTIC NET MODELS ###
###########################################################


# train one model for each of the datasets, then apply it to the other datasets (holdout testing functionality of SIAMCAT)
# block the cross-validation of subjects with repeated samples to not bias the results

# create tibble to store all of the predictions
auroc.all.enet<- tibble(study.train = character(0), 
                          study.test = character(0),
                          AUC = double(0))

# create a list to save the trained SIAMCAT objects
sc.list.enet<- list()

for (i in datasets){
  # restrict to a single study
  meta.train <- meta.all %>% 
    filter(Study == i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  
  ## take into account repeated sampling by blocking CV by Individual_ID (for studies with repeated samples)
  block <- NULL
  if (i %in% c("metaHIT", "Lewis_2015", "HMP2")){
    block <- "Individual_ID"
    if (i == "HMP2"){ 
      # for HMP2 dataset, reduce number of repeated samples per subject (some have up to 20 samples) 
      meta.train <- meta.all %>% 
        filter(Study == "HMP2") %>% 
        group_by(Individual_ID) %>% 
        sample_n(5, replace = TRUE) %>% 
        distinct() %>% 
        as.data.frame()
      rownames(meta.train) <- meta.train$Sample_ID
    }
  }
  # create SIAMCAT object
  sc.obj.train <- siamcat(feat = feat, meta = meta.train, 
                          label = "Group", case = "CD")
  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = "log.std",
                                     norm.param = list(log.n0 = 1e-05, sd.min.q = 0),
                                     feature.type = "original")
  # split data
  sc.obj.train <- create.data.split(sc.obj.train,
                                    num.folds = 10, num.resample = 10, inseparable = block)
  # train ELASTIC NET model
  sc.obj.train <- train.model(sc.obj.train, method = "enet")
  
  ## apply trained models to other datasets
  # loop through datasets again
  for (i2 in datasets){
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all.enet<- auroc.all.enet%>% 
        add_row(study.train = i, study.test = i,
                AUC = eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate cross-validation predictions (external datasets)
      # use meta.ind here, since we want only one sample per subject
      meta.test <- meta.ind %>% 
        filter(Study == i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Sample_ID
      sc.obj.test <- siamcat(feat = feat, meta = meta.test,
                             label = "Group", case = "CD")
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all.enet<- auroc.all.enet%>% 
        add_row(study.train = i, study.test = i2,
                AUC = eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list.enet[[i]] <- sc.obj.train
}

# calculate the test average for each dataset
test.average.enet<- auroc.all.enet%>% 
  filter(study.train != study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC = mean(AUC), .groups = "drop") %>% 
  mutate(study.train = "Average")

# combine AUROC values with test averages
bind_rows(auroc.all.enet, test.average.enet) %>% 
  # highlight cross validation versus transfer results
  mutate(CV = study.train == study.test) %>%
  # for faceting
  mutate(split = case_when(study.train == "Average" ~ "Average", TRUE ~ "none")) %>% 
  mutate(split = factor(split, levels = c("none", "Average"))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train = factor(study.train, levels = c(datasets, "Average"))) %>% 
  mutate(study.test = factor(study.test, levels = c(rev(datasets), "Average"))) %>% 
  ggplot(aes(y = study.test, x = study.train, fill = AUC, size = CV, color = CV)) +
  geom_tile() + theme_minimal() +
  geom_text(aes_string(label = "format(AUC, digits = 2)"), col = "white", size = 2) +
  scale_fill_gradientn(colours = rev(c("darkgreen","forestgreen", "chartreuse3", "lawngreen", "yellow"))) +
  scale_x_discrete(position = "top") + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle = 45, hjust = .1), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank()) + 
  xlab("Training Set") + ylab("Test Set") + 
  scale_color_manual(values = c("#FFFFFF00", "grey"), guide = "none") + 
  scale_size_manual(values = c(0, 1), guide = "none") + 
  facet_grid(~split, scales = "free", space = "free")


### INVESTIGATE FEATURE WEIGHTS
# trained models are saved in sc.list
# extract the models weights and compare them to the associations calculated above

weight.list.enet<- list()

for (d in datasets){
  sc.obj.train <- sc.list[[d]]
  # extract the feature weights out of the SIAMCAT object
  temp <- feature_weights(sc.obj.train)
  temp$genus <- rownames(temp)
  # save selected info in the weight.list.enet
  weight.list.enet[[d]] <- temp %>% 
    select(genus, median.rel.weight, mean.rel.weight, percentage) %>% 
    mutate(Study = d) %>% 
    mutate(r.med = rank(-abs(median.rel.weight)), 
           r.mean = rank(-abs(mean.rel.weight)))
}

# combine all feature weights into a single tibble
df.weights.enet<- bind_rows(weight.list.enet)
sum(df.weights.enet$genus == "unclassified")
df.weights.enet<- df.weights.enet%>% filter(genus != "unclassified")


### plot a heatmap with the weights for the genera of interest (those that had their associations plotted above)
# compute absolute feature weights
abs.weights.enet<- df.weights.enet%>% 
  group_by(Study) %>% 
  summarise(sum.median = sum(abs(median.rel.weight)),
            sum.mean = sum(abs(mean.rel.weight)),
            .groups = "drop")

df.weights.enet%>% 
  full_join(abs.weights.enet) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight = median.rel.weight/sum.median) %>% 
  # only include genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med = case_when(r.med > 20~NA_real_, TRUE ~ r.med)) %>%
  # enforce the correct ordering by converting to factors 
  mutate(genus = factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  mutate(Study = factor(Study, levels = datasets)) %>% 
  ggplot(aes(y = genus, x = Study, fill = median.rel.weight)) + geom_tile() + 
  scale_fill_gradientn(colours=rev(c("#007A53", "#009F4D", "#6CC24A", "white",
                                     "#EFC06E", "#FFA300", "#BE5400")), 
                       limits = c(-0.15, 0.15)) + theme_minimal() + 
  geom_text(aes(label = r.med), col = "black", size = 2) +
  theme(panel.grid = element_blank()) + xlab("") + ylab("") +
  theme(axis.text = element_text(size = 6))


#####################################################
### MACHINE LEARNING META-ANALYSIS - RIDGE MODELS ###
#####################################################


# train one model for each of the datasets, then apply it to the other datasets (holdout testing functionality of SIAMCAT)
# block the cross-validation of subjects with repeated samples to not bias the results

# create tibble to store all of the predictions
auroc.all.ridge<- tibble(study.train = character(0), 
                         study.test = character(0),
                         AUC = double(0))

# create a list to save the trained SIAMCAT objects
sc.list.ridge<- list()

for (i in datasets){
  # restrict to a single study
  meta.train <- meta.all %>% 
    filter(Study == i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  
  ## take into account repeated sampling by blocking CV by Individual_ID (for studies with repeated samples)
  block <- NULL
  if (i %in% c("metaHIT", "Lewis_2015", "HMP2")){
    block <- "Individual_ID"
    if (i == "HMP2"){ 
      # for HMP2 dataset, reduce number of repeated samples per subject (some have up to 20 samples) 
      meta.train <- meta.all %>% 
        filter(Study == "HMP2") %>% 
        group_by(Individual_ID) %>% 
        sample_n(5, replace = TRUE) %>% 
        distinct() %>% 
        as.data.frame()
      rownames(meta.train) <- meta.train$Sample_ID
    }
  }
  # create SIAMCAT object
  sc.obj.train <- siamcat(feat = feat, meta = meta.train, 
                          label = "Group", case = "CD")
  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = "log.std",
                                     norm.param = list(log.n0 = 1e-05, sd.min.q = 0),
                                     feature.type = "original")
  # split data
  sc.obj.train <- create.data.split(sc.obj.train,
                                    num.folds = 10, num.resample = 10, inseparable = block)
  # train RIDGE model
  sc.obj.train <- train.model(sc.obj.train, method = "enet")
  
  ## apply trained models to other datasets
  # loop through datasets again
  for (i2 in datasets){
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all.ridge<- auroc.all.ridge%>% 
        add_row(study.train = i, study.test = i,
                AUC = eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate cross-validation predictions (external datasets)
      # use meta.ind here, since we want only one sample per subject
      meta.test <- meta.ind %>% 
        filter(Study == i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Sample_ID
      sc.obj.test <- siamcat(feat = feat, meta = meta.test,
                             label = "Group", case = "CD")
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all.ridge<- auroc.all.ridge%>% 
        add_row(study.train = i, study.test = i2,
                AUC = eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list.ridge[[i]] <- sc.obj.train
}

# calculate the test average for each dataset
test.average.ridge<- auroc.all.ridge%>% 
  filter(study.train != study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC = mean(AUC), .groups = "drop") %>% 
  mutate(study.train = "Average")

# combine AUROC values with test averages
bind_rows(auroc.all.ridge, test.average.ridge) %>% 
  # highlight cross validation versus transfer results
  mutate(CV = study.train == study.test) %>%
  # for faceting
  mutate(split = case_when(study.train == "Average" ~ "Average", TRUE ~ "none")) %>% 
  mutate(split = factor(split, levels = c("none", "Average"))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train = factor(study.train, levels = c(datasets, "Average"))) %>% 
  mutate(study.test = factor(study.test, levels = c(rev(datasets), "Average"))) %>% 
  ggplot(aes(y = study.test, x = study.train, fill = AUC, size = CV, color = CV)) +
  geom_tile() + theme_minimal() +
  geom_text(aes_string(label = "format(AUC, digits = 2)"), col = "white", size = 2) +
  scale_fill_gradientn(colours = rev(c("darkgreen","forestgreen", "chartreuse3", "lawngreen", "yellow"))) +
  scale_x_discrete(position = "top") + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle = 45, hjust = .1), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank()) + 
  xlab("Training Set") + ylab("Test Set") + 
  scale_color_manual(values = c("#FFFFFF00", "grey"), guide = "none") + 
  scale_size_manual(values = c(0, 1), guide = "none") + 
  facet_grid(~split, scales = "free", space = "free")


### INVESTIGATE FEATURE WEIGHTS
# trained models are saved in sc.list
# extract the models weights and compare them to the associations calculated above

weight.list.ridge<- list()

for (d in datasets){
  sc.obj.train <- sc.list[[d]]
  # extract the feature weights out of the SIAMCAT object
  temp <- feature_weights(sc.obj.train)
  temp$genus <- rownames(temp)
  # save selected info in the weight.list.ridge
  weight.list.ridge[[d]] <- temp %>% 
    select(genus, median.rel.weight, mean.rel.weight, percentage) %>% 
    mutate(Study = d) %>% 
    mutate(r.med = rank(-abs(median.rel.weight)), 
           r.mean = rank(-abs(mean.rel.weight)))
}

# combine all feature weights into a single tibble
df.weights.ridge<- bind_rows(weight.list.ridge)
sum(df.weights.ridge$genus == "unclassified")
df.weights.ridge<- df.weights.ridge%>% filter(genus != "unclassified")


### plot a heatmap with the weights for the genera of interest (those that had their associations plotted above)
# compute absolute feature weights
abs.weights.ridge<- df.weights.ridge%>% 
  group_by(Study) %>% 
  summarise(sum.median = sum(abs(median.rel.weight)),
            sum.mean = sum(abs(mean.rel.weight)),
            .groups = "drop")

df.weights.ridge%>% 
  full_join(abs.weights.ridge) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight = median.rel.weight/sum.median) %>% 
  # only include genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med = case_when(r.med > 20~NA_real_, TRUE ~ r.med)) %>%
  # enforce the correct ordering by converting to factors 
  mutate(genus = factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  mutate(Study = factor(Study, levels = datasets)) %>% 
  ggplot(aes(y = genus, x = Study, fill = median.rel.weight)) + geom_tile() + 
  scale_fill_gradientn(colours=rev(c("#007A53", "#009F4D", "#6CC24A", "white",
                                     "#EFC06E", "#FFA300", "#BE5400")), 
                       limits = c(-0.15, 0.15)) + theme_minimal() + 
  geom_text(aes(label = r.med), col = "black", size = 2) +
  theme(panel.grid = element_blank()) + xlab("") + ylab("") +
  theme(axis.text = element_text(size = 6))


#############################################################
### MACHINE LEARNING META-ANALYSIS - RANDOM FOREST MODELS ###
#############################################################


# train one model for each of the datasets, then apply it to the other datasets (holdout testing functionality of SIAMCAT)
# block the cross-validation of subjects with repeated samples to not bias the results

# create tibble to store all of the predictions
auroc.all.rf<- tibble(study.train = character(0), 
                      study.test = character(0),
                      AUC = double(0))

# create a list to save the trained SIAMCAT objects
sc.list.rf<- list()

for (i in datasets){
  # restrict to a single study
  meta.train <- meta.all %>% 
    filter(Study == i) %>% 
    as.data.frame()
  rownames(meta.train) <- meta.train$Sample_ID
  
  ## take into account repeated sampling by blocking CV by Individual_ID (for studies with repeated samples)
  block <- NULL
  if (i %in% c("metaHIT", "Lewis_2015", "HMP2")){
    block <- "Individual_ID"
    if (i == "HMP2"){ 
      # for HMP2 dataset, reduce number of repeated samples per subject (some have up to 20 samples) 
      meta.train <- meta.all %>% 
        filter(Study == "HMP2") %>% 
        group_by(Individual_ID) %>% 
        sample_n(5, replace = TRUE) %>% 
        distinct() %>% 
        as.data.frame()
      rownames(meta.train) <- meta.train$Sample_ID
    }
  }
  # create SIAMCAT object
  sc.obj.train <- siamcat(feat = feat, meta = meta.train, 
                          label = "Group", case = "CD")
  # normalize features
  sc.obj.train <- normalize.features(sc.obj.train, norm.method = "log.std",
                                     norm.param = list(log.n0 = 1e-05, sd.min.q = 0),
                                     feature.type = "original")
  # split data
  sc.obj.train <- create.data.split(sc.obj.train,
                                    num.folds = 10, num.resample = 10, inseparable = block)
  # train RANDOM FOREST model
  sc.obj.train <- train.model(sc.obj.train, method = "enet")
  
  ## apply trained models to other datasets
  # loop through datasets again
  for (i2 in datasets){
    if (i == i2){
      # make and evaluate cross-validation predictions (same dataset)
      sc.obj.train <- make.predictions(sc.obj.train)
      sc.obj.train <- evaluate.predictions(sc.obj.train)
      auroc.all.rf<- auroc.all.rf%>% 
        add_row(study.train = i, study.test = i,
                AUC = eval_data(sc.obj.train)$auroc %>% as.double())
    } else {
      # make and evaluate cross-validation predictions (external datasets)
      # use meta.ind here, since we want only one sample per subject
      meta.test <- meta.ind %>% 
        filter(Study == i2) %>%
        as.data.frame()
      rownames(meta.test) <- meta.test$Sample_ID
      sc.obj.test <- siamcat(feat = feat, meta = meta.test,
                             label = "Group", case = "CD")
      # make holdout predictions
      sc.obj.test <- make.predictions(sc.obj.train, 
                                      siamcat.holdout = sc.obj.test)
      sc.obj.test <- evaluate.predictions(sc.obj.test)
      auroc.all.rf<- auroc.all.rf%>% 
        add_row(study.train = i, study.test = i2,
                AUC = eval_data(sc.obj.test)$auroc %>% as.double())
    }
  }
  # save the trained model
  sc.list.rf[[i]] <- sc.obj.train
}

# calculate the test average for each dataset
test.average.rf<- auroc.all.rf%>% 
  filter(study.train != study.test) %>% 
  group_by(study.test) %>% 
  summarise(AUC = mean(AUC), .groups = "drop") %>% 
  mutate(study.train = "Average")

# combine AUROC values with test averages
bind_rows(auroc.all.rf, test.average.rf) %>% 
  # highlight cross validation versus transfer results
  mutate(CV = study.train == study.test) %>%
  # for faceting
  mutate(split = case_when(study.train == "Average" ~ "Average", TRUE ~ "none")) %>% 
  mutate(split = factor(split, levels = c("none", "Average"))) %>% 
  # convert to factor to enforce ordering
  mutate(study.train = factor(study.train, levels = c(datasets, "Average"))) %>% 
  mutate(study.test = factor(study.test, levels = c(rev(datasets), "Average"))) %>% 
  ggplot(aes(y = study.test, x = study.train, fill = AUC, size = CV, color = CV)) +
  geom_tile() + theme_minimal() +
  geom_text(aes_string(label = "format(AUC, digits = 2)"), col = "white", size = 2) +
  scale_fill_gradientn(colours = rev(c("darkgreen","forestgreen", "chartreuse3", "lawngreen", "yellow"))) +
  scale_x_discrete(position = "top") + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle = 45, hjust = .1), 
        panel.grid = element_blank(), panel.border = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank()) + 
  xlab("Training Set") + ylab("Test Set") + 
  scale_color_manual(values = c("#FFFFFF00", "grey"), guide = "none") + 
  scale_size_manual(values = c(0, 1), guide = "none") + 
  facet_grid(~split, scales = "free", space = "free")


### INVESTIGATE FEATURE WEIGHTS
# trained models are saved in sc.list
# extract the models weights and compare them to the associations calculated above

weight.list.rf<- list()

for (d in datasets){
  sc.obj.train <- sc.list[[d]]
  # extract the feature weights out of the SIAMCAT object
  temp <- feature_weights(sc.obj.train)
  temp$genus <- rownames(temp)
  # save selected info in the weight.list.rf
  weight.list.rf[[d]] <- temp %>% 
    select(genus, median.rel.weight, mean.rel.weight, percentage) %>% 
    mutate(Study = d) %>% 
    mutate(r.med = rank(-abs(median.rel.weight)), 
           r.mean = rank(-abs(mean.rel.weight)))
}

# combine all feature weights into a single tibble
df.weights.rf<- bind_rows(weight.list.rf)
sum(df.weights.rf$genus == "unclassified")
df.weights.rf<- df.weights.rf%>% filter(genus != "unclassified")


### plot a heatmap with the weights for the genera of interest (those that had their associations plotted above)
# compute absolute feature weights
abs.weights.rf<- df.weights.rf%>% 
  group_by(Study) %>% 
  summarise(sum.median = sum(abs(median.rel.weight)),
            sum.mean = sum(abs(mean.rel.weight)),
            .groups = "drop")

df.weights.rf%>% 
  full_join(abs.weights.rf) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight = median.rel.weight/sum.median) %>% 
  # only include genera of interest
  filter(genus %in% genera.of.interest$genus) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med = case_when(r.med > 20~NA_real_, TRUE ~ r.med)) %>%
  # enforce the correct ordering by converting to factors 
  mutate(genus = factor(genus, levels = rev(genera.of.interest$genus))) %>% 
  mutate(Study = factor(Study, levels = datasets)) %>% 
  ggplot(aes(y = genus, x = Study, fill = median.rel.weight)) + geom_tile() + 
  scale_fill_gradientn(colours=rev(c("#007A53", "#009F4D", "#6CC24A", "white",
                                     "#EFC06E", "#FFA300", "#BE5400")), 
                       limits = c(-0.15, 0.15)) + theme_minimal() + 
  geom_text(aes(label = r.med), col = "black", size = 2) +
  theme(panel.grid = element_blank()) + xlab("") + ylab("") +
  theme(axis.text = element_text(size = 6))


### RESULTS OF MACHINE LEARNING MODELS

combined_test.avg <- data.frame(
  study.test = test.average.lasso$study.test,
  AUC.lasso = test.average.lasso$AUC,
  AUC.enet = test.average.enet$AUC,
  AUC.ridge = test.average.ridge$AUC,
  AUC.rf = test.average.rf$AUC
)

#      study.test AUC.lasso  AUC.enet AUC.ridge    AUC.rf
# 1 Franzosa_2019 0.8098620 0.8491274 0.7980418 0.8638393
# 2          HMP2 0.6598077 0.6921154 0.6928846 0.6926923
# 3       He_2017 0.8926646 0.8978629 0.8787062 0.8948787
# 4    Lewis_2015 0.8040000 0.8384706 0.8157647 0.8381176
# 5       metaHIT 0.8846154 0.9155802 0.8797262 0.9260104


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
#   [1] future_1.58.0   lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4    
# [7] readr_2.1.5     tidyr_1.3.1     tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0 SIAMCAT_2.12.0 
# [13] phyloseq_1.52.0 mlr3_1.0.0     
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            beanplot_1.3.1          pROC_1.18.5             gridExtra_2.3          
# [5] permute_0.9-7           rlang_1.1.6             magrittr_2.0.3          gridBase_0.4-7         
# [9] ade4_1.7-23             matrixStats_1.5.0       compiler_4.5.0          mgcv_1.9-3             
# [13] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         shape_1.4.6.1          
# [17] crayon_1.5.3            backports_1.5.0         XVector_0.48.0          labeling_0.4.3         
# [21] PRROC_1.4               utf8_1.2.6              tzdb_0.5.0              UCSC.utils_1.4.0       
# [25] nloptr_2.2.1            bit_4.6.0               glmnet_4.1-9            mlr3misc_0.18.0        
# [29] GenomeInfoDb_1.44.0     jsonlite_2.0.0          biomformat_1.36.0       progress_1.2.3         
# [33] rhdf5filters_1.20.0     uuid_1.2-1              Rhdf5lib_1.30.0         mlr3measures_1.0.0     
# [37] parallel_4.5.0          prettyunits_1.2.0       cluster_2.1.8.1         R6_2.6.1               
# [41] stringi_1.8.7           RColorBrewer_1.1-3      parallelly_1.45.0       boot_1.3-31            
# [45] numDeriv_2016.8-1.1     Rcpp_1.0.14             iterators_1.0.14        future.apply_1.20.0    
# [49] IRanges_2.42.0          timechange_0.3.0        Matrix_1.7-3            splines_4.5.0          
# [53] igraph_2.1.4            tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1      
# [57] mlr3tuning_1.4.0        vegan_2.7-1             codetools_0.2-20        curl_6.3.0             
# [61] listenv_0.9.1           lattice_0.22-7          lmerTest_3.1-3          plyr_1.8.9             
# [65] Biobase_2.68.0          withr_3.0.2             survival_3.8-3          Biostrings_2.76.0      
# [69] infotheo_1.2.0.1        pillar_1.10.2           corrplot_0.95           checkmate_2.3.2        
# [73] foreach_1.5.2           stats4_4.5.0            reformulas_0.4.1        generics_0.1.4         
# [77] vroom_1.6.5             bbotk_1.6.0             S4Vectors_0.46.0        hms_1.1.3              
# [81] scales_1.4.0            minqa_1.2.8             globals_0.18.0          glue_1.8.0             
# [85] LiblineaR_2.10-24       tools_4.5.0             data.table_1.17.4       lme4_1.1-37            
# [89] rhdf5_2.52.1            grid_4.5.0              ape_5.8-1               rbibutils_2.3          
# [93] colorspace_2.1-1        paradox_1.0.1           nlme_3.1-168            GenomeInfoDbData_1.2.14
# [97] palmerpenguins_0.1.1    cli_3.6.5               gtable_0.3.6            digest_0.6.37          
# [101] BiocGenerics_0.54.0     lgr_0.4.4               farver_2.1.2            multtest_2.64.0        
# [105] lifecycle_1.0.4         mlr3learners_0.12.0     httr_1.4.7              bit64_4.6.0-1          
# [109] MASS_7.3-65

