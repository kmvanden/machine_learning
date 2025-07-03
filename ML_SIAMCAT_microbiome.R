# Machine Learning Workflow of Microbiome Data using SIAMCAT

# load libraries
library(tidyverse)
library(ggplot2)
library(SIAMCAT)
library(ranger)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### load data
# sample names need to be rownames in the metadata
meta <- read.table("metadata.txt", header = TRUE)
head(meta)
rownames(meta) <- meta$sample_name
meta$sample_name <- NULL

# feature table needs to be in relative abundances
feat <- read.table("feature_table.txt", header = TRUE)
feat[1:5, 1:5]
feat <- apply(feat, 2, function(x) x / sum(x))
feat["Faecalibacterium_prausnitzii", 1:5]

# rownames of metadata need to match the column names of the feature table
all(rownames(meta) == colnames(feat))

# create label object from condition column in meta
label <- create.label(meta = meta, label = "condition", case = "disease")


### create SIAMCAT object
obj <- siamcat(feat = feat, label = label, meta = meta)
show(obj)
dim(obj@phyloseq@otu_table) # 2302   70

### perform unsupervised feature selection
# remove features whose maximum abundance is never above 0.00001
obj <- filter.features(obj, filter.method = "abundance", cutoff = 0.00001)
# remove features that are not detected in at least 20% of the samples
obj <- filter.features(obj, cutoff = 0.20, filter.method = "prevalence", 
                       feature.type = "filtered")
dim(obj@filt_feat$filt.feat) # 681  70


### calculate associations between features and the label
obj <- check.associations(obj, alpha = 0.1, log.n0 = 1e-07)
show(obj)
association.plot(obj, sort.by = "fc", panels = c("fc", "prevalence", "auroc"))
# association.plot(obj, plot.type = "quantile.rect", fn.plot = "./association_plot.pdf")

assoc <- associations(obj) # extract associations
assoc <- assoc %>%
  arrange(p.adj) # arrange by p.adj

# top 20 species by padj
top <- rownames(assoc)[1:20]

### check for confounders
check.confounders(obj, meta.in = NULL, feature.type = "filtered", fn.plot = "confounder_plots.pdf")

### factorize.metadata()
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

# get label, meta and feat
label <- label(obj)
meta <- meta(obj)
feat <- get.filt_feat.matrix(obj)

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

### get variance explained by sex
sex <- meta$sex
sex <- factorize.metadata(data.frame(sex = sex), verbose = 0)$sex
names(sex) <- rownames(meta)

# ensure names stil match
all(names(sex) %in% colnames(feat))
all(colnames(feat) %in% names(sex))

var.sex <- vapply(rownames(feat), FUN = function(x){
  x_vals <- rank(feat[x, names(sex)]) / length(sex)
  ss.tot <- sum((x_vals - mean(x_vals))^2) / length(x_vals)
  ss.o.i <- sum(vapply(levels(sex), function(s){
    sum((x_vals[sex == s] - mean(x_vals[sex == s]))^2)
  }, FUN.VALUE = double(1))) / length(x_vals)
  return(1 - ss.o.i / ss.tot)
}, FUN.VALUE = double(1))


# plot size as average log abundance values
r.mean   <- rowMeans(log10(feat + 1e-05))
mean.size <- (r.mean + 5) * 8/5 + 1

# define maximum values for axes
lim <- max(c(var.label, var.sex), na.rm = TRUE) 

plot(var.label, var.sex, type = "n",
     xlab = "Variance explained by label",
     ylab = "Variance expained by sex",
     xlim = c(0, lim), ylim = c(0, lim))
symbols(x = var.label, y = var.sex, circles = mean.size, inches = 1/9,
        bg = alpha("darkgrey", 0.4), fg = alpha("black", 0.7), add = TRUE)
abline(0, 1, lty = 3, col = "black")

# sort variance explained by label
sort_var <- sort(var.label, decreasing = TRUE)
sort_var[1:10]

# sort variance explained by sex
sort_sex <- sort(var.sex, decreasing = TRUE)
sort_sex[1:10]

# data.frame with variance explained by label and sex
variance <- data.frame(
  sort_var = sort_var[names(sort_var)],
  sort_sex = sort_sex[names(sort_var)]
)
variance["Eubacterium_sp._c-25", ]


# highlight specific species on the above variance plot
point.names <- names(var.label) # vector of names
target <- c("Eubacterium_sp._c-25", "Lachnoclostridium_sp._YL32", "Faecalicatena_sp._Marseille-Q4148") # names to highlight
highlight.index <- match(target, names(var.label))

bg.colors <- rep(alpha("darkgrey", 0.4), length(var.label)) # fill color
fg.colors <- rep(alpha("black", 0.7), length(var.label)) # border color

# length needs to be the same as target length
bg.colors[highlight.index] <- alpha(c("red", "blue", "purple"), 0.7) # fill color for highlighted point
fg.colors[highlight.index] <- c("red3", "blue3", "purple3") # border color for highlighted point

# plot symbols with highlighted species
symbols(x = var.label, y = var.sex, circles = mean.size, inches = 1/9,
        bg = bg.colors, fg = fg.colors, add = TRUE)


##### pre-process data for machine learning models
### log normalize features (pseudocount = 1e-06) | normalize over all features
obj <- normalize.features(obj, norm.method = "log.unit",
                          norm.param = list(log.n0 = 1e-06, n.p = 2, norm.margin = 1))

### split data into 5 cross-validation folds with 50 resampling rounds
obj <- create.data.split(obj, num.folds = 5, num.resample = 50)


############# train a lasso model on training folds
obj.lasso <- train.model(obj, method = "lasso")
show(obj.lasso)

### apply model to testing folds
obj.lasso <- make.predictions(obj.lasso)
pred_matrix.lasso <- pred_matrix(obj.lasso)
head(pred_matrix.lasso)

### evaluate predictions
obj.lasso <-  evaluate.predictions(obj.lasso)
model.evaluation.plot(obj.lasso, show.all = TRUE) 
# AUROC: 0.638
# precision-recall AUC: 0.728

# consens.thres = minimal ratio of models incorporating a feature in order to include it into the heatmap
model.interpretation.plot(obj.lasso, fn.plot = "interpretation_lasso.pdf",
                          consens.thres = 0.5, limits = c(-3, 3), 
                          heatmap.type = "zscore") # interpretation plot


############# train a random forest model on training folds
parallel::detectCores()
options(ranger.num.threads = 8)
obj.rf <- train.model(obj, method = "randomForest")
show(obj.rf)

### apply model to testing folds
obj.rf <- make.predictions(obj.rf)
pred_matrix.rf <- pred_matrix(obj.rf)
head(pred_matrix.rf)

### evaluate predictions
obj.rf <-  evaluate.predictions(obj.rf)
model.evaluation.plot(obj.rf, show.all = TRUE) 

# AUROC: 0.631
# precision-recall AUC: 0.679

# consens.thres = for random forest, specifies the minimum median Gini coefficient for a feature to be included
model.interpretation.plot(obj.rf, fn.plot = "interpretation_rf.pdf",
                          consens.thres = 0.01, limits = c(-3, 3), 
                          heatmap.type = "zscore") # interpretation plot


############# train a ridge model on training folds
obj.ridge <- train.model(obj, method = "ridge")
show(obj.ridge)

### apply model to testing folds
obj.ridge <- make.predictions(obj.ridge)
pred_matrix.ridge <- pred_matrix(obj.ridge)
head(pred_matrix.ridge)

### evaluate predictions
obj.ridge <-  evaluate.predictions(obj.ridge)
model.evaluation.plot(obj.ridge, show.all = TRUE) 
# AUROC: 0.643
# precision-recall AUC: 0.696

# consens.thres = minimal ratio of models incorporating a feature in order to include it into the heatmap
model.interpretation.plot(obj.ridge, fn.plot = "interpretation_ridge.pdf",
                          consens.thres = 0.5, limits = c(-3, 3), 
                          heatmap.type = "zscore") # interpretation plot


########## train a elastic net model on training folds
obj.enet <- train.model(obj, method = "enet")
show(obj.enet)

### apply model to testing folds
obj.enet <- make.predictions(obj.enet)
pred_matrix.enet <- pred_matrix(obj.enet)
head(pred_matrix.enet)

### evaluate predictions
obj.enet <- evaluate.predictions(obj.enet)
model.evaluation.plot(obj.enet, show.all = TRUE) 
# AUROC: 0.618
# precision-recall AUC: 0.706

# consens.thres = minimal ratio of models incorporating a feature in order to include it into the heatmap
model.interpretation.plot(obj.enet, fn.plot = "interpretation_enet.pdf",
                          consens.thres = 0.5, limits = c(-3, 3), 
                          heatmap.type = "zscore") # interpretation plot


### investigate feature weights
# list of SIAMCAT objects
siamcat_list <- list(obj.lasso, obj.rf, obj.enet, obj.ridge)
model_names <- c("lasso", "rf", "enet", "ridge")

# extract selected feature weights from the objects
featWeights_list <- lapply(seq_along(siamcat_list), function(i) {
  obj <- siamcat_list[[i]]
  temp <- feature_weights(obj)
  temp$species <- rownames(temp)
  
  temp %>%
    select(species, median.rel.weight, mean.rel.weight, percentage) %>%
    mutate(model = model_names[i]) %>%
    mutate(
      r.med = rank(-abs(median.rel.weight)),
      r.mean = rank(-abs(mean.rel.weight))
    )
})


# absolute median relative weight across lasso models
feat.weights.lasso <- featWeights_list[[1]]
feat.weights.lasso <- feat.weights.lasso %>%
  arrange(r.med)
head(feat.weights.lasso, 15)

# absolute median relative weight across random forest models
feat.weights.rf <- featWeights_list[[2]]
feat.weights.rf <- feat.weights.rf %>%
  arrange(r.med)
head(feat.weights.rf, 15)

# absolute median relative weight across enet models
feat.weights.enet <- featWeights_list[[3]]
feat.weights.enet <- feat.weights.enet %>%
  arrange(r.med)
head(feat.weights.enet, 15)

# absolute median relative weight across ridge models
feat.weights.ridge <- featWeights_list[[4]]
feat.weights.ridge <- feat.weights.ridge %>%
  arrange(r.med)
head(feat.weights.ridge, 15)


# combine all feature weights into a single tibble
feat.weights.all <- bind_rows(featWeights_list)
table(feat.weights.all$model)

# compute absolute feature weights
feat.weights.abs <- feat.weights.all %>% 
  group_by(model) %>% 
  summarise(sum.median = sum(abs(median.rel.weight)),
            sum.mean = sum(abs(mean.rel.weight)),
            .groups = "drop")

# plot heatmap of feature weights
feat.weights.all %>% 
  full_join(feat.weights.abs) %>% 
  # normalize by the absolute model size
  mutate(median.rel.weight = median.rel.weight/sum.median) %>% 
  # only include top 20 species by associations (padj)
  filter(species %in% top) %>% 
  # highlight feature rank for the top 20 features
  mutate(r.med = case_when(r.med > 20~NA_real_, TRUE ~ r.med)) %>%
  # enforce the correct ordering by converting to factors 
  mutate(species = factor(species, levels = rev(top))) %>% 
  mutate(model = factor(model, levels = model_names)) %>% 
  ggplot(aes(y = species, x = model, fill = median.rel.weight)) + geom_tile() + 
  scale_fill_gradientn(colours = rev(c("#007A53", "#009F4D", "#6CC24A", "white",
                                     "#EFC06E", "#FFA300", "#BE5400")), 
                       limits = c(-0.15, 0.15)) + theme_minimal() + 
  geom_text(aes(label = r.med), col = "black", size = 2) +
  theme(panel.grid = element_blank()) + xlab("") + ylab("") +
  theme(axis.text = element_text(size = 6))


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
#   [1] future_1.58.0   ranger_0.17.0   SIAMCAT_2.12.0  phyloseq_1.52.0 mlr3_1.0.0      lubridate_1.9.4
# [7] forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
# [13] tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            beanplot_1.3.1          pROC_1.18.5             gridExtra_2.3          
# [5] permute_0.9-7           rlang_1.1.6             magrittr_2.0.3          gridBase_0.4-7         
# [9] ade4_1.7-23             matrixStats_1.5.0       compiler_4.5.0          mgcv_1.9-3             
# [13] vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         shape_1.4.6.1          
# [17] crayon_1.5.3            backports_1.5.0         XVector_0.48.0          labeling_0.4.3         
# [21] PRROC_1.4               utf8_1.2.6              tzdb_0.5.0              UCSC.utils_1.4.0       
# [25] nloptr_2.2.1            bit_4.6.0               glmnet_4.1-9            mlr3misc_0.18.0        
# [29] GenomeInfoDb_1.44.0     jsonlite_2.0.0          progress_1.2.3          biomformat_1.36.0      
# [33] rhdf5filters_1.20.0     uuid_1.2-1              Rhdf5lib_1.30.0         mlr3measures_1.0.0     
# [37] prettyunits_1.2.0       parallel_4.5.0          cluster_2.1.8.1         R6_2.6.1               
# [41] stringi_1.8.7           RColorBrewer_1.1-3      parallelly_1.45.0       boot_1.3-31            
# [45] numDeriv_2016.8-1.1     Rcpp_1.0.14             iterators_1.0.14        future.apply_1.20.0    
# [49] IRanges_2.42.0          Matrix_1.7-3            splines_4.5.0           igraph_2.1.4           
# [53] timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1      
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

