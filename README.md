# Machine Learning

<!-- it is important to estimate how well a trained model generalizes to an independent test set. this is typically the objective of microbial biomarker discovery -->

<!-- common pitfalls leading to poor generalization of machine learning models -->
<!-- machine learning workflows that are incorrectly set-up can lead overoptomistic accuracy estimates (overfitting) -->

<!-- 1. naive combination of feature selection on the whole dataset and subsequent cross-validation on the same data --> 
<!-- 2. arises when samples not taken independently (not replicates or at multiple time points from the same subject) are randomly partitioned in CV with the aim to assess the cross-subject generalization error -->
<!-- external validation shows that overfitting occurs when feature selection and cross-validation are incorrectly combined in a sequential manner, rather than in a nested approach -- the fewer the features selected, the more pornounced the issue becomes -->
<!-- supervised feature selection should always be nested into CV: supervised feature selection has to be applied to each training fold of the CV seperately -->

<!-- when dependent observations are randomly assigned to CV partitions, the ability of the model to generalize across subject is not assessed (generalized across timepoints) -- repeated measurements need to be blocked (either all of them into the training set or all of them into the test set) -->

<!-- 1. data preprocessesing: unsupervised abundance and prevalance filtering -->
<!-- 2. univariate associations of single species with the disease using the non-parametric Wilcoxon test (shown to reliably control the FDR rate in metagenomics data) -->
<!-- 3. normalization, set-up CV scheme, ML algorithms (LASSO, Elastic Net and RF) -- models trained and applied to test data -- performance of model is assessed using AUROC and interprability plots of the importance of individual features in the classification model are presented -->
