# Machine Learning
## :deciduous_tree::evergreen_tree::deciduous_tree: Random Forest
<!-- A random forest averages multiple decision trees that were trained on different parts of the same training set in order to overcome the over-ftting observed with individual decision trees. -->
<!-- overfitting: explaining your training data instead of finding generalizable patterns -->
<!-- classification: categorical dependent variable and regression: continuous dependent variable -->
<!-- each tree is trained on 2/3rd of the total training data - cases drawn at random with replacement -->
<!-- m (default in the square root of the total number of all predictors for classification | for regression, m is the total number of all predictors divided by 3) predictor variables are selected at random out of all the predictor variables and the best split on these m is used to split the node -->
<!-- for each tree, the leftover 1/3 of the data is used to calculate the misclassification rate : out of bag (OOB) error rate -- aggregate error from all trees is used to determine overall OOB error rate for the classification -- eah tree gives a classification on leftover data (OOB): the tree votes for that class
<!-- imbalanced data set: if one class contains significantly more samples than the other (makes it hard to create appropriate testing and training data sets - most classifiers are built with the assumption that hte test data is dran from the same distribution as the training data) -->
<!-- two parameters that are important in the random forest algorithm are the number of trees (mtree) used in the forest and the number of random varaibles used in each tree (mtry): where does the OOB error rate reach a minimum and stabilize -->
<!-- for datasets where there is a higher number of features present, a good idea is to use cross-validation to perform feature selection using the OOB error rate -->


<!-- we can measure how accurate the random forest model is by the proportion of out-of-bag samples that were correctly classified by the random forest - the proportion of samples that were incorrectly classified is the out of bag error -->

<!-- Create Many Decision Trees: The algorithm makes many decision trees each using a random part of the data. So every tree is a bit different. -- Pick Random Features: When building each tree it doesn’t look at all the features (columns) at once. It picks a few at random to decide how to split the data. This helps the trees stay different from each other. -- Each Tree Makes a Prediction: Every tree gives its own answer or prediction based on what it learned from its part of the data. -- Combine the Predictions: --For classification we choose a category as the final answer is the one that most trees agree on i.e majority voting.
For regression we predict a number as the final answer is the average of all the trees predictions. -- Why It Works Well: Using random data and features for each tree helps avoid overfitting and makes the overall prediction more accurate and trustworthy. -->




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

