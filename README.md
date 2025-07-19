# Machine Learning
## :evergreen_tree::deciduous_tree::evergreen_tree: Random Forest
Random forests are an ensemble learning method that aggregates the predictions of multiple decision trees built using bootstrap sampling and random feature selection in order to improve classification (categorical outcomes) or regression (continuous outcomes) performance.

### Overview of Random Forest Model Construction
Data is split into training data and testing data. 
  - **Training data**: the portion of the data used to train (build) the models. The models learn the patterns/relationships between the input features and the target labels.
  - **Testing data**: the portion of the data that is used to evaluate the models (determine how well the models generalizes to unseen data).
  - **Cross-validation**: resampling technique used to evaluate the performance of a model on unseen data and avoid overfitting (by estimating how well the model might perform on unseen data).
    - **Overfitting**: when a model performs well on training data, but poorly on testing data (the model has learned noise or specific patterns only present in the training set that are not generalizable to other datasets).
    - **k-fold cross-validation**: splits the data into k equal parts (folds), trains the model on k-1 folds and tests the model on the remaining fold. The data used in the folds is then rotated, and the process is repeated k times so that each fold is used once as the test dataset. The k-fold process can be repeated multiple times with different fold splits for more reliable estimates.
    - **Stratification**: ensures that the class distribution in each fold of the cross-validation split reflects the overall distribution of the dataset so that each training and testing fold has a representative mix of classes.
    
Multiple subsets of the training data are created using bootstrap sampling and each subset is used to build an individual tree. 
  - **Bootstrap sampling**: generates multiple new datasets from the training data using random sampling with replacement. Each bootstrap sample is the same size as the original dataset, but some of the samples are repeated and some are omitted. 

At each node split in a decision tree, only a random subset of features is considered. 
  - **Random feature selection**: each decision tree is trained on a bootstrap sample, but at each node split, only a random subset of the features (```mtry```) is selected, thereby decreasing correlation between the trees and improving performance.

The samples in the testing data are then dropped down the tree until they reach a terminal node and then are assigned the label of the majority class in that node (classification) or the mean of the target variable (regression). The predictions made by each tree are then aggregated by majority vote (classification) or average prediction (regression) to make the overall prediction.
  - **Majority vote**: each model votes for a class label and the class with the most votes wins.
  - **Average prediction**: each model outputs a numeric value and the average number is the final prediction.

### Random Forest Hyperparameters
Hyperparameters are tunable settings that are set before the training process begins. They control how the model is trained and thus the performance of the model.

Each decision tree within a random forest is trained on a bootstrap sample of the training data. This introduces variation, so that each tree doesn’t learn the same patterns. On average, about two-thirds of the training data is included in any given boostrap sample (the remaining one-third of the data are the out-of-bag (OOB) samples).

Features within the data (especially in high dimensional datasets) can be irrelevant to class separation. Models trained with many irrelevant features are more likely to learn patterns specific to the training dataset and thus not generalize well to unseen data. 
  - **Feature selection**: the process of identifying and using only the most informative features in the dataset for training a model.
    - **Recursive feature elimination**: a model is iteratively trained, removing the least important features (based on mean decrease in Gini impurity or mean decrease in accuracy) at each iteration and comparing model performance for the different numbers of features used to train the model.
    - **Mean decrease in Gini impurity**: measures how much a feature contributes to decreasing Gini impurity (how pure or homogeneous a set of classes is at a given node in a decision tree) when it is used to split data across all of the trees in the forest.
      - High mean decrease in Gini impurity = splitting the data by the feature results in nodes with increased purity
    - **Mean decrease in accuracy**: measures how much model accuracy drops when the feature is randomly permuted. Values of the feature in the OOB samples are permuted and the accuracy of the model before and after the permutation is measured.
      - High mean decrease in accuracy = the feature contributes a lot to the accuracy of the predictions made by the model
    - Both mean decrease in Gini and mean decrease in accuracy can be biased toward continuous features and categorical features with many levels (more unique values = more opportunities to split the data = more chances to fit the target or reduce impurity, even by chance) and highly correlated features (they share predictive information, so the decrease in accuracy/impurity is masked).

Each tree is grown by recursively splitting the data at internal decision nodes, starting at the root node. With random forests, only a random subset of features is considered at each split to help decorrelate the trees and thereby increase the robustness of the ensemble.
  - **Number of features per split** (```mtry```): number of features randomly selected at each
node split. The default is the square root of the number of features for classification tasks and one third of the features for regression tasks. 
    - Lower ```mtry``` = more randomization. 
    - Higher ```mtry``` = potentially better split quality, but higher correlation between the trees
   
At each split, among the selected mtry features, the algorithm chooses the feature and the threshold that best separates the classes (how homogenous the resulting nodes are after the split). For classification, the goal is to reduce impurity (which is typically measured using Gini impurity), and for regression the goal is to minimize variance in the child nodes. However, if the classes are imbalanced, the model will favor majority class predictions and the model will have poor sensitivity for the minority class. 
  - **Imbalanced dataset**: target classes (labels) are not represented equally within the dataset. Standard models tend to optimize for overall accuracy, which favors the majority class, resulting in the minority class being underpredicted.
  - **Class weights** (```classwt```): correct the imbalance between class labels by adjusting the penalty for misclassifying samples of different classes. Essentially, class weights make the model pay more attention to the minority class by increasing the impact (effective contribution) of the minority class on impurity calculations and thereby influencing which splits are considered good.
    - If the class weights are too extreme, the model may overfit the minority class.

The splitting process continues recursively for each new node until a stopping rule is met. Typical stopping conditions include: all samples in a node belonging to the same class or the node containing fewer samples than ```nodesize```.
  - **Minimum node size for a split** (```nodesize```): minimum number of samples required to split an internal node. ```nodesize``` controls tree depth and complexity.
One is the default for classification tasks and five is the default for regression tasks. 
    - Smaller values = deeper trees with more splits that capture more detail, but risk overfitting the data (especially in noisy or small datasets).
    - Larger values = shallower trees with fewer splits that capture less detail and may underfit the model (miss patterns in the data).

The decision tree-building process is then repeated ```ntree``` times to create an ensemble of decision trees, each built on a different bootstrap sample of the training data, and each making independent predictions that are then aggregated.
  - **Number of trees** (```ntree```): determines how many decision trees are in the forest. 
More trees = more stable and consistent predictions, improvement of generalization and accuracy (up to a point), but requires longer computation time

## :cowboy_hat_face::mountain::spider_web: Logistic Regression
Logistic regression is a statistical classification algorithm that models the relationship between predictor variables and outcome labels.

### Overview of Logisitic Regression Model Construction
The input features are combined linearly and coefficients (weights) are applied to each predictor variable, creating the linear predictor.
  - Linear predictor = β0​+β1​x1​+β2​x2​+⋯+βn​xn​ (for a binary logistic regression model)
  - Coefficients are typically initialized at zero (or at very small random values), which leads to initial predicted probabilities near 0.5 (unbiased starting point).

The logistic regression model assumes that the log-odds of the probability of the outcome being in the positive class is a linear combination of the input features.
  - Linear predictor = β0​+β1​x1​+β2​x2​+⋯+βn​xn​ = log(p/1-p​) = log odds
  - Log-odds are used as an intermediary rather than directly using probabilities because the output of a linear combination (like the linear predictor) is unbounded, whereas probabilities have to be between 0 and 1. 

To convert the log-odds into a probabilities, the log-odds transformation is inverted using the logistic (sigmoid) function, which maps any real number to a value between 0 and 1, producing an S-shaped curve.
  - Probability = 1/(1+ e^-(β0​+β1​x1​+⋯+βn​xn​))
  - Small changes in the log-odds around zero result in rapid changes in probability near 0.5, allowing smooth transitions between class probabilities, and the curve flattens at extreme values, reflecting confidence in the classification.

Optimization algorithms, such as gradient descent or coordinate descent, are used to learn the coefficients that minimize the loss function, aiming to produce predicted probabilities that match the true class labels as closely as possible.
  - **Loss function**: a mathematical expression that measures how far off the predictions of the model are from the true values. Log loss (binary cross-entropy) is commonly used in logistic regression.
    - Confident predictions (when probabilities are very close to 0 or 1) are heavily penalized by if they are wrong (large increase in loss) and heavily rewarded if they are correct (large decrease in loss).
  - **Gradient descent**: iteratively calculates the gradient of the loss function with respect to each coefficient (how the loss would change if the coefficient was slightly increased or decreased) and each coefficient is then adjusted in the direction that reduces loss. This process is repeated until loss stabilizes (i.e., stops decreasing significantly).
  - **Coordinate descent**: minimizes the loss function with respect to one parameter at a time. Parameters are updated sequentially, cycling through each coefficient and adjusting it in the direction that reduces loss until convergence is reached. Since all but one of the parameters are fixed during each step, each coefficient can be updated using a closed-form solution (computed directly via a formula), leading to faster convergence in certain problems (see below).
  - **L1-regularized optimization**: L1 regularization adds a penalty term to the loss function equal to the sum of the absolute values of the coefficients, which is used to encourage sparsity by penalizing larger coefficients. The absolute value function is not differentiable at zero (there is a kink in it).
    - Gradient descent relies on computing the gradient (derivative) of the loss function with respect to the parameters. Since the L1 penalty is not differentiable at zero, gradient descent doesn't work well with L1 regularization.
    - Coordinate descent solves a one-variable L1-regularized optimization problem for each coefficient individually. It uses the soft-thresholding function, which shrinks coefficients toward zero and sets them to exactly zero if they’re not large enough to overcome the L1 penalty. This balances the fitting of the data (loss function) with model simplicity (L1 penalty).
  - **Sparse data optimization**: gradient descent computes the full gradient vector over all features at each iteration, which doesn’t take advantage of sparsity (many features or coefficients being zero). In contrast, coordinate descent updates one feature at a time and can efficiently skip computations for samples where that feature is zero, leading to faster computations especially when the data or the coefficients are sparse.

After convergence, the final coefficients are used to make predictions on unseen data.

### Logistic Regression Hyperparameters

### Classification Performance Metrics
Performance metrics are used to evaluate how well the model performs on the test data set (the prediction quality of the model).

**Confusion matrix**: compares predicted labels versus actual labels of the data.
  - True positive (TP): model correctly predicted positive class
  - False positive (FP): model incorrectly predicted positive class
  - False negative (FN): model incorrectly predicted negative class
  - True negative (TN): model correctly predicted negative class

**Accuracy**: the proportion of total correct predictions out of all predictions. 
  - (TP + TN)/(TP + TN + FP + FN)
  - Can be misleading if the dataset is imbalanced (e.g., if 90% of the samples are healthy and the model predicts healthy for all samples, the model will have 90% accuracy). 

**Balanced accuracy**: the average of sensitivity and specificity. 
  - (sensitivity + specificity)/2
  - Accounts for performance on both classes equally (better than accuracy for imbalanced datasets).

**Sensitivity (recall)**: how well a model identifies positive cases (what proportion of positive cases were correctly identified as positive). 
  - TP/(TP+FN)
  - High sensitivity = few false negatives (the model rarely misses actual positives)

**Specificity**: how well a model identifies negative cases (what proportion of negative cases were correctly identified as negative).
  - TN/(TN+FP)
  - High specificity = few false positives (the model rarely misses actual negatives)

**Precision**: what proportion of the predicted positives that are actually positive.
  - TP/(TP +FP)
  - High precision = few false positives (if the model identifies a case as positive, it probably is)

**F1 score**: the harmonic mean of precision and recall (balances how many predicted positives are correct (precision) with how many actual positives were found (sensitivity). 
  - 2x(precision x recall)/(precision + recall)
  - High F1 score = the model has both high precision and high recall (most of the positive cases have probably been identified and the cases that have identified as positive are probably actually positive)

**ROC (receiver operating characteristic) curve**: plots the true positive rate (sensitivity) vs the false positive rate (1-specificity) at all possible classification thresholds. 
  - The area under the curve (AUC) is a measure of the overall ability of the model to discriminate positives versus negatives across all thresholds, and ranges from 0.5 (random guessing) to 1.0 (perfect classification).

### Key Issues with Microbiome Data When Using Random Forests

**Compositional data**: microbiome data are relative abundances (the values are not independent). Random forests do not rely on features being independent, but correlated features can affect performance and interpretability. Correlated features can decrease efficiency of the model, since they carry very similar information, resulting in less true randomness in the ensemble. Correlated features can also decrease the generalizability of the model by amplifying noise or over-representing relationships. Random forests also tend to assign higher importance to features that are correlated, which can decrease interpretability. 
  - **Solution**: break the constant-sum constraint and treat the data in Euclidean space by using log-ratio transformations like CLR, ALR or ILR. Most log-ratio methods require positive values, so zeros need to be dealt with using pseudocounts or imputation methods.

**Sparsity (zero inflation)**: in microbiome data, many taxa are absent in most samples. Sparse features provide very little splitting power (not informative). Sparse features may also act as noise and cause overfitting (especially if the data is high dimensional and there are a lot of sparse features).
  - **Solution**: filter out rare features (those present in < 10% of the samples), perform feature selection based on feature importance, and/or use transformations that explicitly model zero-inflation (model both structural and sampling zeros, like ZINB or metagenomeSeq(ZIG)).

**High dimensionality**: microbiome data often have many more features than samples. High dimensionality can increase the risk of overfitting, degrade feature importance (important features masked by irrelevant features) and can also increase the computational load. 
  - **Solution**: filter out rare features (those present in < 10% of the samples) and/or perform feature selection based on feature importance.

**Imbalanced classes**: microbiome studies often have fewer cases than controls, which can skew performance metrics and bias predictions toward the majority class. 
  - **Solution**: use class weights to penalize misclassification of the minority class, use stratified cross-validation (to ensure that samples from the minority class are included in each fold) and/or evaluate performance using metrics that emphasize minority class performance (balanced accuracy and AUROC).

**Small datasets**: many microbiome studies have relatively small sample sizes, which results in an increased risk of overfitting (since there is limited data to learn from) and reduced statistical power.
  - **Solution**: perform feature selection to reduce the number of features and/or use repeated k-fold cross-validation to get robust performance estimates.

**Non-linear relationships and interactions**: microbiome feature interaction may be non-linear.
  - **Solution**: random forests naturally handle non-linear relationships and interactions well, but using a sufficient number of trees (≥ 500) can improve detection of complex interactions by allowing all relevant combinations to be explored.

**Batch effects**: if non biological variation (for example, variation due to sample processing or sequencing) is confounded with class labels, the model may learn to predict batch rather than the true biological signal, which can reduce generalizability. 
  - **Solution**: use batch correction methods, such as ComBat or Harmony (but this may impact biological signal if it is confounded with batch), add batch as a categorical variable so that the impact of batch can be modeled and discounted, and/or use grouped cross-validation (all samples from the same batch fall into the same fold, mimicking real-world deployment).

