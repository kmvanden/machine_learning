# Machine Learning
## :cowboy_hat_face::mountain::spider_web: Logistic Regression
Logistic regression is a statistical classification algorithm that models the relationship between predictor variables and outcome labels.

### Overview of Logistic Regression Model Construction
The input features are combined linearly and coefficients (weights) are applied to each predictor variable, creating the linear predictor.
  - Linear predictor = β0​+β1​x1​+β2​x2​+⋯+βn​xn​ (for a binary logistic regression model)
  - Coefficients are typically initialized at zero (or at very small random values), resulting in initial predicted probabilities near 0.5 (unbiased starting point).

The logistic regression model assumes that the log-odds of the probability of the outcome being in the positive class is a linear combination of the input features.
  - Linear predictor = β0​+β1​x1​+β2​x2​+⋯+βn​xn​ = log(p/1-p​) = log odds
  - Log-odds are used as an intermediary rather than directly modeling probabilities because:
    - Linearity in parameters: in log-odds (logit) space, the model can treat the relationship between the features and the outcomes as a linear combination of features (i.e., the linear predictor) and linear models are easier to analyze, interpret and optimize.
    - Unbounded domain: log-odds are unbounded, whereas probabilities have to be between 0 and 1. Unbounded values make gradient optimization more effective and stable (gradients behave smoothly across the entire real number line, avoiding saturation near probability bounds where gradients become increasingly small) and facilitate closed-form updates in coordinate descent (there are no box constraints, so constrained optimization is avoided).
    - Convex loss functions: log-odds enable the use of convex loss functions (like log loss), which have a single global minimum, making optimization more reliable (they are guaranteed to converge at the global minimum).

To convert the log-odds into probabilities, the log-odds transformation is inverted using the logistic (sigmoid) function, which maps any real number to a value between 0 and 1, producing an S-shaped curve.
  - Probability = 1/(1+ e<sup>-(β0​+β1​x1​+⋯+βn​xn​)</sup>)
  - Small changes in the log-odds around zero result in rapid changes in probability near 0.5, allowing smooth transitions between class probabilities, and the curve flattens at extreme values, reflecting confidence in the classification.

Optimization algorithms, such as gradient descent or coordinate descent, are used to learn the coefficients that minimize the loss function, aiming to produce predicted probabilities that match the true class labels as closely as possible.
  - **Loss function**: a mathematical expression that measures how far off the predictions of the model are from the true values. Log loss (binary cross-entropy) is commonly used in logistic regression.
    - Confident predictions (when probabilities are very close to 0 or 1) are heavily penalized if they are wrong (large increase in loss) and heavily rewarded if they are correct (large decrease in loss).
  - **Gradient descent**: iteratively calculates the gradient of the loss function with respect to each coefficient (how the loss would change if the coefficient was slightly increased or decreased) and each coefficient is then adjusted in the direction that reduces loss. This process is repeated until loss stabilizes (i.e., stops decreasing significantly).
  - **Coordinate descent**: minimizes the loss function with respect to one parameter at a time. Parameters are updated sequentially, cycling through each coefficient and adjusting it in the direction that reduces loss until convergence is reached. Since all but one of the parameters are fixed during each step, each coefficient can be updated using a closed-form solution (computed directly via a formula), leading to faster convergence in certain problems (see below).
  - **L1-regularized optimization**: L1 regularization adds a penalty term to the loss function equal to the sum of the absolute values of the coefficients, which is used to encourage sparsity by penalizing larger coefficients. The absolute value function is not differentiable at zero (there is a kink in it).
    - Gradient descent relies on computing the gradient (derivative) of the loss function with respect to the parameters. Since the L1 penalty is not differentiable at zero, gradient descent doesn't work well with L1 regularization.
    - Coordinate descent solves a one-variable L1-regularized optimization problem for each coefficient individually. It uses the soft-thresholding function (closed-form update step), which shrinks coefficients toward zero and sets them to exactly zero if they’re not large enough to overcome the L1 penalty. This balances the fitting of the data (loss function) with model simplicity (L1 penalty).
  - **Sparse data optimization**: gradient descent computes the full gradient vector over all features at each iteration, which doesn’t take advantage of sparsity (many features or coefficients being zero). In contrast, coordinate descent updates one feature at a time and can efficiently skip computations for samples where that feature is zero, leading to faster computations especially when the data or the coefficients are sparse.

After convergence, the final coefficients are used to make predictions on unseen data.

In logistic regression, each coefficient represents the log-odds change in the outcome associated with a one-unit increase in the predictor, assuming other features are held constant. A positive coefficient means that an increase in the feature is associated with an increased probability of the positive class and a negative coefficient means that an increase in the feature is associated with a decreased probability of the positive class. 
  - For example, if the log-transformed abundance of a bacterial taxon has a coefficient of 1.5, that implies that for every unit increase in log abundance, the odds of the positive class increase by a factor of approximately 4.5 (e<sup>1.5</sup> ≈ 4.5). 
  - However, the interpretability of logistic regression coefficients can be complicated by several factors.
    - If the data is high dimensional and sparse, the coefficients become unstable and unreliable without appropriate regularization.
    - If the data is highly correlated, they share explanatory power, causing the individual coefficients to lose unique interpretability.
    - And while regularization improves stability, it complicates direct interpretation of the coefficients, as the magnitudes of the coefficients are affected by the penalty term in addition to the underlying effect size. 

### Logistic Regression Hyperparameters
In logistic regression, key hyperparameters relate to regularization (addition of a penalty term to the loss function), which helps to control model complexity and to prevent overfitting by discouraging large or unnecessary coefficients.

**L1 regularization (lasso)**: adds a penalty equal to the sum of the absolute values of the coefficients to the loss function.
  - **Feature selection**: L1 regularization introduces a kink (a non-differentiable point) at zero in the absolute value function (penalty function), which encourages sparse solution (pushes some coefficients to exactly zero), effectively performing feature selection. 
  - **Sparse data**: L1 naturally ignores irrelevant features by assigning zero coefficients.
  - **Prevents overfitting**: the larger the absolute value of the coefficient, the higher the penalty. Thus, only features that significantly improve model accuracy are assigned non-zero weights, which helps to prevent overfitting and improve generalization.
  - **Limitations**: to minimize the L1 penalty, the optimizer is encouraged to set as many coefficients as possible to zero. Thus, if two or more features are correlated, L1 often randomly selects only one of the features and sets the coefficients of the other features to zero, making feature selection for correlated features unstable.

**L2 regularization (ridge)**: adds a penalty equal to the sum of squares of the coefficients to the loss function. 
  - **Prevents overfitting**: the L2 penalty grows quadratically with the size of the coefficient, thus encouraging the model to keep the coefficients small, but without forcing them to be exactly zero, thus retaining all features within the model (no feature selection).
  - **Correlated features**: L2 shares the weights among all correlated features rather than picking one and discarding the others. This tends to make the model more stable and less sensitive to small changes in the data compared to L1.
  - **Improves numerical stability**: the covariance matrices used for coefficient estimation for high dimensional and/or highly correlated data are nearly singular (ill-conditioned) and inversion of a nearly singular matrix is numerically unstable (estimated coefficients are sensitive to noise and small data changes). High dimensional data (p>n) is linearly dependent because you can only have n independent features in n dimensional space, thus p-n features have to be linear combinations of the other features. L2 regularization adds an identity matrix (scaled by regularization strength) to the covariance matrix, which increases the diagonal elements. This increases the eigenvalues and improves the condition of the matrix, thereby making the inversion more stable.
  - **Limitations**: does not perform feature selection, which can be a major limitation when modeling high dimensional data where interpretability is important.

**L1 and L2 regularization (elastic net)**: combines the penalties of L1 and L2 regularization. The combination of L1 and L2 penalties make elastic net more stable with high dimensional and highly correlated data.
  - **Feature selection**: L1 regularization pushes some coefficients to exactly zero.
  - **Improves numerical stability**: L2 regularization improves coefficient stability with data that is high-dimensional and highly correlated.

**Alpha**: determines the balance between L1 and L2 regularization used in elastic net.   
  - When alpha = 1, the model is pure L1. When alpha = 0, the model is pure L2. When 0 < alpha < 1 the model is a mixture of L1 and L2 (elastic net). 
  - Different data sets benefit from different regularization types. Tuning alpha allows the model to find the optimum balance between L1 and L2 regularization for the given dataset.

**Lambda**: controls the strength of the regularization penalty. 
  - Very low values result in a model that resembles the unregularized model and is prone to overfitting, whereas high values can shrink the coefficients to zero (especially with L1 regularization), resulting in an underfitted model.
  - Lambda is typically tuned over a logarithmic scale because small lambda values can lead to large shifts in model coefficients, whereas large values often result in flat model behavior (all coefficients are near or at zero).

**Class weights**: when classes are imbalanced, logistic regression models tend to minimize loss by favoring the majority class, often at the expense of the minority class. Class weights correct for this by assigning higher weights to samples from the minority class in the loss function, thus the penalty for misclassifying minority samples is amplified.

## :evergreen_tree::deciduous_tree::evergreen_tree: Random Forest
Random forests are an ensemble learning method that aggregates the predictions of multiple decision trees built using bootstrap sampling and random feature selection in order to improve classification (categorical outcomes) or regression (continuous outcomes) performance.

### Overview of Random Forest Model Construction
Data is split into training data and testing data. 
  - **Training data**: the portion of the data used to train (build) the model. The model learns the patterns/relationships between the input features and the target labels.
  - **Testing data**: the portion of the data that is used to evaluate the model (determine how well the model generalizes to unseen data).
  - **Cross-validation**: resampling technique used to evaluate the performance of a model on unseen data and avoid overfitting (by estimating how well the model might perform on unseen data).
    - **Overfitting**: when a model performs well on training data, but poorly on testing data (the model has learned noise or specific patterns only present in the training set that are not generalizable to other datasets).
    - **k-fold cross-validation**: splits the data into k equal parts (folds), trains the model on k-1 folds and tests the model on the remaining fold. The data used in the folds is then rotated, and the process is repeated k times so that each fold is used once as the test dataset. The k-fold process can be repeated multiple times with different fold splits for more reliable estimates.
    - **Stratification**: ensures that the class distribution in each fold of the cross-validation split reflects the overall distribution of the dataset so that each training and testing fold has a representative mix of classes.
    
Multiple subsets of the training data are created using bootstrap sampling and each subset is used to build an individual tree. 
  - **Bootstrap sampling**: generates multiple new datasets from the training data using random sampling with replacement. Each bootstrap sample is the same size as the original dataset, but some of the samples are repeated and some are omitted. 

At each node split in a decision tree, only a random subset of features is considered. 
  - **Random feature selection**: each decision tree is trained on a bootstrap sample, but at each node split, only a random subset of the features (```mtry```) is selected. This decreases the correlation between the trees, which reduces variance and improves generalization.

The samples in the testing data are then dropped down the tree until they reach a terminal node and then are assigned the label of the majority class in that node (classification) or the mean of the target variable (regression). The predictions made by each tree are then aggregated by majority vote (classification) or average prediction (regression) to make the overall prediction.
  - **Majority vote**: each model votes for a class label and the class with the most votes wins.
  - **Average prediction**: each model outputs a numeric value and the average number is the final prediction.

### Random Forest Hyperparameters
Each decision tree within a random forest is trained on a bootstrap sample of the training data. This introduces variation, so that each tree doesn’t learn the same patterns. On average, about two-thirds of the training data is included in any given bootstrap sample (the remaining one-third of the data are the out-of-bag (OOB) samples).

Features within the data (especially in high dimensional datasets) can be irrelevant to class separation. Models trained with many irrelevant features are more likely to learn patterns specific to the training dataset and thus not generalize well to unseen data. 
  - **Feature selection**: the process of identifying and using only the most informative features in the dataset for training a model.
    - **Recursive feature elimination**: a model is iteratively trained, removing the least important features (based on mean decrease in Gini impurity or mean decrease in accuracy) at each iteration and comparing model performance for the different numbers of features used to train the model.
    - **Mean decrease in Gini impurity**: measures how much a feature contributes to decreasing Gini impurity (how pure or homogeneous a set of classes is at a given node in a decision tree) when it is used to split data across all of the trees in the forest.
      - High mean decrease in Gini impurity = splitting the data by the feature results in nodes with increased purity
    - **Mean decrease in accuracy**: measures how much model accuracy drops when the feature is randomly permuted (permutation importance). Values of the feature in the OOB samples are permuted and the accuracy of the model before and after the permutation is measured.
      - High mean decrease in accuracy = the feature contributes a lot to the accuracy of the predictions made by the model
    - Both mean decrease in Gini and mean decrease in accuracy can be biased toward continuous features and categorical features with many levels (more unique values = more opportunities to split the data = more chances to fit the target or reduce impurity, even by chance) and highly correlated features (they share predictive information, so the decrease in accuracy/impurity is masked).

Each tree is grown by recursively splitting the data at internal decision nodes, starting at the root node. With random forests, only a random subset of features is considered at each split to help decorrelate the trees and thereby increase the robustness of the ensemble.
  - **Number of features per split** (```mtry```): number of features randomly selected at each node split. The default is the square root of the number of features for classification tasks and one third of the features for regression tasks. 
    - Lower ```mtry``` = increases tree diversity, which reduces correlation between trees (helps reduce overfitting) 
    - Higher ```mtry``` = potentially better split quality, but higher correlation between the trees
   
At each split, among the selected mtry features, the algorithm chooses the feature and the threshold that best separates the classes (how homogenous the resulting nodes are after the split). For classification, the goal is to reduce impurity (which is typically measured using Gini impurity), and for regression the goal is to minimize variance in the child nodes. However, if the classes are imbalanced, the model will favor majority class predictions and the model will have poor sensitivity for the minority class. 
  - **Imbalanced dataset**: target classes (labels) are not represented equally within the dataset. Standard models tend to optimize for overall accuracy, which favors the majority class, resulting in the minority class being underpredicted.
  - **Class weights** (```classwt```): mitigate the effects of the imbalance between class labels by adjusting the penalty for misclassifying samples of different classes. Essentially, class weights make the model pay more attention to the minority class by increasing the impact (effective contribution) of the minority class on impurity calculations and thereby influencing which splits are considered good.
    - If the class weights are too extreme, the model may overfit the minority class.

The splitting process continues recursively for each new node until a stopping rule is met. Typical stopping conditions include: all samples in a node belonging to the same class or the node containing fewer samples than ```nodesize```.
  - **Minimum node size for a split** (```nodesize```): minimum number of samples required to split an internal node. ```nodesize``` controls tree depth and complexity.
One is the default for classification tasks and five is the default for regression tasks. 
    - Smaller values = deeper trees with more splits that capture more detail, but risk overfitting the data (especially in noisy or small datasets).
    - Larger values = shallower trees with fewer splits that capture less detail and may underfit the model (miss patterns in the data).

The decision tree-building process is then repeated ```ntree``` times to create an ensemble of decision trees, each built on a different bootstrap sample of the training data, and each making independent predictions that are then aggregated.
  - **Number of trees** (```ntree```): determines how many decision trees are in the forest. 
More trees = more stable and consistent predictions, improvement of generalization and accuracy (up to a point), but requires longer computation time

## :evergreen_tree::deciduous_tree::rocket: XGBoost
XGBoost (extreme gradient boosting) is an ensemble learning method that uses a second-order Taylor expansion of the loss function to sequentially optimize shallow decision trees. 

### Overview of XGBoost Model Construction
Similarly to random forests, the data is split into training data and testing data. However, unlike random forests, XGBoost trees see the full training dataset (no sampling), unless ```subsample``` is set to a value less than one, and all features are considered (no random feature selection), unless one of the following subsampling parameters is set to a value less than one: ```colsample_bytree```, ```colsample_bylevel```, or ```colsample_bynode```.
  - **Training data**: the portion of the data used to train (build) the model. The model learns the patterns/relationships between the input features and the target labels.
Testing data: the portion of the data that is used to evaluate the model (determine how well the model generalizes to unseen data).
  - **Cross-validation**: resampling technique used to evaluate the performance of a model on unseen data and avoid overfitting (by estimating how well the model might perform on unseen data).
    - **Overfitting**: when a model performs well on training data, but poorly on testing data (the model has learned noise or specific patterns only present in the training set that are not generalizable to other datasets).
    - **k-fold cross-validation**: splits the data into k equal parts (folds), trains the model on k-1 folds and tests the model on the remaining fold. The data used in the folds is then rotated, and the process is repeated k times so that each fold is used once as the test dataset. The k-fold process can be repeated multiple times with different fold splits for more reliable estimates.
    - **Stratification**: ensures that the class distribution in each fold of the cross-validation split reflects the overall distribution of the dataset so that each training and testing fold has a representative mix of classes.

XGBoost sequentially improves its predictions by adding new trees that minimize the objective function. Within the objective function, the loss function is approximated as a quadratic function of the leaf weight using a second-order Taylor expansion. This allows for a closed-form solution for the optimal logit correction, avoiding the need for iterative optimization within every tree node. 
  - **Objective function** = loss function + regularization
  - **Loss function**: measures how close the predicted probabilities match the true class labels (0 or 1) and how confident the model is in those predictions. Log loss is typically used by XGBoost for binary classification. 
    - The loss is low when the predictions are both accurate and confident, and high when the predictions are confident but wrong. 
    - Loss ≈ loss(current logit) + logit correction x gradient + ½ x (logit correction)<sup>2</sup> x Hessian
      - **Gradient** = the first derivative of the loss function with respect to the current logit. It represents the slope of the loss function and determines the direction and magnitude of the adjustment needed to reduce the loss.
      - **Hessian** = the second derivative of the loss function with respect to the current logit. It represents the curvature (confidence) of the loss function and qualifies how aggressively the correction should be applied (step size)).
  - **Regularization**: penalties that control complexity and prevent overfitting.

The model starts with an initial logit based on class distribution (e.g., if 60% of the training labels are in class 1, then the logit = log(0.6/0.4) = 0.405), which is converted into a probability using the sigmoid function (e.g., probability = 1/(1 + e<sup>-0.405</sup>) = 0.60). The probability is then used to compute the gradients and the Hessians for each sample. The gradient and the Hessian for each sample are fixed before each subsequent tree is built. 

  For binary classification using log loss:
  - **Gradient (current sample)** = predicted probability (current sample) - true label (current sample)
    - For example: if the gradient = 0.6 - 1 = - 0.4, the prediction is too low and needs to be adjusted upward (magnitude and direction of the correction).
  - **Hessian (current sample)** = predicted probability (current sample) * (1 - predicted probability (current sample)) 
    - For example, if the Hessian = 0.95(1 - 0.95) = 0.0475, the model is confident, so smaller updates should be made to avoid overshooting (confidence of the correction).

For each parent node, the split (feature and threshold) that maximizes the gain function (i.e., minimizes the objective function) is chosen. The algorithm uses a greedy strategy (the locally optimal split is chosen without regard for future splits). After the tree is fully grown, weak splits (i.e., where the gain is less than ```gamma```) are pruned from the leaf to the root. L1 regularization is not included in the gain function, because it involves an absolute value term, which is not differentiable at zero, and therefore requires soft-thresholding to allow a closed-form solution (which would be computationally expensive).
  - **Gain** = objective function (parent node) - (objective function (left child node) + objective function (right child node)) - ```gamma```

The decision tree effectively groups together samples with similar gradients and Hessians, as these samples require similar logit corrections (direction and magnitude) to their current predictions in order to minimize the loss function. This process continues recursively until the tree reaches a stopping condition (e.g., ```max_depth``` or ```min_child_weight```). 

Within each leaf (terminal node), the optimal leaf value (logit correction), the one that minimizes the objective function, is found using a closed-form solution. 
  - If L1 regularization is applied, soft-thresholding is used to handle the non-differentiability at zero, enabling a closed-form solution.

The calculated logit correction is scaled by the learning rate (```eta```) and is prevented from exceeding the value set by ```max_delta_step```. The logit correction is then added to the current logit prediction of every sample that ended up in that leaf. The new logit predictions are then converted into probabilities using the sigmoid function and used to build the subsequent tree. This allows each tree to learn from the mistakes of the previous tree (boosting). 

Trees are added sequentially until a stopping condition is reached, such as hitting the predefined number of rounds (```nrounds```) or satisfying early stopping criteria (```early_stopping_rounds```). The prediction for a sample is the sum of the initial logit and all of the logit corrections from the sequence of trees. For binary classification, the sigmoid function is applied to this final logit to produce a probability.

### Feature Importance for XGBoost
Model-agnostic methods, like ```fastshap```, need to consider all possible subsets of features to compute Shapley values, which is exponential in complexity (2<sup>number of features</sup>) and therefore typically rely on approximation methods like Monte Carlo sampling.
  - **Shapley values**: measure the average marginal contribution of each feature across all possible combinations of features (subsets) to the model prediction. They explain how each feature contributes the output of the model (prediction) for a specific sample (local interpretability) and can be aggregated across samples to explain how each feature influences the model overall (global feature importance).

Tree SHAP (SHapley Additive exPlanations) is an exact algorithm for computing SHAP values specifically for tree-based models. It uses the structure of decision trees (i.e., tree splits, leaf values, tree structure/sample paths and node cover) to compute Shapley values in polynomial time.
  - **Tree splits**: which features and thresholds are used at each node to route samples.
  - **Leaf values**: output value at terminal nodes (e.g., logit correction for binary classification).
  - **Tree structure/sample paths**: how each sample flows from root to leaf based on splits.
  - **Node cover**: the number of training samples that pass through each node.

For each tree, Tree SHAP determines the change in prediction (leaf value) for a given sample caused by each feature used in a split along the path taken by that sample. The change in prediction is the difference between the prediction obtained using the actual path taken by the sample and the expected prediction without the feature (weighted by child node cover). 

The difference in the prediction value is weighted by the probability that this split would be encountered by a subset of features that includes the feature used in the split (i.e., that the probability that a random subset of features contains the feature of interest and that this subset of features reaches that node). Tree SHAP derives this probability using node cover and the structure of the path from the root of that tree to the current node, specifically how many nodes occur before the current node and how many unique features are used in split decisions in these nodes. The SHAP values are summed across all trees in order to get the final SHAP values for per feature for each sample.

### XGBoost Hyperparameters
In XGBoost, important hyperparameters relate to weight regularization (L1 and L2), structural regularization and loss-guided post-tree pruning are used to prevent overfitting by discouraging large leaf weights and overly complex tree structure.

**Subsampling** (```subsample```): controls the fraction of the training data that is randomly selected (without replacement) to build each tree. This introduces randomness into the training process, which helps to prevent overfitting. The default is to use all of the data (1.0). Values between 0.5 - 0.8 can potentially reduce the risk of overfitting.

**Fraction of features** (```colsample_bytree```): controls the fraction of features (columns) randomly selected to train each tree. This decreases the correlation across the trees, which helps to reduce the risk of overfitting. The default is to use all of the data (1.0). Values between 0.5 - 0.8 can potentially reduce the risk of overfitting.
  - If stronger regularization is needed, the following hyperparameters can also be modified: 
    - ```colsample_bylevel```: controls the fraction of features randomly selected per level (often used when building deeper trees).
    - ```colsample_bynode```: controls the fraction of features randomly selected per node (typically used if the dataset is high dimensional).

**L1 regularization** (```alpha```): adds a penalty equal to the sum of the absolute values of the leaf weights scaled by ```alpha``` to the objective function. L1 encourages sparsity and prevents overfitting by pushing small leaf weights to exactly zero. This reduces overfitting by pruning splits that contribute little to no model performance. Start by setting ```alpha``` to zero (no L1 regularization), then increase the values gradually to observe whether the performance metrics improve.

**L2 regularization** (```lambda```): adds a penalty equal to the sum of the squared leaf weights scaled by ```lambda``` to both the objective function and the gain function. L2 shrinks all leaf weights toward zero, discouraging large logit corrections, which helps to prevent overfitting. L2 also shrinks the magnitude of the gain for all splits, which can result in irrelevant or noisy splits being pruned from the tree if the gain is below the threshold set by ```gamma```. Start by setting ```lambda``` to zero (no L2 regularization), then increase the values gradually to observe whether the performance metrics improve.

**Minimum split gain** (```gamma```): minimum gain needed for a split to be kept. After the entire tree is built in a greedy fashion, XGBoost performs loss-guided post-pruning starting from the last non-leaf node. If a split with a gain less than ```gamma``` leads to a split with a gain greater than ```gamma```, the subtree will be kept. ```gamma``` removes weak splits that don’t significantly improve the model, but doesn’t pre-prune trees before seeing their full potential. Start by setting it to zero (no pruning) and increase the values to reduce tree complexity.

**Learning rate** (```eta```): scales the contribution of the logit correction of each new tree to the overall model prediction. Lower values of ```eta``` (e.g., decreasing the default value of 0.3 to 0.1, 0.05 or 0.01) make each step more conservative, but reduces the risk of overfitting. Lower values of ```eta``` generally require more trees (```nrounds```) to reach optimal performance and should be used with early stopping (```early_stopping_rounds```) during training to avoid wasting computation on unnecessary trees.

**Maximum rounds** (```nrounds```): sets the maximum number of trees (boosting iterations) to perform (typically set between 500 and 1000). If it is set too low, the model will be underfit (not enough trees to capture the patterns), and if it is set too high, the model can be overfit (learns noise), especially if regularization and early stopping are not used. A smaller learning rate (```eta```) requires more trees in order to converge to a good solution.

**Early stopping** (```early_stopping_rounds```): stops training early if the performance of the model on the testing dataset has not improved (based on a chosen performance metric) after a specified number of consecutive rounds (e.g., 10 - 50, depending on dataset size and training noise). This helps to prevent overfitting and speeds up training by skipping unnecessary tree building. The final model only includes the trees up to the best iteration found during training. Using this parameter also avoids the need to guess how many trees (```nrounds```) would be ideal. Instead, set ```nrounds``` generously high and allow early stopping to determine the actual number of useful trees.

**Maximum depth** (```max_depth```): sets the maximum depth of each individual decision tree built during the boosting process. Lower values result in shallower trees that improve generalization, but have an increased risk of underfitting. Whereas higher values result in deeper trees that capture more patterns, but have an increased risk of overfitting. In boosting, each tree is meant to make small corrections, not learn the full structure like in random forests. Therefore, shallow trees (3 - 10) often work better.
  - **Maximum leaves** (```max_leaves```): sets the maximum number of terminal nodes (leaves) that a decision tree is allowed to have. It is an alternative to ```max_depth``` for controlling the depth of the tree. It is more flexible and often preferred for highly imbalanced or high-dimensional datasets (prevents over-specialized deep branches that fit noisy data and/or the majority class). For small trees, set ```max_leaves``` to 16, 32 or 64.

**Minimum sum of Hessians** (```min_child_weight```): the minimum sum of Hessians needed in a child node to allow a split. It prevents the model from creating nodes that only represent a very small number of observations, even if it reduces loss a little bit (which often leads to overfitting). For classification, the sum of Hessians reflects confidence in the predictions. Thus, a small sum of Hessians indicates that the model doesn’t trust the gain from the split enough to justify the complexity. Start by setting it to one (default) and increase the values if you observe overfitting. ```min_child_weight``` interacts with class imbalance. If most samples are confidently predicted early (majority class), their Hessians become small. Therefore, setting this hyperparameter too high may prevent valid splits. 

**Maximum logit correction** (```max_delta_step```): limits the maximum logit corrections that can be taken during each boosting iteration. If the sum of the Hessians is very small (which can be caused by imbalanced datasets (probability set close to 0 or 1 early in training) or by saturated sigmoid outputs (probability approaching 0 or 1)), the leaf values end up being very large and learning becomes destabilized. The default is zero (no constraint). For imbalanced data, try values between 1 - 10.

**Class imbalance** (```scale_pos_weight```): scales the gradients and the Hessians for the positive class (label = 1) during training. An imbalanced dataset can be biased toward the majority class and underperform on the minority class. To counter this, ```scale_pos_weight``` is commonly set to the ratio of the number of negative samples divided by the number of positive samples. This number is used to scale the gradients and the Hessians of the positive class. Larger gradients and Hessians contribute more to the gain function, resulting in increased gain from splits that separate the minority class from the majority class, and thus encouraging the model to put more importance on the minority class samples.

## :bar_chart::medal_sports: Classification Performance Metrics
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

## :compass::wrench: Bayesian Optimization of Hyperparameters
**Initial scoring function**: ```initPoints``` points are chosen uniformly at random from the predefined hyperparameter bounds (i.e., each hyperparameter is independently sampled from a uniform distribution over its respective range). Because of this independent random sampling, especially in high-dimensional spaces, the initial points may exhibit clustering or sparse coverage (gaps), particularly when initPoints is small relative to the number of dimensions. Each of the sampled hyperparameter combinations (input) are then evaluated with the scoring function to get their respective scores (output). These input-output pairs form the training data for a Gaussian Process (GP), a probabilistic model used to approximate the unknown objective function. The GP assumes that the function values (scores) at all input locations are jointly Gaussian distributed (i.e., they follow a multivariate normal distribution). It uses a kernel function (covariance function) to model how correlated the scores of two points are based on how close they are within the input (hyperparameter) space. Using this information, the GP is then able to fit a posterior distribution over all possible functions that are consistent with the observed data, thus estimating the mean prediction (expected score) and the standard deviation (uncertainty of the model) for any new point in the input space. 

**Local optimum search**: The GP gives a posterior distribution over the objective function and the acquisition function (e.g., expected improvement) uses the GP posterior to determine how promising any given point would be to evaluate next. Expected improvement balances exploration (points with high uncertainty) and exploitation (points with high predicted performance). Since the acquisition function is a non-convex, nonlinear function in high-dimensional space that may have multiple local maxima, it has to be numerically (iteratively) optimized in order to find the point in the hyperparameter space where it is maximized. The mean and variance predictions from the GP are smooth functions of the input and the hyperparameters have defined bounds. L-BFGS-B is a commonly used, gradient-based optimization algorithm for this type of problem. The optimizer searches the local hyperparameter space via multiple starting points to find the locations with the highest acquisition values. The number of candidate points that are searched depends on the dimensionality of the space, and the acquisition function and optimizer method used. 

**Scoring function**: The top ```iters.k``` number of points (those with the highest acquisition values) are selected from the local optimum search and evaluated using the scoring function. The results of the scoring function are then used to update the GP for the next epoch.

## :thinking::memo: Key Issues with Microbiome Data When Using Machine Learning Models
**Compositional data**: microbiome data consists of relative abundances (the features are not independent), which introduces spurious correlations between the features due to the constant-sum constraint. Random forests and logistic regression do not assume that features are independent, but their performance and interpretability can still be impacted by the presence of correlated features. Correlated features (both real biological relationships and those caused by the closure effect) provide redundant information to the model, which can cause the model to overfit small patterns or noise, decreasing generalizability on unseen data. Random forests and logistic regression (with L1 regularization) tend to randomly assign importance to one feature (or sometimes a few features in the case of random forests) out of a group of correlated features, even if they are equally informative, which can lead to unstable or misleading interpretations. Additionally, this redundant information results in less true randomness (created during bootstrap sampling and random features selection) in the ensemble, which decreases the diversity of decision trees within a random forest.
  - **Solution**: break the constant-sum constraint and treat the data in Euclidean space by using log-ratio transformations like CLR, ALR or ILR. Most log-ratio methods require positive values, so zeros need to be dealt with using pseudocounts or imputation methods. For logistic regression, the use of L2 regularization can help with the stability and interpretability of correlated features. For random forests, mean decrease in accuracy (permutation importance) is a better measure of feature importance for correlated features than mean decrease in Gini impurity (correlated features compete to be selected at the splits). 

**Sparsity (zero inflation)**: in microbiome data, many taxa are absent in most samples, resulting in sparse (zero-dominated) features. Sparse features provide very little splitting power (random forests) and very little weight to the linear predictor (logistic regression). Moreover, sparse features may act as noise and increase the risk of overfitting, especially if the data is high dimensional and there are a lot of sparse features).
  - **Solution**: filter out rare features (those present in < 10% of the samples) and/or perform feature selection (apply L1 regularization in logistic regression or eliminate features based on feature importance in random forest).

**High dimensionality**: microbiome data often have many more features than samples. High dimensionality can degrade feature importance (important features masked by irrelevant features/noise), increase the risk of overfitting (more likely to model spurious patterns, which decreases generalizability), and can increase the computational load.
  - **Solution**: filter out rare features (those present in < 10% of the samples) and perform feature selection (apply L1 regularization in logistic regression or eliminate features based on feature importance in random forest). For random forests, mean decrease in accuracy (permutation importance) is a better measure of feature importance for high dimensional data than mean decrease in Gini impurity (biased toward more frequent appearances in splits, which increases with dimensionality). 

**Imbalanced classes**: in microbiome studies, it is common to have fewer cases than controls, leading to class imbalance. When classes are imbalanced, logistic regression models tend to minimize loss by favoring the majority class, since misclassifying minority samples contributes less to the overall loss. Similarly, classification tasks in random forest also tend to favor the majority class, because the splitting criteria (Gini impurity) are dominated by the majority class.
  - **Solution**: use class weights to correct for this by assigning higher weights to samples from the minority class in the loss function in logistic regression or at split decisions and class probabilities in leaf nodes in random forest. This increases the penalty for misclassifying minority samples and forces the model to place more importance on minority class predictions. Stratified cross-validation should be used to ensure that samples from the minority class are included in each fold. Additionally, performance of the model should be evaluated using metrics that emphasize minority class (precision, recall and F1 score) or balanced (balanced accuracy and AUROC) performance.

**Small datasets**: many microbiome studies have relatively small sample sizes. In unregularized logistic regression, the model estimates one coefficient per feature, and if the number of features exceeds the number of samples, the optimization problem has no unique minimum. Since there are infinite solutions, logistic regression will find coefficients that fit the training data exactly, resulting in models that capture noise and spurious patterns (overfitted). In random forests, the decision trees are high variance learners: a greedy algorithm is used to find the best splits. Thus very different tree structures can result from small changes in the data that cause small changes in split impurity. Smaller datasets can increase this instability (fewer samples in the training set to learn from). Additionally, smaller datasets increase the overlap of observations shared between bootstrap samples, which decreases the diversity of the decision trees and thereby increases the risk of overfitting.
  - **Solution**: regularization during logistic regression constrains the coefficient values, breaking the infinite solution problem by selecting a simpler unique solution. Feature selection (L1 regularization in logistic regression or elimination of features based on feature importance in random forest) filters irrelevant features and prevents the model from fitting noise. Additionally, repeated cross-validation should be used with small sample sizes to avoid overfitting and to get robust performance estimates.

**Non-linear relationships and interactions**: microbiome datasets frequently exhibit complex, nonlinear relationships between microbial features and outcomes, due to microbial interactions, threshold effects and compositionality. Logistic regression models the outcome as a linear combination of the input features and therefore cannot naturally capture nonlinear patterns or interactions. In contrast, random forests make no assumptions about linearity and are able to capture nonlinear thresholds and interactions between taxa.
  - **Solution**: random forests naturally handle non-linear relationships and interactions, but stability and the detection of complex combinations can be improved by using a sufficient number of trees to capture subtle interactions. 

**Batch effects**: non biological variation introduced during sample collection, processing or sequencing can lead to batch effects. If batch is confounded with class labels, models may learn to predict batch artifacts rather than true biological differences, leading to poor generalizability and spurious associations.
  - **Solution**: Apply batch correction methods, such as ComBat or Harmony, to adjust the data prior to modeling. However, if biological signal is confounded with batch, these corrections can remove meaningful variation along with technical noise. Include batch as a covariate (logistic regression) or as a stratification variable (random forest) to allow the model to adjust for barch influence without needing to remove it from the input features. Use grouped cross-validation (all samples from the same batch are assigned to the same fold), which mimics real-world deployment, where the model encounters unseen data.
