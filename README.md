# Machine Learning
## :evergreen_tree::deciduous_tree::evergreen_tree: Random Forest
Random forests are an ensemble learning method that aggregates the predictions of multiple decision trees built using bootstrap sampling and random feature selection in order to improve classification (categorical outcomes) or regression (continuous outcomes) performance.
### Overview of Random Forest Construction
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
  - **Random feature selection**: each decision tree is trained on a bootstrap sample, but at each node split, only a random subset of the features (mtry) is selected, thereby decreasing correlation between the trees and improving performance.

The samples in the testing data are then dropped down the tree until they reach a terminal node and then are assigned the label of the majority class in that node (classification) or the mean of the target variable (regression). The predictions made by each tree are then aggregated by majority vote (classification) or average prediction (regression) to make the overall prediction.
  - **Majority vote**: each model votes for a class label and the class with the most votes wins.
  - **Average prediction**: each model outputs a numeric value and the average number is the final prediction.
### Hyperparameters
Hyperparameters are tunable settings that are set before the training process begins. They control how the model is trained, and therefore control the performance of the model.
  - **Feature selection**: the process of identifying and using only the most informative features in the dataset for training a model. 






