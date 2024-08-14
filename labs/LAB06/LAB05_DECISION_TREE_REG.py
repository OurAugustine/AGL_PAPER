#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def RMSE(ypred, ytrue):
    return np.sqrt(np.sum((ypred-ytrue)**2)/ypred.shape[0])


import numpy as np

class Node:
    def __init__(self, feature_index=None, threshold=None, left=None, right=None, value=None):
        self.feature_index = feature_index  # Index of feature to split on
        self.threshold = threshold  # Threshold value for the split
        self.left = left  # Left child (subtree)
        self.right = right  # Right child (subtree)
        self.value = value  # Value of the node for leaf nodes (average of target values)


# In[ ]:


class my_DecisionTreeReg:
    def __init__(self, max_depth=None):
        self.max_depth = max_depth

    def fit(self, X, y):
        self.root = self._build_tree(X, y, depth=0)

    def _build_tree(self, X, y, depth):
        num_samples, num_features = X.shape
        if depth == self.max_depth or num_samples <= 1:
            return Node(value=np.mean(y))

        # Find the best split
        best_feature_index, best_threshold = self._find_best_split(X, y)

        # Split the data
        left_indices = np.where(X[:, best_feature_index] <= best_threshold)[0]
        right_indices = np.where(X[:, best_feature_index] > best_threshold)[0]

        # Recursively build left and right subtrees
        left = self._build_tree(X[left_indices], y[left_indices], depth + 1)
        right = self._build_tree(X[right_indices], y[right_indices], depth + 1)

        return Node(feature_index=best_feature_index, threshold=best_threshold, left=left, right=right)

    def _find_best_split(self, X, y):
        num_samples, num_features = X.shape
        best_mse = float('inf') #used to keep track of the best mean MSE found during the split search.
        best_feature_index = None #store the index of the best feature found for splitting the data
        best_threshold = None #store the threshold value of the best split found

        for feature_index in range(num_features):#loops over every feature in the data set 
            thresholds = np.unique(X[:, feature_index])#computes unique threshold values for the current feature

            for threshold in thresholds:#loop iterates over each unique threshold value for the current feature
                # computes the indices of samples where the feature value is <= threshold
                left_indices = np.where(X[:, feature_index] <= threshold)[0]
                right_indices = np.where(X[:, feature_index] > threshold)[0] # greater than the threshold

                if len(left_indices) == 0 or len(right_indices) == 0:
                    continue #skip the rest of the loop iteration n move to nxt threshold
                
                
                left_y = y[left_indices]
                right_y = y[right_indices]

                mse = self._mse(left_y) + self._mse(right_y)
                if mse < best_mse:
                    best_mse = mse
                    best_feature_index = feature_index
                    best_threshold = threshold

        return best_feature_index, best_threshold

    def _mse(self, y):
        return np.mean((y - np.mean(y))**2)

    def predict(self, X):
        return np.array([self._predict_tree(x, self.root) for x in X])

    def _predict_tree(self, x, node):
        if node.value is not None:
            return node.value

        if x[node.feature_index] <= node.threshold:
            return self._predict_tree(x, node.left)
        else:
            return self._predict_tree(x, node.right)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




