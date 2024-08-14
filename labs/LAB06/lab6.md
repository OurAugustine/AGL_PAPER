# 6. Gradient Boosting Decision Trees (GBDT)
This lab is about the Gradient Boosting Decision Trees. We will demonstrate the GBDT algorithm using the `numpy` library with the Hubor loss function.
Your task of this lab is to have a folder named `lab6` in your repository as `labs/lab6`. Inside the folder, you should have the python file named `decision_tree_regression.py` and `gradient_boosting_decision_tree.py`. The `decision_tree_regression.py` file should contain the decision tree regression class and is from your previous lab 5, and the `gradient_boosting_decision_tree.py` file should contain the gradient boosting decision tree class.
* create the `decision_tree_regression.py` file and copy the decision tree regression class from your previous lab 5.
* use the following code and complete the function `mae_gradient_descent` and `mae_loss_gradient` in the `gradient_boosting_decision_tree.py` file.
* the goal is under the folder `labs/lab6`, you should be able to run the code by `python gradient_boosting_decision_tree.py` and get the output of the RMSE value.

```python
import numpy as np
from decision_tree_regression import decisiontree


class GradientBoostingDecisionTree:
    def __init__(self, max_depth=8, n_estimators=10, lr=0.01, max_features=5, delta=1.0):
        self.max_depth = max_depth
        self.n_estimators = n_estimators
        self.lr = lr # Learning rate (shrinkage)
        self.delta = delta # Huber loss delta parameter

        self.gradient_coeff = []
        self.stop = n_estimators

        self.tree0 = decisiontree(max_depth=self.max_depth, max_features=max_features)
        self.trees = [decisiontree(max_depth=self.max_depth, max_features=max_features) for _ in range(self.n_estimators)]
        
    def fit(self, X, y):
        self.M, self.N = X.shape

        self.tree0.fit(X, y)
        y_pred = self.tree0.predict(X).reshape((self.M, 1))
        residue = y - y_pred

        for idx, itree in enumerate(self.trees):
            if np.linalg.norm(residue) < 1e-4:
                self.stop = idx
                break
                
            itree.fit(X, residue)
            ipred = itree.predict(X).reshape((self.M, 1))

            alpha = self.huber_gradient_descent(ipred, residue, self.lr, 150)
            self.gradient_coeff.append(alpha)

            y_pred = np.add(y_pred, self.lr * alpha * ipred)
            residue = y - y_pred
            # Update residue based on Huber loss gradient
            residue = self.huber_loss_gradient(residue, self.delta)

    def predict(self, X_test):
        y_pred = self.tree0.predict(X_test).reshape((-1, 1))
        for idx, itree in enumerate(self.trees):
            if idx == self.stop:
                break
            coeff = self.lr * self.gradient_coeff[idx]
            y_pred = np.add(y_pred, coeff * itree.predict(X_test).reshape((-1, 1)))
        return y_pred

    def huber_gradient_descent(self, a, b, lr, epochs):
        alpha = np.random.randn(1)[0]
        for epoch in range(epochs):
            # update alpha using the gradient of the Huber loss
            grad = np.sum(self.huber_loss_gradient(b - a * alpha, self.delta) * (-a))
            alpha -= lr * grad
        return alpha

    def huber_loss_gradient(self, a, delta):
        # Gradient of the Huber loss
        return np.where(np.abs(a) <= delta, -a, -delta * np.sign(a))
```

```
if __name__ == '__main__':
    import pandas as pd
    X_train = pd.read_csv('../datasets/airfoil/airfoil_self_noise_X_train.csv').values
    y_train = pd.read_csv('../datasets/airfoil/airfoil_self_noise_y_train.csv').values
    X_test  = pd.read_csv('../datasets/airfoil/airfoil_self_noise_X_test.csv').values
    y_test  = pd.read_csv('../datasets/airfoil/airfoil_self_noise_y_test.csv').values

    GBDT = GradientBoostingDecisionTree(n_estimators=200, max_depth=15)
    GBDT.fit(X_train, y_train)

    y_pred = GBDT.predict(X_test)
    print(RMSE(y_pred, y_test))
```
