from decision_tree_regression import *
from scipy import stats

class randomforest:
    def __init__(self, n_estimators=9, max_depth=8, max_features=8):
        self.n_estimators = n_estimators
        self.forest = []
        for _ in range(self.n_estimators):
            self.forest.append(decisiontree(max_depth=max_depth, max_features=max_features))

        self.X = None
        self.y = None
        self.N = 0
        self.M = 0

    def fit(self, X, y):
        self.M, self.N = X.shape

        self.trees_idx = np.random.randint(0, self.M, size=(self.n_estimators, self.M))
        for i, itree in enumerate(self.forest):
            itree.fit(X[self.trees_idx[i]], y[self.trees_idx[i]])

    def predict(self, X_test):
        n_test = X_test.shape[0]
        y_pred = np.zeros((self.n_estimators, n_test))
        for i, itree in enumerate(self.forest):
            y_pred[i, :] = itree.predict(X_test)
        return np.mean(y_pred, axis=0)

if __name__ == '__main__':
    import pandas as pd
    X_train = pd.read_csv('airfoil_self_noise_X_train.csv').values
    y_train = pd.read_csv('airfoil_self_noise_y_train.csv').values
    X_test  = pd.read_csv('airfoil_self_noise_X_test.csv').values
    y_test  = pd.read_csv('airfoil_self_noise_y_test.csv').values

    RF = randomforest()
    RF.fit(X_train, y_train)
    y_pred = RF.predict(X_test)
    print(RMSE(y_pred.ravel(), y_test.ravel()))
