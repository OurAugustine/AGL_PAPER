from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
data = pd.read_csv('../datasets/Boston.csv')

X = data.drop('medv', axis=1)
y = data['medv']

class linear_regression:
    def __init__(self, N, M):
        # N = number of features
        # M = number of examples
        self.intercept_ = 0
        self.coef_ = np.zeros((N,), dtype=float)
        self.scaler = StandardScaler()

    def forward(self, X):
        y = np.dot(X, self.coef_) + self.intercept_
        return y

    def fit(self, X, y, lr = 0.001, epochs = 20000):
        self.scaler.fit(X)
        X = self.scaler.transform(X)

        N, M = X.shape
        error = 100
        for epoch in range(epochs):
            y_pred = self.forward(X)
            error = y_pred - y
            grad_coff = np.dot(X.T, error) / M
            grad_intercept = np.sum(error) / M
            self.coef_ -= lr * grad_coff
            self.intercept_ -= lr * grad_intercept
        return

    def predict(self, X):
        X = self.scaler.transform(X)
        return self.forward(X)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
lin_reg = LinearRegression()
lin_reg.fit(X_train, y_train)
y_pred = lin_reg.predict(X_test)
print(mean_squared_error(y_test, y_pred))

my_lin_reg = linear_regression(X_train.shape[1], X_train.shape[0])
my_lin_reg.fit(X_train, y_train)
y_pred = my_lin_reg.predict(X_test)
print(mean_squared_error(y_test, y_pred))
