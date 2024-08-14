from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import pandas as pd
data = pd.read_csv('../datasets/Boston.csv')

X = data.drop('medv', axis=1)
y = data['medv']

import torch # install pytorch

class my_linear_regression:
    def __init__(self, X_shape, lr=0.01, epochs=10000):
        self.lr = lr
        self.epochs = epochs
        self.coef_ = torch.zeros(X_shape, 1, requires_grad=True)
        self.intercept_ = torch.zeros(1, requires_grad=True)

    def forward(self, X):
        y_hat = X @ self.coef_ + self.intercept_
        return y_hat

    def fit(self, X, y):
        for epoch in range(self.epochs):
            y_hat = self.forward(X)
            loss = ((y_hat - y)**2).mean()
            loss.backward()
            with torch.no_grad():
                self.coef_ -= self.lr * self.coef_.grad
                self.intercept_ -= self.lr * self.intercept_.grad
                '''
                The reason to zero the gradients is to ensure that we only have 
                the gradients for the current batch of data when we perform the update step. 
                If we didn't zero the gradients, we would be accumulating gradients over time 
                and our update step would be using incorrect information.
                '''
                self.coef_.grad.zero_() # Zero the gradients after updating
                self.intercept_.grad.zero_()

        return self


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

lin_reg = LinearRegression()
lin_reg.fit(X_train, y_train)
y_pred = lin_reg.predict(X_test)
print(mean_squared_error(y_test, y_pred))

my_lin_reg = my_linear_regression(X_train.shape[1])
y_train = y_train.values.reshape(-1, 1)
my_lin_reg.fit(torch.tensor(X_train, dtype=torch.float32), torch.tensor(y_train, dtype=torch.float32))

y_pred = my_lin_reg.forward(torch.tensor(X_test, dtype=torch.float32)).detach().numpy()
print(mean_squared_error(y_test, y_pred))
