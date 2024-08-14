## Homework 1 
The first homework is about the basic concepts of machine learning. Questions are selected from the book *An introduction to statistical learning*.
The homework is due on 12th of February.
* For undergraduate students, answer 1-3 and 5. Question 4 is optional.
* For graduate students, answer all questions.

### Question 1 
For the following questions, indicate whether we would generally expect the performance of a flexible statistical learning method to be better or worse than an inflexible method. Justify your answer.
The predictor is also known as an independent variable, feature, or input variable. The response is also known as a dependent variable or output variable.
1. The sample size $n$ is extremely large, and the number of predictors $p$ is small.
2. The number of predictors $p$ is extremely large, and the number of observations $n$ is small.
3. The relationship between the predictors and response is highly non-linear.
4. The variance of the error terms, i.e. $\sigma^2 = Var(\epsilon)$, is extremely high.

### Question 2 
Explain whether each scenario is a classification or regression problem, and indicate whether we are most interested in inference or prediction. Finally, provide $n$ and $p$.
1. We collect a set of data on the top 500 firms in the US. For each firm we record profit, number of employees, industry and the CEO salary. We are interested in understanding which factors affect CEO salary.
2. We are considering launching a new product and wish to know whether it will be a *success* or a *failure*. We collect data on 20 similar products that were previously launched. For each product we have recorded whether it was a success or failure, price charged for the product, marketing budget, competition price, and ten other variables.
3. We are interested in predicting the % change in the US dollar in relation to the weekly changes in the world stock markets. Hence we collect weekly data for all of 2012. For each week we record the % change in the dollar, the % change in the US market, the % change in the British market, and the % change in the German market.

### Question 3
You will now think of some real-life applications for machine/statistical learning.
1. Describe three real-life applications in which classification might be useful. Describe the response, as well as the predictors. Is the goal of each application inference or prediction? Explain your answer.
2. Describe three real-life applications in which regression might be useful. Describe the response, as well as the predictors. Is the goal of each application inference or prediction? Explain your answer.
3. Describe three real-life applications in which cluster analysis might be useful.

### Question 4 
Describe the differences between a parametric and a non-parametric statistical learning approach. What are the advantages of a parametric approach to regression or classification (as opposed to a non-parametric approach)? What are its disadvantages?

### Question 5 
This exercies involves the `boston` housing dataset. Have a jupyter notebook ready to answer the following questions.
The `boston` dataset is in the folder `datasets` of this repository and describes housing values and other information about Boston census tracts.
1. Load the dataset and explore the variables by `pandas`
2. Use the `matplotlab` to plot the two variables/predictors `indus` vs `nox` by the following code:
```python
import matplotlib.pyplot as plt
plt.scatter(data.loc[:,'indus'], data.loc[:,'nox']) # data = pd.read_csv('datasets/boston.csv')
plt.show()
```
3. Use the `pandas` to list all the variables/predictors (columns).
4. After exploring the dataset, select the `medv` as the response/label 'y' and the rest of the variables as the predictors `X`.
5. Split the dataset into a training set and a testing set. The training set should be 70% of the original dataset and the test set is the rest 30%.
6. Train a linear regression model on the training set and evaluate its performance on the test set.
```python
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
lin_reg = LinearRegression()
lin_reg.fit(X_train, y_train)
y_pred = lin_reg.predict(X_test)
mean_squared_error(y_test, y_pred)
```
