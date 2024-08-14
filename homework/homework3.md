## Homework 3
This homework is due on 04/15/2024 by 11:59 PM. Homework 3 is about resampling methods and decision trees.

### Question 1 
The bootstrap is one part of the bagging algorithm. 
At its core, bootstrap involves repeatedly resampling, with replacement, from the observed data set and recalculating the statistic of interest on each resample. This process builds an empirical distribution of the statistic, providing insight into its variability without requiring stringent assumptions about the underlying population distribution from which the sample was drawn.

The bootstrap is a statistical technique that allows you to estimate the distribution of a statistic (such as the mean, median, or standard deviation) of a population from a sample of data. Its beauty lies in its simplicity and versatility, making it applicable across a wide range of situations in statistical inference, where understanding the variability of an estimator is crucial. Here's how it works and why it's so powerful:

### Step-by-Step Explanation

1. **Original Sample**: You start with an original sample of \(n\) observations, which serves as a proxy for the larger population. This sample is presumed to contain enough information to model the variability of the population.

2. **Bootstrap Samples**: From this original sample, you generate a large number of bootstrap samples. Each bootstrap sample is the same size as the original and is created by randomly selecting observations from the original data with replacement. This means any observation can appear more than once or not at all in each bootstrap sample.

3. **Calculate Statistics**: For each bootstrap sample, calculate the statistic of interest (e.g., mean, variance). Because you're working with a high number of bootstrap samples, you'll get a wide array of values for the statistic.

4. **Estimate Distribution**: The collection of calculated statistics from all the bootstrap samples forms an empirical distribution. This distribution approximates the sampling distribution of the statistic if you were able to draw infinite samples from the population.

5. **Quantify Uncertainty**: With the empirical distribution, you can now quantify the uncertainty of the estimator. For example, you can estimate the standard error, confidence intervals, or prediction intervals of the statistic, providing a robust measure of its reliability.

### Questions
We will now derive the probability that a given observation is part of a boostrap sample. Suppose that we obtain a boostrap sample from a set of $n$ observations. 

1. What is the probability that the first bootstrap observation is not the $j$th observation from the original sample? Justify your answer.

2. What is the probability that the second bootstrap observation is not the $j$th observation from the original sample? Justify your answer.

3. Argue that the probability that the $j$th observation is not in the bootstrap sample is $(1-1/n)^n$.

4. When $n = 5$, what is the probability that the $j$th observation is in the bootstrap sample? 

5. When $n = 100$, what is the probability that the $j$th observation is in the bootstrap sample?

6. When $n = 10000$, what is the probability that the $j$th observation is in the bootstrap sample?

7. Create a plot that shows how the probability that the $j$th observation is in the bootstrap sample changes as $n$ increases from 1 to 10000. Comment on what you observe. (Need coding)

## Question 2 
a. Explain how $k$-fold cross-validation is implemented.

b. What are the advantages and disadvantages of $k$-fold cross-validation relative to LOOCV?

- Leave-one-out cross-validation (LOOCV) is a special case of $k$-fold cross-validation where $k$ is equal to the number of observations in the data set.
The learning method is fit on the $n-1$ training observations and a prediction is made for the excluded observation. This process is repeated $n$ times, with each of the $n$ observations used exactly once as the validation data. The $n$ validation errors are then averaged to obtain the LOOCV estimate of the test error.

## Question 3 
Consider the Gini index, classification error, and entropy in a simple classification setting with two classes. 
Create a single plot that displays each of these quantities as a function of \(\hat{p}_{mk}\). The x-axis should display \(\hat{p}_{mk}\), ranging from 0 to 1, and the y-axis should display the value of the Gini index, classification error, and entropy.

Hint: In a setting with two classes, \(\hat{p}_{m1} = 1 - \hat{p}_{m2}\). The Gini index is given by \(2\hat{p}_{m1}\hat{p}_{m2}\), the classification error is given by \(1 - \max(\hat{p}_{mk}, 1 - \hat{p}_{mk})\), and the entropy is given by \(-\hat{p}_{mk}\log_2(\hat{p}_{mk}) - (1 - \hat{p}_{mk})\log_2(1 - \hat{p}_{mk})\).

## Question 4 
Suppose we produce ten bootstrapped samples from a data set containing red and green classes. We then apply a classification tree to each bootstrapped sample and, for a specific value of $X$, produce 10 estimates of $P(\text{Class} = \text{Green} | X)$

\[0.1, 0.15, 0.2, 0.2, 0.55, 0.6, 0.6, 0.65, 0.7, 0.75\]

There are two common ways to combine these results together to produce a single class prediction. One is the majority vote approach, where the most commonly occurring class is predicted. The other is the average probability approach, where the average of the estimated probabilities is taken.
In this example, what is the final prediction under each of these two approaches?

## Question 5 
This is a semi report of your project. You need to write a brief summary of your project and the steps you have taken so far.

If you choose the project 1, you need to write 
1. how you have split the data, 
2. what models you have used, 
3. what are the results of your models,
4. how you process the grid search,
5. and what are the best parameters you have found.

Also, you need to write what you are planning to do next.

If you choose the project 2, you need to write how you understand the method, what you have done so far, and what you are planning to do next.
