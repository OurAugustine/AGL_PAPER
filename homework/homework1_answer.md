## Homework 1
Since questions are selected from the book *An introduction to statistical learning*, the figures in the following answers are also from the book.

### Question 1 Answer
Some <u>key</u> points note before answering questions. <br>
- When we measure the quality of fitting, we use MSE, which consists of variance of the model $\text{Var}(\hat{f})$, bias of the model $\text{Bias}(\hat{f})$, and variance of noise $\text{Var}(\epsilon)$.
- Lower the <u>test</u> MSE is what we really concern. It is a $U$-shaped pattern with the increase of flexibility (e.g. Figure 2.9), resulted from the Bias-Variance trade-off.
- As opposed to the inflexible method, the flexible method is more expensive to implement, is more difficult to provide statistic inferences, has higher variance $\text{Var}(\hat{f})$ and lower bias $\text{Bias}(\hat{f})$, fits the (training) data pattern better and will cause overfitting with the increase of flexibility, 
- Which method is preferred also depends on the pattern of data $f(x)$, the noise, and the available observations $n$. The optimal flexibility level differs considerably among data (e.g. Figure 2.12).

(a) Flexible method is preferred since it in general fits the data better than the inflexible method and the large $n$ reduces the risk of overfitting. <br>
(b) This is the situation both methods will encounter difficulties. Inflexible method is preferred because when $n$ is small, the flexible method will more likely cause overfitting to the pattern of the noise and it is more expensive to implement. Note: flexible method (e.g. non-parametric) in general requires more observations. <br>
(c) Being highly nonlinear indicates the complexity of the model thus we need flexible method to pick up the pattern to the reduce the bias. (compare figures 2.10 and 2.11)<br>
(d) This is the situation both methods will encounter difficulties. An inflexible method is preferred since it is easier to implement and less likely to overfit to the noise from randomness with high variance.<br>

### Question 2 Answer
(a) This is a regression problem, because CEO salary is a continuous variable. We are interested in inference, because we want to know how the predictors impact the response CEO salary. n = 500 (500 firms provide our data points) and p = 3 (profit, number of employees, and industry). <br>
(b) This is a classification problem, as we are trying to classify success vs. failure. This is also a prediction problem because we want to predict success vs. failure with less concern for understanding the underlying relationships. n = 20 (20 similar products previously launched) and p = 13 (price, marketing budget, competition price, and ten other variables).<br>
(c) This is a regression problem because \%change in the US dollar is a quantitative dependent variable. This is also a prediction problem because we are interested in predicting the \% change. n = 52 (weekly data over 2012) and p = 3 (\% change in US, \% change in British, \% change in German).<br>

### Question 3 Answer
Describe the differences between a parametric and a non-parametric statistical learning approach. What are the advantages of a parametric approach to regression or classification (as opposed to a non-parametric approach)? What are its disadvantages?

### Question 4 Answer
A non-parametric approach is better at fitting highly nonlinear patterns, because less assumptions are made about the structure of the data. However, it carries the risk of overfitting thus we need to optimize the level of flexibility. In addition, nonparametric method is preferred when observation is abundant and/or the irreducible error is low. 
Parametric method is easier to implement and can interpret the data using the parameters of the model and provide inference. Because of the bias-variance tradeoff, neither approach is absolutely better than the other. 

