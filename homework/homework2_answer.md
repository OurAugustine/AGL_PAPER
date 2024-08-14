## Homework 2

### Question 1 Answer
(a) Compute the Euclidean distance between each observation and the test point, $X_1 = X_2 = X_3 = 0$.

Obs 1 has Euclidean Distance $\sqrt{(0 - 0)^2 + (3 - 0)^2 + (0 - 0)^2} = 3$.

Obs 2 has Euclidean Distance $\sqrt{(2 - 0)^2 + (0 - 0)^2 + (0 - 0)^2} = 2$.

Obs 3 has Euclidean Distance $\sqrt{(0 - 0)^2 + (1 - 0)^2 + (3 - 0)^2}  \approx 3.16$.

Obs 4 has Euclidean Distance $\sqrt{(0 - 0)^2 + (1 - 0)^2 + (2 - 0)^2}  \approx 2.24$.

Obs 5 has Euclidean Distance $\sqrt{(-1 - 0)^2 + (0 - 0)^2 + (1 - 0)^2} \approx 1.41$.

Obs 6 has Euclidean Distance $\sqrt{(1 - 0)^2 + (1 - 0)^2 + (1 - 0)^2}  \approx1.73$.

(b) What is our prediction with K = 1? Why?

The nearest neighbor to test point (0, 0, 0) is Obs 5 (-1, 0, 1) with euclidean distance ~1.41. Since Obs 5 was Green, we predict (K = 1) that the test point will also be Green.

(c) What is our prediction with K = 3? Why?

The nearest three neighbors to test point (0, 0, 0) are Obs 5, Obs 6 (with distance ~1.73), and Obs 2 (with distance 2). Since Obs 5 was Green, Obs 6 was Red, and Obs 2 was Red, we predict (K = 3) the test point will be the majority -- Red.

(d) If the Bayes decision boundary in this problem is highly nonlinear, then would we expect the best value for K to be large or small? Why?

With the increase of $K$, the KNN method becomes less flexible thus with lower variance and higher bias, which is more difficult to match the nonlinear Bayes decision boundary.  So the best value for K would be small, but not too small, which might cause overfitting.

### Question 4 Answer
$\hat{y}_i=\hat{\beta} x_i = \displaystyle \frac{\sum\limits_{j=1}^n x_j y_j}{\sum\limits_{j=1}^n x_j^2}x_i = \sum\limits_{j=1}^n \left(\frac{x_i x_j}{\sum\limits_{j=1}^n x_j^2} \right) y_j = \sum\limits_{j=1}^n  a_j y_j$. <br>
Note $x_i$ and $\sum\limits_{j=1}^n x_j^2$ does not depend on $j$ thus can get into the sum with index $j$.
