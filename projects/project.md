## Final Project

**NOTE**: Only select one of the following projects to be your final project. Or, you can do your own project, but you must get approval from the instructor first.<br>
The final project should be turned in as a Jupyter notebook. It should be well-documented and organized. It should be clear how to run the code and what the code does. The code should be well-written and efficient. The project should be your own work. You can discuss the project with others, but the code and write-up should be your own work. You should cite any sources you use. You should also include a list of references at the end of the project. The project should be turned in on your GitLab repository. The project is due on the last day of class. 
The penality of late projects is 10% per day.

### Project 1
There are 29 datasets generated from real applications. Choose one dataset as the target traing/test set and use the methods covered in this semster to study the dataset. 

#### Datasets
``` text
1. Wine Quality, 2. Bike Data, 3. Kobe Bryant, 4. West Nile virus, 5. Year Music prediction, 
6. Income classification, 7. Bank marketing, 8. Credit card clients, 9. Drug consumption, 
10. Geographical original of music, 11. HTRU2, 12. Letter recognition, 13. Mice protein expression, 
14. Occupancy detection, 15. Online news popularity, 16. Ozone level detection, 
17. Phishing websites, 18. Polish companies bankruptcy, 19. Spambase, 20. Taxi Newyork, 
21. Allstate, 22. House prices, 23. Shelter Animal outcomes, 24. Weather prediction, 
25. Binding energy, 26. Toxicity, 27. Slvation, 28. Coil-20, 29. USPS
```

#### Tasks
1. After the 1-on-1 meeting with the instructor, you need to know what type of problem you are going to solve, regression or classification or clustering.
   - Classification methods: Logistic Regression, Random Forest, Gradient Boosting Decision Tree, SVM, KNN
   - Regression methods: Linear Regression (Ridge Regression, Lasso Regression), SVM, Random Forest, Gradient Boosting Decision Tree
2. Check your dataset and think about each feature. What if the feature is word or text? What if the feature is image? What if the feature is time series?
3. To start the project, under your gitlab repository folder `projects`, create the following files/folders:
   - `projects/readme.md`: A file to describe your project, your target dataset, and your plan. In the end, you should also put your results in this markdown file.
   - `projects/data/`: A folder to put your dataset.
   - `projects/src/`: A folder to put your code. The code loads the dataset, preprocesses the dataset, trains the model, and analyze the models.
4. The project follows the following steps:
   - Load the dataset: you need to load the dataset and check the basic information of the dataset.
   - Preprocess the dataset: you need to preprocess the dataset. The preprocessing includes checking word or text features, and time series.
   - Time series: if your dataset is time series, we will cover the Long Short-Term Memory (LSTM) model. You need to use the LSTM model to predict the future. This is optional.
   - Preprocess the dataset also includes splitting the training dataset by 10-fold cross-validation.
   - Train the model: you need to train the model using the methods we covered in this semester. You need to compare the results of different methods. For each methods, you need to employ the 10-fold cross-validation on the training set. Using the grid search to find the optimized parameters (for-loop). Once the parameters are fixed, evaluate the model on the test set and compare the methods with other methods.
   - Analyze the model: you need to analyze the model. You need to check the feature importance, the confusion matrix, the ROC curve, the precision-recall curve, and the accuracy or RMSE.
5. Write a short report to explain your project in `projects/readme.md`. The report should include the following:
   - Introduction: describe the dataset and the problem you are going to solve.
   - Preprocessing: describe the preprocessing of the dataset.
   - Methods: describe the methods you used in the project.
   - Results: describe the results of the project.
   - Conclusion: describe the conclusion of the project.

### Project 2 
There are some topics that we don't cover in this semester. You can choose one of the following topics and do a project on it. 
The requirement of this project is a little bit different since you don't work on a dataset specifically. 
You need to write a short report and code to explain the topic. The topic includes:
1. Generative Adversarial Networks (GANs)
2. Recurrent Neural Networks (RNNs)
3. Boltzmann Machines
4. ResNet
5. Graph Attention Networks (GATs) 
6. Transformer
7. Reinfocement Learning
