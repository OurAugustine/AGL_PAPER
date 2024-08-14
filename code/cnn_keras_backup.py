import pandas as pd

X_train = pd.read_csv('../datasets/digits/Digits_X_train.csv').values
y_train = pd.read_csv('../datasets/digits/Digits_y_train.csv').values
X_test  = pd.read_csv('../datasets/digits/Digits_X_test.csv').values
y_test  = pd.read_csv('../datasets/digits/Digits_y_test.csv').values

# ======================= Gradient Boosting Classifier ======================
# Call greadeint boosting classifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score

# Create the model
model = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)
# Train the model
model.fit(X_train, y_train.ravel())
# Predict the model
predicted = model.predict(X_test)
# Calculate the accuracy
accuracy = accuracy_score(y_test, predicted.reshape((-1, 1)))
print("Accuracy: ", accuracy)
# ======================= Gradient Boosting Classifier ======================

# ================================= CNN =====================================
from keras.models import Sequential
from keras.layers import Dense, Convolution2D, Flatten, MaxPooling2D

# Reshape the data to 8 * 8 * 1
X_train = X_train.reshape(X_train.shape[0], 8, 8, 1)
X_test = X_test.reshape(X_test.shape[0], 8, 8, 1)

# Create the model
model = Sequential()
# model.add(Convolution2D(32, (3, 3), activation='relu', input_shape=(8, 8, 1)))
model.add(Convolution2D(32, 3, 3, input_shape=(8, 8, 1)))
model.add(MaxPooling2D((2, 2)))
model.add(Flatten())
model.add(Dense(100, activation='relu'))
model.add(Dense(10, activation='softmax'))

# Print the model summary
print(model.summary())

# Compile the model
model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(X_train, y_train, epochs=10, batch_size=32)

# Evaluate the model
loss, accuracy = model.evaluate(X_test, y_test)
print("Accuracy: ", accuracy)
# ================================= CNN =====================================
