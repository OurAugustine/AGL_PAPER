import pandas as pd    # to load dataset
import nltk
from nltk.corpus import stopwords   # to get a collection of stopwords

data = pd.read_csv('../datasets/IMDB.csv')

custom_path = '../datasets/'

# Append your custom path to the NLTK data path
nltk.data.path.append(custom_path)

nltk.download('stopwords', download_dir=custom_path)
english_stops = set(stopwords.words('english'))

x_data = data['review']       # Reviews/Input
y_data = data['sentiment']    # Sentiment/Output
# PRE-PROCESS REVIEW
x_data = x_data.replace({'<.*?>': ''}, regex = True)          # remove html tag
x_data = x_data.replace({'[^A-Za-z]': ' '}, regex = True)     # remove non alphabet
x_data = x_data.apply(lambda review: [w for w in review.split() if w not in english_stops])  # remove stop words
x_data = x_data.apply(lambda review: [w.lower() for w in review])   # lower case

from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
tokenizer = Tokenizer(num_words=10000)    # num_words is the number of words to keep based on word frequency
tokenizer.fit_on_texts(x_data)            # fit tokenizer to our training text data

# retrieve the word index
word_index = tokenizer.word_index

x_data = tokenizer.texts_to_sequences(x_data)  # convert our text data to sequence of numbers

import numpy as np
from keras.models import Sequential
from keras.layers import Embedding, SimpleRNN, LSTM, Dense, GRU
from sklearn.model_selection import train_test_split
# from keras.utils import np_utils

# Pad sequences to ensure uniform input size
max_length = 100  # You can choose a different length
x_data = pad_sequences(x_data, maxlen=max_length)

# Convert sentiments to numerical labels
y_data = np.where(y_data == 'positive', 1, 0)

# Split data into training and test sets
x_train, x_test, y_train, y_test = train_test_split(x_data, y_data, test_size=0.2, random_state=42)

# Build the RNN model
model = Sequential()
model.add(Embedding(input_dim=10000, output_dim=128, input_length=max_length))  # Embedding layer
model.add(SimpleRNN(64))  # You can also try using LSTM or GRU layers
# model.add(SimpleRNN(64))  # You can also try using LSTM or GRU layers
# model.add(LSTM(64))  # You can also try using LSTM or GRU layers
# model.add(GRU(64))  # You can also try using LSTM or GRU layers
model.add(Dense(1, activation='sigmoid'))  # Output layer
# print(model.summary())

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(x_train, y_train, epochs=5, batch_size=64, validation_split=0.2)

# Evaluate the model on the test set
test_loss, test_acc = model.evaluate(x_test, y_test)
print(f"Test Accuracy: {test_acc}")

