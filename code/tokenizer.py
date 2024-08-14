from keras.preprocessing.text import Tokenizer

# Example text data
x_data = ['This is a sample text, with some words.']

# Initialize a tokenizer and fit on the text data
tokenizer = Tokenizer()
tokenizer.fit_on_texts(x_data)

# Print original length of the first item
print("Original length:", len(x_data[0].split()))

# Convert text data to sequence of numbers
x_data_seq = tokenizer.texts_to_sequences(x_data)

# Print length after tokenization
print("Length after tokenization:", len(x_data_seq[0]))


from keras.preprocessing.sequence import pad_sequences

sequences = [
    [1, 2, 3],
    [4, 5, 6, 7, 8, 9],
    [10, 11]
]

padded_sequences = pad_sequences(sequences, maxlen=5)
print(padded_sequences)

