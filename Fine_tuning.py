from datasets import load_dataset
from transformers import AutoTokenizer, AutoModelForSequenceClassification, TrainingArguments, Trainer, DataCollatorForLanguageModeling, DistilBertForMaskedLM, DistilBertTokenizer, DistilBertConfig
from tokenizers.implementations import BertWordPieceTokenizer
from tokenizers.processors import BertProcessing
import pandas as pd
import torch
import numpy as np
import evaluate


# Load an existing dataset for sentiment analysis
imdb = load_dataset("imdb")

# Define label mappings
id2label = {0: "NEGATIVE", 1: "POSITIVE"}
label2id = {"NEGATIVE": 0, "POSITIVE": 1}

# Load tokenizer and model for sequence classification
tokenizer = AutoTokenizer.from_pretrained("distilbert/distilbert-base-uncased", padding=True, truncation=True, max_length=512)
model = AutoModelForSequenceClassification.from_pretrained("distilbert/distilbert-base-uncased", num_labels=2, id2label=id2label, label2id=label2id)
model

# Define preprocessing function
def preprocess_function(examples):
    return tokenizer(examples["text"], truncation=True, padding=True, max_length=512)

list(imdb.keys())

# Tokenize the IMDb dataset
tokenized_imdb = imdb.map(preprocess_function, batched=True)

# Define evaluation metric
accuracy = evaluate.load("accuracy")

# Define compute_metrics function
def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    predictions = np.argmax(predictions, axis=1)
    return accuracy.compute(predictions=predictions, references=labels)

# Define training arguments
training_args = TrainingArguments(
    output_dir="my_awesome_model",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=16,
    num_train_epochs=2,
    weight_decay=0.01,
    evaluation_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
)

# Initialize Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_imdb["train"],
    eval_dataset=tokenized_imdb["test"],
    tokenizer=tokenizer,
    compute_metrics=compute_metrics
)

# Fine tune the model
trainer.train()

# Evaluate the model on the training dataset
train_eval_result = trainer.evaluate(eval_dataset=tokenized_imdb["train"])
print("Training evaluation results:", train_eval_result)

# Evaluate the model on the test dataset
test_eval_result = trainer.evaluate(eval_dataset=tokenized_imdb["test"])
print("Test evaluation results:", test_eval_result)
