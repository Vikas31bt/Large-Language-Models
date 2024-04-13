import numpy as np
from transformers import AutoTokenizer, AutoModelForMaskedLM, TrainingArguments, Trainer, DataCollatorForLanguageModeling
from datasets import load_dataset, Dataset

# Load tokenizer and model
tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-250m-multi-species", trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-250m-multi-species", trust_remote_code=True)

# Load dataset
dataset = load_dataset("text", data_files="GRCh38_latest_genomic_output.txt")
train_data = Dataset.from_dict(dataset["train"][:100])  # Adjust the number of samples as needed
eval_data = Dataset.from_dict(dataset["train"][100:120])  # Adjust the number of samples as needed for evaluation

# Preprocess function
preprocess_function = lambda examples: tokenizer(examples["text"], truncation=True, padding=True, max_length=1000)
train_dataset = train_data.map(preprocess_function)
eval_dataset = eval_data.map(preprocess_function)

# Compute metrics function
def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    predictions = predictions.argmax(axis=-1) 
    print("Predictions shape:", predictions.shape)
    print("Labels shape:", labels.shape)
    return {"accuracy": (predictions == labels).mean().item()}

# Data collator
data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=True, mlm_probability=0.15)

# Training arguments
training_args = TrainingArguments(
    output_dir="my_awesome_model",
    learning_rate=2e-5,
    do_eval=True,
    do_train=True,
    per_device_train_batch_size=1,
    num_train_epochs=2,
    weight_decay=0.01,
    evaluation_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
)

# Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    data_collator=data_collator,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset,
    compute_metrics=compute_metrics,
)

# Start training
trainer.train()

# Print the evaluation results
eval_results = trainer.evaluate()
print("Evaluation results:", eval_results)
