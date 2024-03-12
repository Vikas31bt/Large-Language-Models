import os
import argparse
from Bio import SeqIO

def read_reference_genome(file_path):
    with open(file_path, 'r') as file:
        reference_genome_sequence = SeqIO.read(file, 'fasta').seq
    return reference_genome_sequence

def generate_input_sequences(sequence, token_length, overlap):
    input_sequences = []
    start = 0
    end = token_length
    while end <= len(sequence):
        input_sequence = sequence[start:end]
        input_sequences.append(input_sequence)
        start += token_length - overlap
        end = start + token_length
    return input_sequences

def save_input_sequences(input_sequences, output_folder):
    for i, sequence in enumerate(input_sequences):
        output_file_path = os.path.join(output_folder, "input_sequence_%d.txt" % i)
        with open(output_file_path, 'w') as file:
            file.write(str(sequence) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Generate input sequences for tokenization from a reference genome")
    parser.add_argument("input_folder", help="Path to the input folder containing FASTA files")
    parser.add_argument("--token_length", type=int, default=1000, help="Length of each input sequence")
    parser.add_argument("--overlap", type=int, default=50, help="Overlap between input sequences")
    parser.add_argument("--output_folder", default="tokenized_sequences", help="Path to the output folder")
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # Iterate over files in the input folder
    for file_name in os.listdir(args.input_folder):
        if file_name.endswith('.fasta'):
            input_file_path = os.path.join(args.input_folder, file_name)
            # Read reference genome sequence from input FASTA file
            reference_genome_sequence = read_reference_genome(input_file_path)
            # Generate input sequences for tokenization
            input_sequences = generate_input_sequences(reference_genome_sequence, args.token_length, args.overlap)
            # Save input sequences to output folder
            save_input_sequences(input_sequences, args.output_folder)
            print("Input sequences for tokenization saved from %s to %s" % (file_name, args.output_folder))

if __name__ == "__main__":
    main()

