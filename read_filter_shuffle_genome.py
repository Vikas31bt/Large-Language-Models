import os
import argparse
import pandas as pd
import numpy as np
import random
from Bio import SeqIO


def read_reference_genome(file_path):
    """
    Read the reference genome from a FASTA file, shuffle the chromosomes, and return the genome.

    Args:
        file_path (str): Path to the FASTA file containing the reference genome.

    Returns:
        str: Shuffled reference genome.
    """
    chromosomes = []
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            chromosomes.append(str(record.seq))  # Append each chromosome sequence

    random.shuffle(chromosomes)  # Shuffle the list of chromosomes

    reference_genome_sequence = ''.join(chromosomes)  # Concatenate the shuffled chromosomes
    return reference_genome_sequence  # Return the shuffled reference genome


def generate_input_sequences(sequence, token_length, overlap):
    """
    Generate input sequences for tokenization from a reference genome sequence,
    applying filters to select sequences with less than 5% 'N's.

    Args:
        sequence (str): Reference genome sequence.
        token_length (int): Length of each input sequence.
        overlap (int): Overlap between input sequences.

    Returns:
        list: List of input sequences.
    """
    input_sequences = []
    start = 0
    end = token_length
    while end <= len(sequence):
        input_sequence = sequence[start:end].upper()  # Convert to uppercase
        n_count = input_sequence.count('N')
        sequence_length = len(input_sequence)
        if n_count / sequence_length < 0.05:  # Less than 5% Ns
            input_sequences.append(input_sequence)
        start += token_length - overlap
        end = start + token_length

    # Shuffle the list of input sequences
    random.shuffle(input_sequences)

    return input_sequences


def save_input_sequences(input_sequences, output_folder, file_name):
    """
    Save input sequences to a file.

    Args:
        input_sequences (list): List of input sequences.
        output_folder (str): Path to the output folder.
        file_name (str): Name of the output file.

    Returns:
        None
    """
    output_file = os.path.join(output_folder, file_name.replace('.fna', '_output.txt'))
    with open(output_file, 'a') as file:  # Use 'a' mode to append to the file
        for sequence in input_sequences:
            file.write(str(sequence) + '\n')


def main():
    """
    Main function to generate input sequences for tokenization from reference genomes.

    Parses command-line arguments, reads reference genomes from FASTA files, generates input sequences,
    and saves them to an output file.

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description="Generate input sequences for tokenization from a reference genome")
    parser.add_argument("input_path", help="Path to the input folder or file containing FASTA files")
    parser.add_argument("--token_length", type=int, default=1000, help="Length of each input sequence")
    parser.add_argument("--overlap", type=int, default=50, help="Overlap between input sequences")
    parser.add_argument("--output_folder", default="tokenized_sequences", help="Path to the output folder")
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    if os.path.isdir(args.input_path):
        # Input is a directory
        input_files = [f for f in os.listdir(args.input_path) if f.endswith('.fna')]
        for file_name in input_files:
            input_file_path = os.path.join(args.input_path, file_name)
            reference_genome_sequence = read_reference_genome(input_file_path)
            input_sequences = generate_input_sequences(reference_genome_sequence, args.token_length, args.overlap)
            save_input_sequences(input_sequences, args.output_folder, file_name)
            print("Input sequences for tokenization saved from %s to %s" % (file_name, args.output_folder))
    elif os.path.isfile(args.input_path):
        # Input is a single file
        file_name = os.path.basename(args.input_path)
        reference_genome_sequence = read_reference_genome(args.input_path)
        input_sequences = generate_input_sequences(reference_genome_sequence, args.token_length, args.overlap)
        save_input_sequences(input_sequences, args.output_folder, file_name)
        print("Input sequences for tokenization saved from %s to %s" % (file_name, args.output_folder))
    else:
        print("Error: Invalid input path.")

if __name__ == "__main__":
    main()

