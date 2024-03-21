import os
import argparse
import random
from Bio import SeqIO


def read_reference_genome(file_path):
    """
    Read the reference genome

    Args:
    file_path (str): Path to the FASTA file containing the reference genome.

    Returns:
    str: Reference genome.
    """
    reference_genome_sequence = ""
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            reference_genome_sequence += str(record.seq)  # Append each chromosome sequence

    return reference_genome_sequence  # Return the reference genome


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


def save_input_sequences(input_sequences, output_file):
    """
    Save input sequences to a file.

    Args:
        input_sequences (list): List of input sequences.
        output_file (str): Path to the output file.

    Returns:
        None
    """
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
    parser.add_argument("input_folder", help="Path to the input folder containing FASTA files")
    parser.add_argument("--token_length", type=int, default=1000, help="Length of each input sequence")
    parser.add_argument("--overlap", type=int, default=50, help="Overlap between input sequences")
    parser.add_argument("--output_folder", default="tokenized_sequences", help="Path to the output folder")
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    # Iterate over files in the input folder
    for file_name in os.listdir(args.input_folder):
        if file_name.endswith('.fna'):
            input_file_path = os.path.join(args.input_folder, file_name)
            # Read reference genome sequence from input FASTA file
            reference_genome_sequence = read_reference_genome(input_file_path)
            # Generate input sequences for tokenization
            input_sequences = generate_input_sequences(reference_genome_sequence, args.token_length, args.overlap)
            # Construct the output file path
            output_file = os.path.join(args.output_folder, file_name.replace('.fna', '_output.txt'))
            # Save input sequences to output file
            save_input_sequences(input_sequences, output_file)
            print("Input sequences for tokenization saved from %s to %s" % (file_name, output_file))

if __name__ == "__main__":
    main()

