#!/bin/bash
#SBATCH -N 1 -c 1
#SBATCH --mem=200g
#SBATCH -t 1-0
#SBATCH --tmp 20g
#SBATCH -o tokenization_slurm_output.log

# Define variables
INPUT_FOLDER="."  # Default to current directory
TOKEN_LENGTH=6000
OVERLAP=100
OUTPUT_FOLDER="output"

#Python Installation:
module load python/3.6.8

# Biopython Installation for Python 3.6.8
python -m pip install --user biopython==1.77

# Check if Biopython installation was successful
if [ $? -ne 0 ]; then
    echo "Error: Biopython installation failed."
    exit 1
fi

# Loop over all .fna files in the input folder
for file in "$INPUT_FOLDER"/*.fna; do
    # Get the filename without extension
    filename=$(basename "$file" .fna)
    # Run input_for_tokenization.py for each file
    echo "Running tokenization process for file: $filename"
    python input_for_tokenization.py "$file" --token_length "$TOKEN_LENGTH" --overlap "$OVERLAP" --output_folder "$OUTPUT_FOLDER"
    # Check if the script executed successfully
    if [ $? -ne 0 ]; then
        echo "Error: Tokenization process failed for file: $filename"
        exit 1
    fi
    echo "Tokenization process completed successfully for file: $filename"
done

echo "Tokenization process completed successfully for all files."

