#!/bin/bash
#SBATCH -N 1 -c 1
#SBATCH --mem=200g
#SBATCH -t 4-0
#SBATCH --tmp 20g
#SBATCH -o snp_sv_slurm_output.log

# Biopython Installation for Python 3.6.8
python -m pip install --user biopython==1.77


# PyVCF Installation for Python 3.6.8
python -m pip install --user PyVCF

# Run Genome_reconstruction.py
echo "Running genome reconstruction with Genome_reconstruction.py:"
python snp_sv.py

