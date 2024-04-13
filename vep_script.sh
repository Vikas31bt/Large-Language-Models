#!/bin/bash
#SBATCH -N 1 -c 12
#SBATCH --mem=50g 
#SBATCH -t 4-0 
#SBATCH --tmp 20g
#SBATCH -o slurm_output.log

#Call the vep module
module load vep/94

# Specify the directory containing  VCF files(SNP)
vcf_directory="/hpf/largeprojects/ccmbio/vikas_kaushik/vcfs/snp"

# Specify the output directory for annotated VCF files
output_directory="/hpf/largeprojects/ccmbio/vikas_kaushik/vcfs/snp"

# VEP job options
vep_options="--force_overwrite --no_stats --cache --offline --assembly GRCh37 --dir /hpf/tools/centos6/vep/cache94/ --everything --fork 12"

# Iterate over VCF files in the directory
for vcf_file in "$vcf_directory"/*.vcf.gz; do
    # Extract the file name without extension
    file_name=$(basename "$vcf_file" .vcf.gz)

    # Run VEP on each VCF file and save the annotated output
    vep -i "$vcf_file" $vep_options -o "$output_directory/${file_name}_annotated.vcf" --vcf
done

