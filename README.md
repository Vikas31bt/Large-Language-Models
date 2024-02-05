# DNA Sequence Generator

## Project Overview

The DNA Sequence Generator project aims to generate individualized DNA sequences from existing patient data. This includes handling various types of genetic variations, such as SNPs, indels, and structural variants. The goal is to create a holistic representation of an individual's DNA and leverage large, annotated resources for comprehensive understanding.

## Project Goals

1. Generate individualized DNA sequences considering all variations for a specific individual.
2. Benefit from large, annotated resources to re-train models for a more holistic view.

## Types of Genetic Variations

### 1. SNPs/Indels
- Extract sequences from the FASTA file.
- Replace existing sequences with variants.
- Consider homozygous and heterozygous cases.
- Pre-filter variants to reduce complexity.

### 2. Structural Variants
- Handle insertions, deletions, duplications, and translocations.
- Address SNPs/Indels inside structural variants.
- Make decisions on implementing changes.

### 3. Genome Versions
- Support different genome versions (hg37, hg38, t2t consortium sequence).

### 4. Homozygotes/Heterozygotes
- Handle both homozygotes and heterozygotes.
- Limit the number of outputs for heterozygotes in larger sequence intervals.

## Script Requirements

### Input
- VCF file(s)
- FASTA file
- (Optional) BED file of intervals for multiple regions

### Output
- Single multifasta file
- (Initial iteration) Number of possible outputs for a single interval
- (Final iteration) Multifasta files for specified regions from BED file

### Runtime Considerations
- Fast execution for processing a set of VCFs in a matter of hours.
- Utilize multiprocessing for bedfile options to enhance speed.
