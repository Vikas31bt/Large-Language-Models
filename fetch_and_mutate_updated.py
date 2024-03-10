import os
import logging
import re
from math import floor, ceil
import argparse as arg

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord, write
import pandas as pd

class Genome:
    def __init__(self, fasta, index=None):
        """
        This will generate the genome instance where you can extract a given region, it will call samtools in the 
        backend to index the fasta file if the index is not there, it will also have a dict of chrom names to be 
        replaced when the fetching is happening if they are different e.g 1 vs chr1
        :param fasta: uncompressed fasta file otherwise samtools will give an error
        :param index: fasta index, if none will be created using pysam
        """
        if fasta is not None:
            self.fasta_file=os.path.abspath(fasta)
        elif not os.path.exists(fasta):
            raise FileNotFoundError("{} does not exist".format(fasta))
        else:
            raise ValueError("No fasta provided")

        if index is None and os.path.exists("{}.fai".format(fasta)):
            logging.info("Index file {}.fai found will be using this one".format(fasta))
            index=self.fasta_file+".fai"
            self.index=os.path.abspath(index)
        elif index is not None and os.path.exists(index):
            logging.info("Index file specified will be using {}".format(index))
            self.index=os.path.abspath(index)
        elif index is not None and not os.path.exists(index):
            raise FileNotFoundError("{} does not exist".format(index))
        else:
            logging.info("No index provided will be generating one now this may take couple of minutes")
            #TODO check for compression
            try:
                pysam.faidx(os.path.abspath(fasta))
                index = self.fasta_file + ".fai"
                self.index=os.path.abspath(index)
                logging.info("Indexing done the index file is {}.fai".format(os.path.abspath(fasta)))
            except Exception as e:
                logging.error("Indexing Failed with error {}".format(repr(e)))

        self.fasta=pysam.FastaFile(self.fasta_file, self.index)


    def fetch(self, chrom, start, end):
        seq=self.fasta.fetch(chrom, start, end)
        return seq

    def close(self):
        self.fasta.close()
        print("Done!")

    def closed(self):
        print(self.fasta.closed())

    
class Variant:
    def __init__(self, chrom, pos, ref, alt, genome, chrom_replace=None, regex=False):
        """
        Simple variant class, this could potentially store SVs, but I am not working on that right now.
        :param chrom: Chromosome
        :param ref: Reference sequence
        :param alt: Alternate sequence
        :param chrom_replace: Tuple for string replace of chrom names. If regex is True, it will be passed to the re module,
                             where item 0 will be replaced with item 1. For example, ("chr", "") will take "chr1" and return "1"
        :param genome: Genome class instance (from above)
        """
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.genome = genome

        if chrom_replace is not None:
            if regex:
                # This is regex
                self.chrom = str(re.sub(chrom_replace[0], chrom_replace[1], chrom))
            else:
                self.chrom = chrom.replace(chrom_replace[0], chrom_replace[1])
        else:
            self.chrom = chrom

    def mutate(self, total=10, ratio=0.5, return_as='str'):
        """
        Get the mutated sequence for a variant, only SNPs are supported.
        :param total: Total length of the sequence
        :param ratio: Ratio of the up and downstream sequences. If it's even, it will be 1 nt longer on the upstream end
        :param return_as: A string of DNA sequence or a Biopython Seq object
        """
        # Calculate lengths based on the window size
        len_to_get = total+len(self.ref)-len(self.alt)
        upstream = ceil(len_to_get * ratio)
        downstream = floor(len_to_get * ratio)

        try:
            if len(self.alt) > len(self.ref) and all(char.isalpha() and 'A' <= char <= 'Z' for char in self.alt):
                # Handle insertion or duplication
                up_seq = self.genome.fasta.fetch(self.chrom, self.pos - upstream, self.pos)
                down_seq = self.genome.fasta.fetch(self.chrom, self.pos + len(self.ref), self.pos + downstream)
            elif len(self.alt) < len(self.ref):
                # Handle deletion
                up_seq = self.genome.fasta.fetch(self.chrom, self.pos - upstream, self.pos)
                down_seq = self.genome.fasta.fetch(self.chrom, self.pos + len(self.ref), self.pos + downstream)
            elif len(self.alt) == len(self.ref) and len(self.alt) > 1:
                # Check if the reverse complement of the alternative allele is equal to the reference sequence
                reverse_complement_alt = str(Seq(self.alt).reverse_complement())
                if reverse_complement_alt == self.ref:
                    # Handle inversion
                    up_seq = self.genome.fasta.fetch(self.chrom, self.pos - upstream, self.pos)
                    down_seq = self.genome.fasta.fetch(self.chrom, self.pos + len(self.ref), self.pos + downstream)
                    # Handle translocation
                    alt_chrom, alt_pos = self.fetch_translocation_info()
                    if alt_chrom and alt_pos:
                        up_seq = self.genome.fasta.fetch(self.chrom, self.pos - upstream, self.pos)
                        down_seq = self.genome.fasta.fetch(alt_chrom, alt_pos, alt_pos + downstream)
                else:
                    raise NotImplementedError("Handling other structural variants is not implemented.")
            else:
                raise NotImplementedError("Handling other structural variants is not implemented.")
        except Exception as e:
            print("Error fetching sequences for variant: {} {} {} {}".format(self.chrom, self.pos, self.ref, self.alt))
            raise e

        # Combine the sequences to form the mutated sequence
        full_seq = up_seq + self.alt + down_seq

        # Update the seq attribute with the mutated sequence
        self.seq = full_seq

        if return_as == "str":
            mutated_seq = full_seq
        elif return_as == "seq":
            mutated_seq = Seq(full_seq)
        else:
            raise NotImplementedError("{} instance is not implemented; the current options are 'str' for string and 'seq' for Bio.Seq")

        return mutated_seq


class VariantList:
    def __init__(self, variant_list):
        """
        This is a list of Variant instances.
        
        :param variant_list: A list of variant instances. If any of the items in the list is not of Variant class, 
                             it raises ValueError.
        """
        if all(isinstance(item, Variant) for item in variant_list):
            self.variants = variant_list
        else:
            raise ValueError("Not all items in the list are of class Variant")

    def append(self, variant):
        """Append a variant to the list."""
        self.variants.append(variant)

    def write(self, path, file_type='tsv', **kwargs):
        """
        Write results to a file.
        
        :param path: Path of the file.
        :param file_type: Either 'tsv' or 'fasta'. If 'tsv', columns will be chrom, pos, ref, alt, seq.
                          If 'fasta', record names will be "chrom_pos_ref_alt".
        :param kwargs: These are passed to pd.to_csv for tsv file type.
        :return: Nothing but writes to the file.
        """
        if os.path.exists(path):
            raise OSError("{} already exists".format(path))

        if file_type == "tsv":
            data = {
                "chrom": [item.chrom for item in self.variants],
                "pos": [item.pos for item in self.variants],
                "ref": [item.ref for item in self.variants],
                "alt": [item.alt for item in self.variants],
                "seq": [str(item.seq) if hasattr(item, 'seq') else '' for item in self.variants]
            }

            df = pd.DataFrame.from_dict(data)
            df.to_csv(path, **kwargs)

        elif file_type == "fasta":
            records = []
            for item in self.variants:
                name = "_".join([str(item.chrom), str(item.pos), str(item.ref), str(item.alt)])
                if hasattr(item, 'seq'):
                    seq = Seq(str(item.seq))
                else:
                    seq = Seq('')
                records.append(SeqRecord(seq=seq, id=name, name="", description=""))
            SeqIO.write(records, path, "fasta")

        else:
            raise NotImplementedError("{} is not implemented; the current options are 'tsv' and 'fasta'".format(file_type))


class VCFReader:
    def __init__(self, vcf_file, genome_fasta_path):
        self.vcf_file = vcf_file
        self.genome = Genome(fasta=genome_fasta_path)

    def read_vcf(self):
        """
        Read variant information from a VCF file and return a list of Variant instances.
        """
        variants = []

        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip comment/header lines

                # Parse VCF line and extract relevant information
                fields = line.strip().split()

                # Ensure the line has enough fields to proceed
                if len(fields) < 5:
                    print("Skipping line: {}".format(line.strip()))
                    continue

                # Print relevant fields for debugging
                print("Chrom:", fields[0], "Pos:", fields[1], "Ref:", fields[3], "Alt:", fields[4])

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]

                # Create a Variant instance and append it to the variants list
                variant = Variant(chrom, pos, ref, alt, self.genome)
                variants.append(variant)

        return variants

    def fetch_translocation_info(self, chrom, pos):
        with open(self.vcf_file, 'r') as vcf_file:
            for line in vcf_file:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    chrom1, pos1, _, ref1, alt1, _, _, info1 = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
                    if 'SVTYPE=BND' in info1 and chrom1 == chrom and pos1 == pos:
                        # Extract information about the other end of the translocation
                        alt_chrom, alt_pos = self.extract_translocation_info(alt1)
                        return alt_chrom, alt_pos
        return None, None

    @staticmethod
    def extract_translocation_info(alt):
        # Implement the logic to extract chromosome and position from the ALT field
        # Adjust this based on how translocations are represented in your VCF file
        # Here is a simple example assuming the ALT field is in the format N[chrX:posX[N
        alt_info = alt.split('[')
        alt_chrom = alt_info[1].split(':')[0]
        alt_pos = int(alt_info[1].split(':')[1].split('[')[0])
        return alt_chrom, alt_pos



vcf_file = "/hpf/largeprojects/ccmbio/vikas_kaushik/1112_JR_J140-wham_annotated.vcf"
genome_fasta_path = "/hpf/largeprojects/ccmbio/vikas_kaushik/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

# Create a Genome object
genome = Genome(fasta=genome_fasta_path)

# Create a VCFReader object
vcf_reader = VCFReader(vcf_file, genome_fasta_path)

# Read variants from the VCF file
variants = vcf_reader.read_vcf()

# Create a VariantList object and append the variants
variant_list = VariantList(variants)

# Specify the output file path for the FASTA format
output_fasta_path = "/hpf/largeprojects/ccmbio/vikas_kaushik/output.fasta"

# Write the variant list to a file in FASTA format
variant_list.write(path=output_fasta_path, file_type='fasta')

# Print a message
print("Variant list written to {}".format(output_fasta_path))

