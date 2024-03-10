import os
import logging
import re
from math import floor, ceil
import argparse as arg

import pysam
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
        simple variant class, this could potentially store svs but I am not working on that right now.
        :param chrom: chromosome
        :param ref: reference sequence
        :param alt: alternate sequence
        :param chrom_replace: tuple string replace for chrom names if regex will be passed to re module item 0 will be
        replaced with item 1 so ("chr", ""), will take "chr1" and return "1"
        :param genome: Genome class instance (from above)
        """
        self.pos=pos
        self.ref=ref
        self.alt=alt
        self.genome=genome

        if chrom_replace is not None:
            if regex:
                #this is regex
                self.chrom=str(re.sub(chrom_replace[0], chrom_replace[1]))
            else:
                self.chrom=chrom.replace(chrom_replace[0], chrom_replace[1])
        else:
            self.chrom=chrom


    def mutate(self, total=1000, ratio=0.5, return_as='str'):
        """
        get the mutated sequence for a variant, only snps are supported
        :param total: total length of the sequence
        :param ratio: ratio of the up and down stream sequences, if it's even it will be 1 nt longer on the upstream end
        :return_as: a string of dna sequence or a biopython str object
        """
        len_to_get=total+len(self.ref)-len(self.alt)
        upstream=ceil(len_to_get*ratio)
        downstream=floor(len_to_get*ratio)

        #TODO this is a slooooooow implementation, I would need to find a way to store chrom sequences as a torch tensor and
        # offload to gpu
        up_seq=self.genome.fasta.fetch(self.chrom, self.pos-upstream, self.pos)
        down_seq = self.genome.fasta.fetch(self.chrom, self.pos+len(self.ref), self.pos+downstream)

        full_seq=up_seq+self.alt+down_seq

        if return_as=="str":
            self.seq= full_seq
        elif return_as=="seq":
            seq=Seq(full_seq)
            self.seq=seq
        else:
            raise NotImplementedError("{} instance is not implemented the current options are 'str' for string and 'seq' for Bio.Seq")

        return self

    #TODO a __call__ method would be helpful here


class VariantList:
    def __init__(self, list):
        """
        this is a list of Variant instances
        :param list: a list of variant instances if any of the items in the list is not of Variant class raises ValueError
        """
        check=all([item==Variant for item in list])
        if not check:
            raise ValueError("not all items in the list are of class Variant")
        else:
            self.variants=list

    def append(self, variant):
        self.variants.append(variant)

    def write(self, path, type='tsv', **kwargs):
        """
        write results to file
        :param path: path of the file
        :param type: either tsv or fasta, if tsv columns will be chrom, pos, ref, alt, seq if fasta record names will be
        "chrom_pos_ref_alt"
        :param kwargs: these are passed to pd.to_csv
        :return: Nothing but writes to file
        """
        if os.path.exists(path):
            raise FileExistsError("{} already exists".format(path))

        #TODO this is a slooooow implementation, I would need to re-write this in numpy or even polars or cudf to speed
        # it up to millions
        if type=="tsv":
            df={"chrom":[],
                "pos":[],
                "ref":[],
                'alt':[],
                'seq':[]}

            for item in self.variants:
                df["chrom"].append(item.chrom)
                df["pos"].append(item.chrom)
                df["ref"].append(item.chrom)
                df["alt"].append(item.chrom)
                df["seq"].append(item.chrom)

            df=pd.DataFrame.from_records(df)
            df.to_csv(path, **kwargs)

        elif type=="fasta":
            records=[]
            for item in self.variants:
                name="_".join([item.chrom, item.pot, item.ref, item.alt])
                if type(item.seq)==str:
                    seq=Seq(item)
                else:
                    seq=item.seq
                records.append(SeqRecord(seq=seq, id=name, name="", description=""))
            write(records, path, "fasta")

        else:
            raise NotImplementedError("{} is not implemented the current options are tsv and fasta")



# Example Test
# 1. Create an Instance of Genome
genome_instance = Genome("/hpf/largeprojects/ccmbio/vikas_kaushik/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa")
# 2. Fetch a Region
sequence = genome_instance.fetch("3", 1000, 2000)
print("Fetched Sequence:", sequence)
# 3. Create an Instance of Variant
variant_instance = Variant(chrom="3", pos=1500, ref="A", alt="G", genome=genome_instance)
# 4. Mutate the Variant
mutated_sequence = variant_instance.mutate(total=100, ratio=0.3, return_as='str')
print("Mutated Sequence:", mutated_sequence.seq)
# 5. Create an Instance of VariantList
variant_list_instance = VariantList([variant_instance])
# 6. Append to VariantList
new_variant_instance = Variant(chrom="3", pos=2000, ref="C", alt="T", genome=genome_instance)
variant_list_instance.append(new_variant_instance)
# 7. Write to File
output_file_path = "path/to/output.tsv"
variant_list_instance.write(output_file_path, type='tsv')
print(f"Variants written to {output_file_path}")
