import os
import sys
import pandas as pd
import subprocess
from Bio import SeqIO

#ybr_locus.bed is a bed file with the coordinates of the YBR locus
# that spans from 612230 to 615856
# chromosome format is chromosome2
#1011Matrix.gvcf.gz is a gvcf file with the genotypes of 1011 yeast strains
#Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa is the reference genome of yeast

subprocess.run("bcftools view -R data/interim/ybr_locus.bed data/raw/1011Matrix.gvcf.gz > data/interim/ybr_locus_intersect_bcftools.gvcf",shell=True)

df = pd.read_csv('data/interim/ybr_locus_intersect_bcftools.gvcf',sep='\t',skiprows=61)
df = df.rename(columns={"ALT.1":"ALT"})
subprocess.run("gatk CreateSequenceDictionary -R data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",shell=True)

subprocess.run("gatk IndexFeatureFile -I data/interim/ybr_locus_intersect_bcftools.gvcf",shell=True)

# for i in df.columns[9:]:   
#     cmd = f"gatk SelectVariants -R data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V data/raw/1011Matrix.gvcf.gz --restrict-alleles-to BIALLELIC -O data/interim/biallelic.vcf"
#     subprocess.run(cmd,shell=True)

# for i in df.columns[9:]:   
#     cmd = f"gatk FastaAlternateReferenceMaker -R data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V sample_vcf/{i}.vcf  -O sample_genome/{i}.fa"
#     subprocess.run(cmd,shell=True)

# for i in df.columns[9:]:   
#     cmd = f"bedtools getfasta -fi sample_genome/{i}.fa -bed ybr_locus.bed -fo sample_genome_ybr/{i}.fa"
#     subprocess.run(cmd,shell=True)

# files = os.listdir('sample_genome_ybr/')
# records = [] 
# for f in files:
#     id_ = f.split('.')[0]
#     for record in SeqIO.parse(f"sample_genome_ybr/{f}",'fasta'):
#         record.id=id_
#         records.append(record.reverse_complement())

# SeqIO.write(records,"1011_ybr_rev.fasta","fasta")
os.mkdir('data/interim/ybr_locus_consensus_per_genome')
for i in df.columns[9:]:   
    cmd = f"samtools faidx data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa chromosome2:612230-615856| bcftools consensus --sample {i} data/raw/1011Matrix.gvcf.gz > data/interim/ybr_locus_consensus_per_genome/{i}.fasta"
    subprocess.run(cmd,shell=True)

files = os.listdir('data/interim/ybr_locus_consensus_per_genome/')
records = [] 
for f in files:
    id_ = f.split('.')[0]
    # if id_ == 'ALT':
    #     id_ = 'ALT.1'
    for record in SeqIO.parse(f"data/interim/ybr_locus_consensus_per_genome/{f}",'fasta'):
        record = record.reverse_complement()
        record.id=id_
        record.description=''
        records.append(record)

SeqIO.write(records,"data/interim/1011_ybr_rev_cons.fasta","fasta")


