import os
import sys
import pandas as pd
import subprocess
from Bio import SeqIO

df = pd.read_csv('ybr_locus_intersect_bcftools.gvcf','\t',skiprows=61)
df = df.rename(columns={"ALT.1":"ALT"})

subprocess.run("bcftools view -R ybr_locus.bed 1011Matrix.gvcf.gz > ybr_locus_intersect_bcftools.gvcf",shell=True)

subprocess.run("gatk CreateSequenceDictionary -R Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",shell=True)

subprocess.run("gatk IndexFeatureFile -I ybr_locus_intersect_bcftools.gvcf",shell=True)
for i in df.columns[9:]:   
    cmd = f"gatk SelectVariants -R Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V ybr_locus_intersect_bcftools.gvcf --sample-name {i} -O sample_vcf/{i}.vcf"
    subprocess.run(cmd,shell=True)

for i in df.columns[9:]:   
    cmd = f"gatk FastaAlternateReferenceMaker -R Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -V sample_vcf/{i}.vcf  -O sample_genome/{i}.fa"
    subprocess.run(cmd,shell=True)

for i in df.columns[9:]:   
    cmd = f"bedtools getfasta -fi sample_genome/{i}.fa -bed ybr_locus.bed -fo sample_genome_ybr/{i}.fa"
    subprocess.run(cmd,shell=True)

files = os.listdir('sample_genome_ybr/')
records = [] 
for f in files:
    id_ = f.split('.')[0]
    for record in SeqIO.parse(f"sample_genome_ybr/{f}",'fasta'):
        record.id=id_
        records.append(record.reverse_complement())

SeqIO.write(records,"1011_ybr_rev.fasta","fasta")

for i in df.columns[9:]:   
    cmd = f"samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa chromosome2:612230-615856| bcftools consensus --sample {i} 1011Matrix.gvcf.gz > sample_genome_ybr_consensus/{i}.fasta"
    subprocess.run(cmd,shell=True)

files = os.listdir('sample_genome_ybr_consensus/')
records = [] 
for f in files:
    id_ = f.split('.')[0]
    for record in SeqIO.parse(f"sample_genome_ybr_consensus/{f}",'fasta'):
        record = record.reverse_complement()
        record.id=id_
        records.append(record)

SeqIO.write(records,"1011_ybr_rev_cons.fasta","fasta")


