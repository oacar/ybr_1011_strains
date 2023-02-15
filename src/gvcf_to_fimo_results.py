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

# ybr_locus_intersect_bcftools.gvcf is only used for the column names of the dataframe
# ybr_locus_scer_positive_strand.fasta is the reference sequence of the YBR locus
# ybr_locus_consensus_per_genome is a folder with the consensus sequences of the YBR locus for each strain
# 1011_ybr_rev_cons.fasta is the consensus sequences of the YBR locus on the negative strand for each strain concatenated



subprocess.run("bcftools view -R data/raw/ybr_locus.bed data/raw/1011Matrix.gvcf.gz > data/interim/ybr_locus_intersect_bcftools.gvcf",shell=True)

subprocess.run("bedtools getfasta -fi data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -bed data/raw/ybr_locus.bed -fo data/interim/ybr_locus_scer_positive_strand.fasta",shell=True)
df = pd.read_csv('data/interim/ybr_locus_intersect_bcftools.gvcf',sep='\t',skiprows=61)
df = df.rename(columns={"ALT.1":"ALT"})
subprocess.run("gatk CreateSequenceDictionary -R data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",shell=True)

subprocess.run("gatk IndexFeatureFile -I data/interim/ybr_locus_intersect_bcftools.gvcf",shell=True)


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
# add the reference sequence
for record in SeqIO.parse('data/interim/ybr_locus_scer_positive_strand.fasta','fasta'):
    record = record.reverse_complement()
    record.id='scer'
    record.description=''
    records.append(record)
SeqIO.write(records,"data/interim/1011_ybr_rev_cons.fasta","fasta")

# run FIMO on the locus consensus sequences to get the motif positions
# First calculate background frequencies
subprocess.run("fasta-get-markov data/interim/1011_ybr_rev_cons.fasta > data/interim/1011_ybr_rev_cons_fimo_background",shell=True)

# Then run FIMO
subprocess.run("fimo -norc -bfile data/interim/1011_ybr_rev_cons_fimo_background -oc data/interim/ybr_locus_fimo /home/oma21/.local/share/meme-5.1.1/db/motif_databases/YEAST/YEASTRACT_20130918.meme data/interim/1011_ybr_rev_cons.fasta", shell=True)



# get YBR ORF sequences
# orf coordinates are from https://www.yeastgenome.org/locus/YBR196c-a
# 614024-614173


os.mkdir('data/interim/ybr_orf_consensus_per_genome')
for i in df.columns[9:]:   
    cmd = f"samtools faidx data/raw/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa chromosome2:614024-614173| bcftools consensus --sample {i} data/raw/1011Matrix.gvcf.gz > data/interim/ybr_orf_consensus_per_genome/{i}.fasta"
    subprocess.run(cmd,shell=True)

files = os.listdir('data/interim/ybr_orf_consensus_per_genome/')
records = [] 
for f in files:
    id_ = f.split('.')[0]
    # if id_ == 'ALT':
    #     id_ = 'ALT.1'
    for record in SeqIO.parse(f"data/interim/ybr_orf_consensus_per_genome/{f}",'fasta'):
        record = record.reverse_complement()
        record.id=id_
        record.description=''
        records.append(record)

SeqIO.write(records,"data/interim/1011_ybr_orf_rev_cons.fasta","fasta")

# translate YBR ORF sequences
from Bio.Seq import Seq
from Bio import SeqIO
records = []
for record in SeqIO.parse('data/interim/1011_ybr_orf_rev_cons.fasta','fasta'):
    # if len(record.seq) % 3 != 0:
    #     print(record.id, len(record.seq))
    #     print(record.seq.translate())
    record.seq = record.seq.translate()
    records.append(record)

SeqIO.write(records,"data/interim/1011_ybr_orf_rev_cons_aa.fasta","fasta")



