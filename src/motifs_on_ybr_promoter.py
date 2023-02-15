import os
import re
import sys
import pandas as pd
import subprocess
from Bio import SeqIO

ybr_block = [i for i in SeqIO.parse('data/interim/ybr_locus_scer_positive_strand.fasta','fasta')][0].reverse_complement()
# S288C_YBR196C-A_YBR196C-A_genomic.fsa is downloaded from https://www.yeastgenome.org/locus/ybr196c-a
ybr_orf =[i for i in SeqIO.parse('data/raw/S288C_YBR196C-A_YBR196C-A_genomic.fsa','fasta')][0]
# the start position of the YBR ORF in the YBR locus sequence
ybr_start = [m.start() for m in re.finditer(str(ybr_orf.seq).upper(),str(ybr_block.seq))][0]
# the end position of the YBR ORF in the YBR locus sequence
ybr_end = ybr_start+len(ybr_orf)

motif_df = pd.read_csv('data/interim/ybr_locus_fimo/fimo.tsv','\t',skipfooter=3)
motif_df['motif_center'] = pd.to_numeric(round((motif_df.start+motif_df.stop)/2),downcast='integer')
ybr_promoter_motifs_500bp= motif_df.loc[(motif_df.start<(ybr_start))&(motif_df.start>(ybr_start-500))]
ybr_promoter_motifs_500bp.to_csv('data/interim/ybr_promoter_motif_df.csv')
