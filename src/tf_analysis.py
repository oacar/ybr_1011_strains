import os
import re
import sys
import pandas as pd
import subprocess
from Bio import SeqIO

ybr_block = [i for i in SeqIO.parse('data/interim/YBR196C-A_scer_block.fa','fasta')][0]
ybr_orf =[i for i in SeqIO.parse('data/interim/YBR196C-A_scer_orf.fa','fasta')][0]

ybr_start = [m.start() for m in re.finditer(str(ybr_orf.seq).upper(),str(ybr_block.seq))][0]
ybr_end = ybr_start+len(ybr_orf)

locus = SeqIO.parse("data/interim/ybr_consensus/1011_ybr_rev_cons.fasta",'fasta')
locus_list= [i for i in locus]
locus_list.sort(key=lambda x:x.id)
ybr = SeqIO.parse("data/interim/ybr_consensus/1011_ybr_orf_rev_cons.fasta",'fasta')
ybr_list= [i for i in ybr]
ybr_list.sort(key=lambda x:x.id)

ybr_translation = [i.translate(id=True,name='',description='') for i in ybr_list]
non_atg = [i for i in ybr_translation if i[0]!='M']
len(non_atg)
len(ybr_list)
stop_pos =[]
SeqIO.write(ybr_translation,'data/interim/sequences/ybr_translation.fa','fasta')
for i in ybr_translation:
   s =  [m.start() for m in re.finditer('\*', str(i.seq))]
   if len(s)!=0:
       stop_pos.append(s[0])
   else:
       stop_pos.append(-1)

df = pd.DataFrame({"strain_name":[i.id for i in ybr_list]})
df['stop_pos'] = stop_pos
df['start']=[i.seq[0] for i in ybr_translation]
df['start']= df.start.astype('category')
df.groupby(['start','stop_pos']).count()
df.groupby('stop_pos').count()
df.to_csv("data/interim/strain_ybr_start_stop.csv")
motif_df = pd.read_csv('ybr_block/fimo.tsv','\t',skipfooter=3)
motif_df['motif_center'] = pd.to_numeric(round((motif_df.start+motif_df.stop)/2),downcast='integer')
motif_df.to_csv('data/interim/motif_df.csv')
scer_motifs = motif_df.loc[(motif_df.start<(ybr_start))&(motif_df.start>(ybr_start-500))]
scer_motifs.to_csv('data/interim/ybr_promoter_motif_df.csv')

df.loc[df.strain_name.isin(
scer_motifs.loc[(scer_motifs.motif_center>1600)& (scer_motifs.motif_alt_id=='Mig1p'),'sequence_name'].unique())]
