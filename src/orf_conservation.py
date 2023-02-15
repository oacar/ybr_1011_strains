import os
import re
import sys
import pandas as pd
import subprocess
from Bio import SeqIO



locus = SeqIO.parse("data/interim/1011_ybr_rev_cons.fasta",'fasta')
locus_list= [i for i in locus]
locus_list.sort(key=lambda x:x.id)
ybr = SeqIO.parse("data/interim//1011_ybr_orf_rev_cons.fasta",'fasta')
ybr_list= [i for i in ybr]
ybr_list.sort(key=lambda x:x.id)

ybr_translation = [i.translate(id=True,name='',description='') for i in ybr_list]
stop_pos =[]
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
df.to_csv("data/interim/strain_ybr_start_stop.csv")


