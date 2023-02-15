Need this raw files in `data/raw/` folder
- 1011Matrix.gvcf.gz
- S288C_YBR196C-A_YBR196C-A_genomic.fsa
- Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
- Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai
- ybr_locus.bed 
    content is : `chromosome2	612230	615856`

Run following in order:
1. `python3 gvcf_to_fimo_results.py`
2. `python3 orf_conservation.py`
3. `python3 motifs_on_ybr_promoter.py`

This will generate the following files in `data/interim/` folder

- strain_ybr_start_stop.csv
- ybr_promoter_motif_df.csv
- ybr_locus_fimo/ (folder)
- ybr_orf_consensus_per_genome/ (folder)
- 1011_ybr_orf_rev_cons.fasta
- 