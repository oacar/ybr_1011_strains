Need this raw files in `data/raw/` folder
- 1011Matrix.gvcf.gz
- S288C_YBR196C-A_YBR196C-A_genomic.fsa
- Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
- Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.fai
- ybr_locus.bed 
    content is : `chromosome2	612230	615856`

YBR locus is supposed to be the region from PGI1 to YBR197C, including these 2.

Run following in order:
1. `python3 gvcf_to_fimo_results.py`
2. `python3 orf_conservation.py`
3. `python3 motifs_on_ybr_promoter.py`
4. `Rscript orf_conservation_motif_ratios.R`

This will generate the following files in `data/interim/` folder

- strain_ybr_start_stop.csv (this is the file that contains the start and stop of the ybr ORF in each strain)
- ybr_promoter_motif_df.csv (this is the file that contains the motifs in the ybr promoter)
- ybr_locus_fimo/ (folder) (this is the folder that contains the fimo results for each strain)
- ybr_orf_consensus_per_genome/ (folder) (this is the folder that contains the consensus sequence for each strain)
- 1011_ybr_orf_rev_cons.fasta (this is the consensus sequence for the ybr ORF in the 1011 strain)
- 1011_ybr_rev_cons.fasta (this is the consensus sequence for the ybr locus in the 1011 strain)