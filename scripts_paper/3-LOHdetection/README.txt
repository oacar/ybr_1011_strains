Script allowing for LOH region detection


This script take several parameters as argument:

vcf	vcf file outputted from GATK Variant Annotator (must contain allele balance information - ABHom/ABHet) (vcf format)
outdir	Output directory
reference	Reference genome initially used to map the reads (fasta format)
density_window	Size of the windows in which SNP frequency will be computed (1000 bp used in the 1002 project)
sliding_window	Number of window to slide. sliding size = sliding_window * density_window (50 used in the 1002 project)
snp_allowed_per_window	Number of SNPs allowed per (sliding window * density_window) (10 used in the 1002 project) 	


To run the script:
Rscript LOH_detection.R $vcf $outdir $reference $density_window $sliding_window $snp_allowed_per_window

