library(VariantAnnotation, quietly=T)
library(data.table, quietly=T)
library(ggplot2, quietly=T)
library(plyr, quietly=T)
library(reshape2, quietly=T)
library(seqinr, quietly=T)
library(tools, quietly=T)
library(zoo, quietly=T)

source("utils.R")

args <- c('1011Matrix.gvcf', 'output/','Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa', 10, 100,10)#(commandArgs(TRUE))
vcf_file <- args[1] # the vcf file
output_dir <- args[2] # full path for output directory
reference_fasta <- args[3] # the reference genome
density_window <- as.numeric(args[4]) # in base pair
sliding_window <- as.numeric(args[5]) # number of density windows
snp_allowed_per_window <- as.numeric(args[6]) # in base pair

print(vcf_file)
print(output_dir)
print(reference_fasta)
print(density_window)
print(sliding_window)
print(snp_allowed_per_window)

# GENOME CREATION
genome <- read.fasta(file = reference_fasta, as.string = TRUE)
organism_name <- basename(file_path_sans_ext(reference_fasta))
chr_length <- getLength(genome)
chr_names <- getName(genome)
chr_count <- length(chr_length)
cum_chr_length <- c(0, cumsum(chr_length[1:chr_count-1]))  # Cumulative sum, required to find "a variant's true position" in the genome (because it starts at position 1 each chromosome in the vcf)
names(chr_length) <- chr_names
names(cum_chr_length) <- chr_names

cols = c("cyan3", "darkorange", "red", "forestgreen", "slateblue3", "honeydew4", "brown", "blue", "darkgoldenrod1", "darkmagenta", "firebrick", "darkorchid2", "hotpink4", "mediumaquamarine", "rosybrown3", "royalblue", "red")

if (chr_count > 16){
  exit("WARNING: Colors are only specified for genome with less than 16 chromosomes.
       Please add colors to the 'cols' variable if you have more than 16 chromosomes.")
} else {
  print(paste("Organism:", organism_name, ". Nr of chromosomes:", chr_count, sep=" "))
}

cols = cols[1:chr_count]

# READ GATK FILE
readGatkVcf <- function(vcf_file, organism_name){
  param <- ScanVcfParam(info=c("ABHet", "ABHom"), format="DP")
  vcf <- readVcf(vcf_file, organism_name, param)
  data <- data.frame(position = rowData(vcf)@ranges@start, coverage = as.vector(geno(vcf)$DP), chromosome = seqnames(vcf), ABHet = info(vcf)$ABHet, ABHom = info(vcf)$ABHom)
  data$chromosome <- factor(data$chromosome, levels=chr_names) 
  data.m <- melt(data, id.vars=c("position", "chromosome", "coverage"), na.rm=TRUE) 
  data <- data.table(data.m)
  data[, true_position := position + cum_chr_length[as.character(chromosome)]]
  data
}

# CREATE A DIRECTORY PER STRAIN
create_strain_dir <- function(ouput_dir, strain_id){
  dir.create(file.path(output_dir, strain_id), showWarnings = FALSE)
}

# LOAD SNPS AND CREATE OUTPUT DIRECTORY
snps = readGatkVcf(vcf_file, organism_name)
snps_heterozygous = snps[variable %in% c("ABHet")]
snps_homozygous = snps[variable %in% c("ABHom")]

strain_id = unlist(strsplit(basename(vcf_file), "\\."))[1]
print(paste("Creating plot for ", strain_id, "...", sep=""))
create_strain_dir(output_dir, strain_id)

# Get heterozygous SNP density per window of n size
density.heterozygous = data.table(ldply(chr_names, function(chromo){
  snps_chromo = snps_heterozygous[chromosome %in% chromo]
  value = CUT(snps_chromo$position, seq(0, chr_length[as.character(chromo)], density_window), include.lowest = TRUE, right = FALSE, dig.lab=10)
  output.data = data.frame(chromosome = chromo, density = table(value$output), value$ranges)
}))

write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".HET.SNP.density.", density_window/1000, "kb_window.txt"), row.names = FALSE, x = density.heterozygous, quote = FALSE, sep = "\t")

# Get homozygous SNP density per window of n size
density.homozygous = data.table(ldply(chr_names, function(chromo){
  snps_chromo = snps_homozygous[chromosome %in% chromo]
  value = CUT(snps_chromo$position, seq(0, chr_length[as.character(chromo)], density_window), include.lowest = TRUE, right = FALSE, dig.lab=10)
  output.data = data.frame(chromosome = chromo, density = table(value$output), value$ranges)
}))

write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".HOM.SNP.density.", density_window/1000, "kb_window.txt"), row.names = FALSE, x = density.homozygous, quote = FALSE, sep= "\t")


all.density <- data.frame(density.heterozygous, "density.homozygous" = density.homozygous$density.Freq)
all.density.m <- melt(all.density)

write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".SNP.density.", density_window/1000, "kb_window.txt"), row.names = FALSE, x = all.density, quote = FALSE, sep= "\t")

# SNP heterozygous
density_sliding_window = data.table(ldply(chr_names, function(chromo){
  density_chromo = density.heterozygous[chromosome == chromo]
  sum_per_sliding_window = rollsum(density_chromo$density.Freq, sliding_window)
  lower = head(as.numeric(as.vector(density_chromo$lower)), -sliding_window+1) # smaller vector 'cause of sliding window 
  upper = head(as.numeric(as.vector(density_chromo$upper))+sliding_window*density_window, -sliding_window+1)
  output.data = data.frame(chromosome = chromo, density.sum = sum_per_sliding_window, start = lower, end = upper)
}))

# SNP homozygous
homo_density_sliding_window = data.table(ldply(chr_names, function(chromo){
  density_chromo = density.homozygous[chromosome == chromo]
  sum_per_sliding_window = rollsum(density_chromo$density.Freq, sliding_window)
  lower = head(as.numeric(as.vector(density_chromo$lower)), -sliding_window+1) # smaller vector 'cause of sliding window 
  upper = head(as.numeric(as.vector(density_chromo$upper))+sliding_window*density_window, -sliding_window+1)
  output.data = data.frame(chromosome = chromo, density.sum = sum_per_sliding_window, start = lower, end = upper)
}))

write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".LOH.density.", snp_allowed_per_window , "snp_allowed_per_", sliding_window, "kb_sliding_window.txt"), x = density_sliding_window, row.names = FALSE, quote = FALSE, sep="\t")
write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".LOH.density.", snp_allowed_per_window , "snp_allowed_per_", sliding_window, "kb_sliding_window.txt"), x = homo_density_sliding_window, row.names = FALSE, quote = FALSE, sep="\t")

# Print LOH regions
loh_regions = data.table(ldply(chr_names, function(chromo){  
  if(nrow(density_sliding_window[density.sum <= snp_allowed_per_window & chromosome == chromo]) > 0){
    ir <- IRanges(as.numeric(density_sliding_window[density.sum <= snp_allowed_per_window & chromosome == chromo]$start), as.numeric(as.character(density_sliding_window[density.sum <= snp_allowed_per_window & chromosome == chromo]$end)))
    data.frame("chromosome" = chromo, reduce(ir))
  }
}))

write.table(file = paste0(output_dir, "/", strain_id, "/", strain_id, ".LOH.regions.", snp_allowed_per_window , "snp_allowed_per_", sliding_window, "kb_sliding_window.txt"), x = loh_regions, row.names = FALSE, quote = FALSE, sep="\t")

###################

p = ggplot(snps_heterozygous) + geom_histogram(aes(position), binwidth=density_window) + facet_wrap(~ chromosome) + geom_rect(data=loh_regions, aes(xmin=start, xmax=end, ymin=0, ymax=+Inf), fill="red")

q = ggplot(snps_heterozygous) + geom_point(aes(position, value, colour = chromosome)) + scale_colour_manual(values=cols) + ylab("Allele balance") + ylim(0, 1) + facet_wrap(~ chromosome)
q = q + geom_rect(data=loh_regions, aes(xmin=start, xmax=end, ymin=0, ymax=+Inf), fill="red", alpha = 0.5)

ggsave(plot = p, paste0(output_dir, "/", strain_id, "/", strain_id, ".LOH.frequence.", snp_allowed_per_window , "snp_allowed_per_", sliding_window, "kb_sliding_window.pdf"), width=20, height=12)
ggsave(plot = q, paste0(output_dir, "/", strain_id, "/", strain_id, ".LOH.allele_balance.", snp_allowed_per_window , "snps_allowed_per_", sliding_window, "kb_sliding_window.pdf"), width=20, height=12)


 
 
 
 
 
 
 
 
