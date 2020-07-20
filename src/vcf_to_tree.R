library(SNPRelate)
library(ape)
library(tidyverse)
setwd('src')

# Create phylogenetic tree from gvcf file -------
##DO NOT RUN
# vcf_file <- '../data/raw/1011Matrix.gvcf.gz'
# 
# snpgdsVCF2GDS(vcf_file,'../data/interim/biallelic.gds', method='biallelic.only')
# 
# snpgdsSummary("../data/interim/biallelic.gds")
# 
# genofile <- snpgdsOpen("../data/interim/biallelic.gds")
# 
# diss <- snpgdsDiss(genofile, autosome.only=FALSE)
# diss.mat <- diss$diss
# 
# phy <- bionj(diss.mat)
load('../data/interim/bionj_image.Rdata')
phy$tip.label <- diss$sample.id
png('../figures/phy_bionj.png')
plot(phy)
dev.off()
#save.image(file='../data/interim/bionj_image.Rdata')

#END DO NO RUN

#library(tidyverse)
#tree = 
# Get binary matrix - ------
data <- read_csv("../data/interim/ybr_promoter_motif_df.csv") # %>%filter(motif_center>500)
strain_df <- read_csv("../data/interim/strain_ybr_start_stop.csv") %>%
  mutate(start_codon = start) %>%
  bind_rows(data.frame(strain_name = "Scer", start_codon = "M", stop_pos = 49))
data <- data %>%
  full_join(strain_df, by = c("sequence_name" = "strain_name")) %>%
  inner_join(data %>% group_by(motif_alt_id) %>% tally()) %>%
  filter(n > 12)
data_orf_cons <- data %>% 
  mutate(orf_conserved = ifelse(start_codon == "M" & (stop_pos == 49 | stop_pos == -1), "ORF retained", "ORF not retained"))

color_df <- data_orf_cons %>%
  filter(sequence_name!='Scer')%>% 
  group_by(sequence_name, orf_conserved) %>%
  distinct(sequence_name,orf_conserved,.keep_all = T)%>%
  ungroup() %>% 
  #filter(motif_alt_id %in% (data %>% filter(sequence_name == "Scer") %>% select(motif_alt_id) %>% pull())) %>%
  select(orf_conserved,sequence_name) %>% mutate(color = ifelse(orf_conserved=='ORF retained','blue','orange'))


binary_mat <- data %>%
  filter(sequence_name!='Scer')%>% 
  group_by(sequence_name, motif_alt_id) %>%
  distinct(sequence_name,motif_alt_id,.keep_all = T)%>%
  ungroup() %>% 
  #filter(motif_alt_id %in% (data %>% filter(sequence_name == "Scer") %>% select(motif_alt_id) %>% pull())) %>%
  select(motif_alt_id,sequence_name) %>% 
  mutate(exists=1) %>%
  pivot_wider(names_from=motif_alt_id,values_from=exists,values_fill = 0) %>% 
  #mutate(sequence_name=factor(sequence_name, levels = (c("Scer", "Spar", "Smik", "Skud", "Sarb", "Suva", "Seub")))) %>% 
  #arrange(sequence_name) %>% 
  column_to_rownames(var="sequence_name") %>% 
  as.matrix() 


# Use strain data from paper -------
strain_data <- readxl::read_xls('../data/raw/strain_data.xls',sheet = 1,skip = 3)

data_with_clades <- data_orf_cons %>% inner_join(strain_data, by=c('sequence_name'='Standardized name'))
orf_not_retained_clades <- data_with_clades %>% filter(orf_conserved=='ORF not retained') %>% distinct(sequence_name, Clades, start_codon,stop_pos) #%>% view()
orf_not_retained_msn1 <- data_with_clades %>% filter(orf_conserved=='ORF not retained' & motif_alt_id=='Msn1p') %>% distinct(sequence_name, Clades, start_codon,stop_pos)
orf_not_retained_gcn4 <- data_with_clades %>% filter(orf_conserved=='ORF not retained' & motif_alt_id=='Gcn4p') %>% distinct(sequence_name, Clades, start_codon,stop_pos)
orf_not_retained_rgt1 <- data_with_clades %>% filter(orf_conserved=='ORF retained' & motif_alt_id=='Rgt1p') %>% distinct(sequence_name, Clades, start_codon,stop_pos)
orf_not_retained_rdr1 <- data_with_clades %>% filter(orf_conserved=='ORF not retained' & motif_alt_id=='Rdr1p') %>% distinct(sequence_name, Clades, start_codon,stop_pos)

# Plot complex heatmap --------
library(DECIPHER)
dend <- ReadDendrogram('../data/interim/bionj.nwk',convertBlanks = FALSE)
#gplots::heatmap.2(binary_mat_rwnames,trace='none',Rowv = FALSE,Colv = FALSE,dendrogram = 'none')

#coords <- locator(1)
# Remove rownames from matrix
binary_mat_rwnames <- binary_mat[labels(dend),]
rownames(binary_mat_rwnames) <- NULL#ifelse(rownames(binary_mat_rwnames)%in%orf_not_retained_clades$sequence_name,rownames(binary_mat_rwnames),'')
ha = ComplexHeatmap::HeatmapAnnotation(YBR = ifelse(labels(dend)%in%orf_not_retained_clades$sequence_name,'ORF broken','ORF retained'),
                                       which='row',col=list(YBR=c('ORF retained'='orange','ORF broken'='blue')))

pdf("../figures/heatmap_complex.pdf", width = 8, height = 8)
hm <- ComplexHeatmap::Heatmap(binary_mat_rwnames[,c('Msn1p','Gcn4p','Rgt1p','Rdr1p')],col = c('white','red'), 
                              cluster_rows = dend %>% set("by_labels_branches_col",    value = orf_not_retained_clades$sequence_name,
                                                          TF_values = c('blue','orange')),
                                #dendextend::set("branches_col", ifelse(rownames(binary_mat_rwnames)%in%c('BND','BEE','BMR'),'blue','orange')),
                           row_dend_width = unit(8, "cm"),
                           row_dend_reorder = FALSE,
                           right_annotation = ha,
                           column_title = "TFBS", 
                           row_title = "S.cerevisiae Strains",
                           #col = colors,
                           border = TRUE,
                           heatmap_legend_param = list(
                             at = c(0, 1),
                             labels = c("False", "True"),
                             title = "TFBS Found",
                             legend_height = unit(4, "cm")
                             #title_position = "lefttop-rot"
                           ),
                           use_raster = TRUE, 
                           raster_device = "png")

hm
dev.off()

# Plot phylogenetic tree --------
library(phylocanvas)
phycanv_stop <- phylocanvas(phy,treetype = "rectangular")
for( i in unique(data$sequence_name[data$start_codon!='M'])){
  phycanv_stop <- phycanv_stop %>% style_node(i,color='blue',fillcolor = 'blue',strokecolor = 'blue')
}

#phycanv_start

#phycanv_stop <- phylocanvas(phy,treetype = "rectangular")
for( i in unique(data$sequence_name[data$stop_pos!=49])){
  phycanv_stop <- phycanv_stop %>% style_node(i,color='green',fillcolor = 'blue')
}

for( i in unique(data$sequence_name[data$stop_pos!=49 & data$start_codon!='M'])){
  phycanv_stop <- phycanv_stop %>% style_node(i,color='purple',fillcolor = 'blue')
}

phycanv_stop


