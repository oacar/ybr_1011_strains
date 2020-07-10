library(tidyverse)
# setwd('~/OneDrive - University of Pittsburgh/ybr/ybr_lyon/motif/')
data <- read_csv("../data/interim/ybr_promoter_motif_df.csv") # %>%filter(motif_center>500)
strain_df <- read_csv("../data/interim/strain_ybr_start_stop.csv") %>%
  mutate(start_codon = start) %>%
  bind_rows(data.frame(strain_name = "Scer", start_codon = "M", stop_pos = 49))
data <- data %>%
  full_join(strain_df, by = c("sequence_name" = "strain_name")) %>%
  inner_join(data %>% group_by(motif_alt_id) %>% tally()) %>%
  filter(n > 12)
# data$sequence_name <- factor(data$sequence_name)
# pdf("figures/motif_map_legend_r.pdf", height = 12, width = 12)
data %>% filter(motif_alt_id %in% c("Put3p")) %>%
  # distinct_at(vars(sequence_name, motif_center, motif_alt_id), .keep_all = T) %>%
  # filter(motif_alt_id %in% (data %>% filter(sequence_name == "Scer") %>% select(motif_alt_id) %>% pull())) %>%
  ggplot(aes(x = motif_center, y = sequence_name, color = motif_alt_id, fill = motif_alt_id)) +
  geom_point(position = ggstance::position_dodgev(height = 0.5)) + # scale_y_discrete(labels=ifelse(sequence_name=='Scer','Scer',''))+
  # theme_void() + #aes(shape = glucose_related),
  # theme(legend.position = "none") +

  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  # xlim(1683,1183)
  scale_x_reverse(limits = c(1700, 1183))
# dev.off()

data %>%
  filter(motif_alt_id %in%
    (data %>% filter(sequence_name == "Scer") %>% select(motif_alt_id) %>% pull()) == F) %>%
  select(motif_alt_id) %>%
  distinct() %>%
  pull()
# [1] "Gsm1p"   "Tod6p"   #"Msn1p"#   "Dot6p"   "Ino4p"   "Ino2p"   "Ste12p"  "Hac1p"   "Cup2p"   #"Pdr3p"#   "Rim101p" "Smp1p"   "Ime1p"   "Aft1p"   "Cbf1p"
# [16] "Rgt1p"   "Mbp1p"   "Ert1p"   "Hsf1p"   #"Pdr1p"#   "Spt15p"  "Put3p"
# data%>%filter(sequence_name%in%c('Scer','Sarb','Suva')&motif_alt_id%in%c('Mig1p'))%>%view
# data%>%filter(sequence_name%in%c('Scer','Spar','Smik')&motif_alt_id%in%c('Mig1p'))%>%view

# data%>%filter(sequence_name=='Suva')%>%view
data %>% filter(sequence_name == "Scer" & motif_alt_id == "Rdr1p") # %>%select(motif_alt_id)%>%distinct()%>%pull
data %>% filter(motif_alt_id == "Rdr1p")
# Motif finding ratio w.r.t start codons-------
start_codon_counts <- data %>%
  distinct(start_codon, sequence_name) %>%
  group_by(start_codon) %>%
  tally()

data %>%
  group_by(start_codon, motif_alt_id) %>%
  distinct(sequence_name, motif_alt_id, .keep_all = T) %>%
  summarise(count = n()) %>%
  inner_join(start_codon_counts) %>%
  mutate(ratio = count / n) %>%
  ggplot(aes(x = motif_alt_id, y = ratio, fill = start_codon, color = start_codon)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Motif finding ratio w.r.t stop codons ----------
stop_pos_counts <- data %>%
  mutate(stop_pos = ifelse(stop_pos == 49, "same", "different")) %>%
  distinct(stop_pos, sequence_name) %>%
  group_by(stop_pos) %>%
  tally()

data %>%
  mutate(stop_pos = ifelse(stop_pos == 49, "same", "different")) %>%
  # filter(stop_pos%in%c(49,24))%>%
  group_by(stop_pos, motif_alt_id) %>%
  distinct(sequence_name, motif_alt_id, .keep_all = T) %>%
  summarise(count = n()) %>%
  inner_join(stop_pos_counts) %>%
  ungroup() %>%
  mutate(ratio = count / n, stop_pos = factor(stop_pos)) %>%
  ggplot(aes(x = motif_alt_id, y = ratio, fill = stop_pos, color = stop_pos)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Motif finding ratio w.r.t orf_conservation-------
data_orf_cons <- data %>% mutate(orf_conserved = ifelse(start_codon == "M" & (stop_pos == 49 | stop_pos == -1), "ORF retained", "ORF not retained"))
orf_cons_counts <- data_orf_cons %>%
  distinct(orf_conserved, sequence_name) %>%
  group_by(orf_conserved) %>%
  tally()

name_order <- data %>%
  select(motif_alt_id, motif_center) %>%
  distinct() %>%
  group_by(motif_alt_id) %>%
  summarise(center_mean = min(motif_center)) %>%
  arrange(desc(center_mean)) %>%
  select(motif_alt_id) %>%
  pull()

data_orf_cons %>%
  group_by(orf_conserved, motif_alt_id) %>%
  distinct(sequence_name, motif_alt_id, .keep_all = T) %>%
  summarise(count = n()) %>%
  inner_join(orf_cons_counts) %>%
  mutate(ratio = count / n, motif_alt_id = factor(motif_alt_id,levels=name_order)) %>%
  # filter(ratio>0.2)%>%
  ggplot(aes(x = motif_alt_id, y = ratio, color = orf_conserved)) +
  geom_point(position = position_dodge(width = 0.8)) +
  theme_bw() +
  theme(
    legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank(), axis.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  xlab("TF name") +
  ylab("%of strains with TFBS identified") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("#E7B800", "#52854C"))

ggsave("figures/orf_conserved.png", dpi = 150, height = 6, width = 10)
data_wide <- data_orf_cons %>%
  select(motif_alt_id, orf_conserved, sequence_name) %>%
  mutate(motif_exists = TRUE) %>%
  distinct() %>%
  pivot_wider(names_from = sequence_name, values_from = motif_exists) # ,values_fill = FALSE)

data_wide[is.na(data_wide)] <- FALSE
data_wide %>%
  filter(motif_alt_id == "Gcn4p") %>%
  pivot_longer(3:ncol(.)) # %>%# select(%>%  table()#%>%filter(value==FALSE)
p_val_list <- list()
for (i in unique(data$motif_alt_id)) {
  tbl <- data_orf_cons %>%
    group_by(orf_conserved, motif_alt_id) %>%
    distinct(sequence_name, motif_alt_id, .keep_all = T) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    complete(orf_conserved, motif_alt_id, fill = list(count = 0)) %>%
    inner_join(orf_cons_counts) %>%
    filter(motif_alt_id == i) %>%
    mutate(gcount = abs(count - n)) %>%
    select(count, gcount) # %>%
  # fisher.test()
  if (all(dim(tbl) == 2)) {
    pval <- tbl %>%
      fisher.test() %>%
      magrittr::use_series(p.value)
    p_val_list[i] <- pval
  }
}

data_orf_cons %>%
  group_by(orf_conserved, motif_alt_id) %>%
  distinct(sequence_name, motif_alt_id, .keep_all = T) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  complete(orf_conserved, motif_alt_id, fill = list(count = 0))

datasetALL %>%
  group_by(YEAR, Region) %>%
  summarise(count_number = n()) %>%
  ungroup() %>%
  complete(Year, Region, fill = list(count_number = 1))


data_wide %>%
  # filter(motif_alt_id == "Gcn4p") %>%
  pivot_longer(3:ncol(.)) %>%
  with(., table(orf_conserved, value)) %>%
  fisher.test() # %>% magrittr::use_series(p.value)
data_orf_cons %>%
  filter(motif_alt_id == "Gcn4p") %>%
  with(., table(motif_alt_id, orf_conserved))
data_orf_cons %>% complete(motif_alt_id, nesting(item_id, item_name), fill = list(value1 = 0))
