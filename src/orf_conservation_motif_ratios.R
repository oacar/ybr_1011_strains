library(tidyverse)
data <- read_csv("data/interim/ybr_promoter_motif_df.csv")
strain_df <- read_csv("data/interim/strain_ybr_start_stop.csv") %>%
    mutate(start_codon = start) %>%
    bind_rows(data.frame(strain_name = "scer", start = "M", start_codon = "M", stop_pos = 49))
data <- data %>%
    full_join(strain_df, by = c("sequence_name" = "strain_name")) %>%
    inner_join(data %>% group_by(motif_alt_id) %>% tally()) %>%
    filter(n > 12)


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
    mutate(ratio = count / n, motif_alt_id = factor(motif_alt_id, levels = name_order)) %>%
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

ggsave("figures/orf_conserved_021523.png", dpi = 150, height = 6, width = 10)
