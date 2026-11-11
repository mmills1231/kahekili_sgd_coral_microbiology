# Kahekili 16S analysis: taxonomic composition
# Author: Makeda I. Mills

source(file.path("scripts", "01_16S_analysis", "02_16S_setup.R"))

# Relative abundance tables ####

class_abun <- mauiseq16S %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  arrange(Class)

genus_abun <- mauiseq16S %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  arrange(Genus)

species_abun <- mauiseq16S %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  arrange(Species)

write.csv(
  class_abun,
  file.path(results_dir, "kahekili_16S_class_abundance.csv"),
  row.names = FALSE
)
write.csv(
  genus_abun,
  file.path(results_dir, "kahekili_16S_genus_abundance.csv"),
  row.names = FALSE
)
write.csv(
  species_abun,
  file.path(results_dir, "kahekili_16S_species_abundance.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Retain genera represented at >1% in at least one sample, then identify the
# 20 genera with the highest mean abundance within that filtered set.
# -----------------------------------------------------------------------------

genus_filtered <- genus_abun %>%
  select(Genus, Site, Sample_type, Sample, Abundance) %>%
  group_by(Genus, Site, Sample_type, Sample) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  filter(avg_abundance > 0.01)

top20_genera <- genus_filtered %>%
  group_by(Genus) %>%
  summarize(total_abundance = mean(avg_abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# -----------------------------------------------------------------------------
# Collapse genera outside the retained set into "Others".
# -----------------------------------------------------------------------------

genus_plot_data <- genus_abun %>%
  mutate(
    Genus = as.character(Genus),
    Genus = if_else(Genus %in% top20_genera, Genus, "Others")
  ) %>%
  group_by(Sample, Site, Sample_type, Genus) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# Preserve the date-based sample ordering used in the working analysis when a
# Date column is available in the metadata.
if ("Date" %in% colnames(genus_abun)) {
  sample_dates <- genus_abun %>%
    select(Sample, Date) %>%
    distinct() %>%
    mutate(Date_only = as.Date(as.character(Date), format = "%m/%d/%y"))

  genus_plot_data <- genus_plot_data %>%
    left_join(sample_dates, by = "Sample") %>%
    arrange(Date_only, Sample) %>%
    mutate(Sample_ordered = factor(Sample, levels = unique(Sample)))
} else {
  genus_plot_data <- genus_plot_data %>%
    mutate(Sample_ordered = factor(Sample, levels = unique(Sample)))
}

write.csv(
  genus_plot_data,
  file.path(results_dir, "kahekili_16S_genus_abundance_collapsed.csv"),
  row.names = FALSE
)

# Taxonomic composition barplot ####

composition_plot <- ggplot(genus_plot_data) +
  geom_col(
    aes(x = Sample_ordered, y = Abundance, fill = Genus),
    position = "fill",
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = coral_colors3_20) +
  facet_grid(Sample_type ~ Site, scales = "free", space = "free_x") +
  ylab("Proportion of Community") +
  xlab("Samples") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(
    axis.text.y.left = element_text(size = 8, family = "Helvetica", color = "black"),
    axis.text.x = element_text(size = 7, hjust = 1, angle = 90,
                               family = "Helvetica", color = "black"),
    axis.title.y = element_text(size = 9, family = "Helvetica", color = "black"),
    strip.background = element_rect(fill = "transparent", color = "black", linewidth = 0.75),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    panel.border = element_rect(color = "black", linewidth = 0.75)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

composition_plot

ggsave(
  file.path(results_dir, "kahekili_16S_taxonomic_composition.pdf"),
  plot = composition_plot,
  width = 20,
  height = 8
)
