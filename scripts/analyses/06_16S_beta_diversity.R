# Kahekili 16S analysis: beta diversity
# Author: Makeda I. Mills

source(file.path("scripts", "01_16S_analysis", "02_16S_setup.R"))

set.seed(123)

# PCoA using Bray-Curtis dissimilarity ####

all_pcoa <- ordinate(
  physeq = rare_mauiseq,
  method = "PCoA",
  distance = "bray"
)

pcoa_plot <- plot_ordination(rare_mauiseq, ordination = all_pcoa) +
  geom_point(aes(fill = Site, shape = Sample_type), size = 5) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_fill_manual(values = a_coral_colors) +
  theme_classic() +
  theme(legend.background = element_rect(fill = "white", color = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  stat_ellipse(
    aes(color = Sample_type),
    level = 0.95,
    linetype = "dashed",
    linewidth = 1
  ) +
  theme(
    axis.text.y.left = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.text.x = element_text(size = 12, hjust = 1, family = "Helvetica", color = "black"),
    axis.title.y = element_text(size = 12, family = "Helvetica", color = "black")
  )

pcoa_plot

ggsave(
  file.path(results_dir, "kahekili_16S_pcoa_bray_curtis.pdf"),
  plot = pcoa_plot,
  width = 8,
  height = 6
)

# PERMANOVA and multivariate dispersion ####

# The working script intended to remove the well sample before these tests but
# compared sample names to the literal string "well". Here the intended Site
# filter is written explicitly.
rare_mauiseq_nowell <- subset_samples(rare_mauiseq, Site != "well")
rare_mauiseq_nowell <- prune_taxa(taxa_sums(rare_mauiseq_nowell) > 0, rare_mauiseq_nowell)

bray_curtis <- phyloseq::distance(rare_mauiseq_nowell, method = "bray")
sampledf <- data.frame(sample_data(rare_mauiseq_nowell))

write.csv(
  as.matrix(bray_curtis),
  file.path(results_dir, "kahekili_16S_bray_curtis_distance_matrix.csv")
)

permanova_sample_type <- adonis2(
  bray_curtis ~ Sample_type,
  data = sampledf,
  permutations = 999
)

permanova_site <- adonis2(
  bray_curtis ~ Site,
  data = sampledf,
  permutations = 999
)

permanova_sample_type
permanova_site

# Homogeneity of multivariate dispersion.
beta_sample_type <- betadisper(bray_curtis, sampledf$Sample_type)
beta_site <- betadisper(bray_curtis, sampledf$Site)

sample_type_dispersion <- permutest(beta_sample_type, permutations = 9999)
site_dispersion <- permutest(beta_site, permutations = 9999)

sample_type_dispersion
site_dispersion
anova(beta_site)
TukeyHSD(beta_site, conf.level = 0.95)

# Distance to centroid by site ####

disp_data <- data.frame(
  Distance = beta_site$distances,
  Site = beta_site$group
)

dispersion_plot <- ggplot(disp_data, aes(x = Site, y = Distance, fill = Site)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_manual(values = a_coral_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = "Site",
    y = "Distance to Centroid",
    title = "Microbial Variability of Samples by Site"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                               family = "Helvetica", color = "black"),
    axis.text.y = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.title.x = element_text(size = 12, family = "Helvetica", color = "black"),
    axis.title.y = element_text(size = 12, family = "Helvetica", color = "black"),
    plot.title = element_text(size = 12, family = "Helvetica", color = "black", hjust = 0.5)
  )

dispersion_plot

ggsave(
  file.path(results_dir, "kahekili_16S_beta_dispersion_site.pdf"),
  plot = dispersion_plot,
  width = 8,
  height = 6
)

