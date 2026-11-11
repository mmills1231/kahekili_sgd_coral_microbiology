# Kahekili 16S analysis: alpha diversity
# Author: Makeda I. Mills

source(file.path("scripts", "01_16S_analysis", "02_16S_setup.R"))

# Shannon diversity ####

alpha_df <- estimate_richness(rare_mauiseq, measures = "Shannon") %>%
  rownames_to_column("SampleID")

enviro_data <- data.frame(sample_data(rare_mauiseq)) %>%
  rownames_to_column("SampleID")

alpha_env <- left_join(alpha_df, enviro_data, by = "SampleID")

write.csv(
  alpha_env,
  file.path(results_dir, "kahekili_16S_shannon_diversity.csv"),
  row.names = FALSE
)

# Shannon diversity by sample type ####

alpha_plot <- plot_richness(
  rare_mauiseq,
  x = "Site",
  measures = "Shannon",
  color = "Sample_type"
) +
  geom_boxplot(
    aes(fill = Sample_type, color = Sample_type),
    alpha = 0.90,
    show.legend = FALSE
  ) +
  geom_point(aes(fill = Sample_type), alpha = 1) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12,
      family = "Helvetica",
      color = "black"
    ),
    panel.background = element_rect(fill = "transparent", color = "black"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = sampletype_colors) +
  scale_fill_manual(values = sampletype_colors)

alpha_plot

ggsave(
  file.path(results_dir, "kahekili_16S_shannon_by_site_sample_type.pdf"),
  plot = alpha_plot,
  width = 8,
  height = 6
)

# Statistical tests from the working analysis ####

shapiro_result <- shapiro.test(alpha_env$Shannon)
kruskal_sample_type <- kruskal.test(Shannon ~ Sample_type, data = alpha_env)
pairwise_sample_type <- pairwise.wilcox.test(
  alpha_env$Shannon,
  alpha_env$Sample_type,
  p.adjust.method = "BH"
)

shapiro_result
kruskal_sample_type
pairwise_sample_type

# The working script also explored a Site + Sample_type linear model. It is
# retained here as an exploratory output rather than substituted for the
# non-parametric tests above.
lm_model <- lm(Shannon ~ Site + Sample_type, data = alpha_env)
anova(lm_model)

