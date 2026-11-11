# Kahekili 16S analysis: environmental relationships
# Author: Makeda I. Mills

source(file.path("scripts", "01_16S_analysis", "02_16S_setup.R"))
require(cowplot)
require(gt)

# Merge Shannon diversity with environmental metadata ####

alpha_df <- estimate_richness(rare_mauiseq, measures = "Shannon") %>%
  rownames_to_column("SampleID")

env_data_corr <- data.frame(sample_data(rare_mauiseq)) %>%
  rownames_to_column("SampleID")

alpha_env <- left_join(alpha_df, env_data_corr, by = "SampleID") %>%
  filter(Site != "well")

# Distribution check and sample-type/site comparisons retained from the working analysis.
shapiro.test(alpha_env$Shannon)
pairwise.wilcox.test(
  alpha_env$Shannon,
  alpha_env$Sample_type,
  p.adjust.method = "BH"
)
pairwise.wilcox.test(
  alpha_env$Shannon,
  alpha_env$Site,
  p.adjust.method = "BH"
)

# Shannon diversity vs environmental variables ####

p1 <- ggplot(alpha_env, aes(x = NO3_NO2, y = Shannon)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_point(aes(shape = Sample_type, fill = Site), size = 5, color = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = a_coral_colors) +
  guides(
    shape = guide_legend(override.aes = list(fill = "#00000000")),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_classic(base_size = 14) +
  labs(x = expression("Nitrate + Nitrite (µmol/L)"), y = "Shannon Diversity")

p2 <- ggplot(alpha_env, aes(x = PO4, y = Shannon)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_point(aes(shape = Sample_type, fill = Site), size = 5, color = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = a_coral_colors) +
  guides(
    shape = guide_legend(override.aes = list(fill = "#00000000")),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_classic(base_size = 14) +
  labs(x = expression("Phosphate (µmol/L)"), y = "Shannon Diversity")

p3 <- ggplot(alpha_env, aes(x = SiO2, y = Shannon)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_point(aes(shape = Sample_type, fill = Site), size = 5, color = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = a_coral_colors) +
  guides(
    shape = guide_legend(override.aes = list(fill = "#00000000")),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_classic(base_size = 14) +
  labs(x = expression("Silicon dioxide"), y = "Shannon Diversity")

p4 <- ggplot(alpha_env, aes(x = Salinity, y = Shannon)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_point(aes(shape = Sample_type, fill = Site), size = 5, color = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = a_coral_colors) +
  guides(
    shape = guide_legend(override.aes = list(fill = "#00000000")),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_classic(base_size = 14) +
  labs(x = "Salinity", y = "Shannon Diversity")

p5 <- ggplot(alpha_env, aes(x = Arag, y = Shannon)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_point(aes(shape = Sample_type, fill = Site), size = 5, color = "black") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = a_coral_colors) +
  guides(
    shape = guide_legend(override.aes = list(fill = "#00000000")),
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = expression("Aragonite saturation state (" * Omega[arag] * ")"),
    y = "Shannon Diversity"
  )

# Apply the Helvetica formatting used in the working figures.
p1 <- p1 + theme_cowplot(font_family = "Helvetica")
p2 <- p2 + theme_cowplot(font_family = "Helvetica")
p3 <- p3 + theme_cowplot(font_family = "Helvetica")
p4 <- p4 + theme_cowplot(font_family = "Helvetica")
p5 <- p5 + theme_cowplot(font_family = "Helvetica")

legend <- get_legend(p1 + theme(legend.position = "right"))

combined_plot <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  p4 + theme(legend.position = "none"),
  p5 + theme(legend.position = "none"),
  legend,
  labels = c("A", "B", "C", "D", "E", ""),
  ncol = 3
)

combined_plot

ggsave(
  file.path(results_dir, "kahekili_16S_shannon_environmental_relationships.pdf"),
  plot = combined_plot,
  width = 14,
  height = 8
)

# Kendall correlations ####

kendall_test <- function(data, variable) {
  test <- cor.test(
    data$Shannon,
    data[[variable]],
    method = "kendall",
    exact = FALSE
  )

  tibble(
    Variable = variable,
    Tau = unname(test$estimate),
    P_value = test$p.value
  )
}

environmental_variables <- c("NO3_NO2", "PO4", "SiO2", "Salinity", "Arag")

kendall_results <- map_dfr(
  environmental_variables,
  ~ kendall_test(alpha_env, .x)
)

write.csv(
  kendall_results,
  file.path(results_dir, "kahekili_16S_kendall_correlations.csv"),
  row.names = FALSE
)

kendall_results

# Correlations by site, matching the variables included in the working script.
cor_results_by_site <- alpha_env %>%
  group_by(Site) %>%
  summarize(
    correlation_NO3 = cor(Shannon, NO3_NO2, method = "kendall", use = "complete.obs"),
    correlation_PO4 = cor(Shannon, PO4, method = "kendall", use = "complete.obs"),
    correlation_SiO2 = cor(Shannon, SiO2, method = "kendall", use = "complete.obs"),
    correlation_Salinity = cor(Shannon, Salinity, method = "kendall", use = "complete.obs"),
    .groups = "drop"
  )

write.csv(
  cor_results_by_site,
  file.path(results_dir, "kahekili_16S_kendall_correlations_by_site.csv"),
  row.names = FALSE
)

cor_results_by_site %>%
  gt() %>%
  tab_header(
    title = "Kendall's Tau Correlation Between Shannon Diversity and Environmental Variables"
  ) %>%
  fmt_number(
    columns = c(correlation_NO3, correlation_PO4, correlation_SiO2, correlation_Salinity),
    decimals = 3
  ) %>%
  cols_label(
    Site = "Site",
    correlation_NO3 = "NO₃⁻ + NO₂⁻",
    correlation_PO4 = "PO₄³⁻",
    correlation_SiO2 = "SiO2",
    correlation_Salinity = "Salinity"
  )

# Environmental differences among sites ####
env <- data.frame(sample_data(rare_mauiseq)) %>%
  filter(Site != "well")

kruskal.test(NO3_NO2 ~ Site, data = env)
kruskal.test(CO3 ~ Site, data = env)

pairwise.wilcox.test(
  env$NO3_NO2,
  env$Site,
  p.adjust.method = "BH"
)
pairwise.wilcox.test(
  env$CO3,
  env$Site,
  p.adjust.method = "BH"
)
