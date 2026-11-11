# Kahekili 16S analysis: sequencing-depth quality control
# Author: Makeda I. Mills

source(file.path("scripts", "01_16S_analysis", "02_16S_setup.R"))

# Sequencing depth before rarefaction ####

asv_matrix <- as(otu_table(mauiseq16S), "matrix")
if (taxa_are_rows(otu_table(mauiseq16S))) {
  asv_matrix <- t(asv_matrix)
}

sample_depths <- sample_sums(mauiseq16S)

qc_summary <- tibble(
  SampleID = names(sample_depths),
  Reads = as.numeric(sample_depths)
)

write.csv(
  qc_summary,
  file.path(results_dir, "kahekili_16S_sequencing_depth.csv"),
  row.names = FALSE
)

print(summary(sample_depths))
cat("Total reads before rarefaction:", sum(sample_sums(mauiseq16S)), "\n")
cat("Total reads after rarefaction to 5,000 reads/sample:",
    sum(sample_sums(rare_mauiseq)), "\n")

# Rarefaction curves ####

pdf(file.path(results_dir, "kahekili_16S_rarefaction_curves.pdf"),
    width = 8, height = 6)
rarecurve(
  asv_matrix,
  step = 50,
  cex = 0.5,
  xlab = "Sequencing depth (reads)",
  ylab = "Observed ASV richness"
)
dev.off()

# Inspect the rarefied phyloseq object.
rare_mauiseq
