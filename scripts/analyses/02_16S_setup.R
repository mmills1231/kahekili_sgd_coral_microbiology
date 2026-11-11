# Kahekili 16S analysis: shared setup
# Author: Makeda I. Mills, adapted from https://rpubs.com/lgschaerer/pstutorial 

# Run scripts from the repository root. Scripts 03–08 source this file.

require(phyloseq)
require(tidyverse)
require(vegan)
require(ggplot2)

# File paths ####

data_dir <- "input_data"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

asv_file <- file.path(data_dir, "kahekili_16S_asv_table.csv")
tax_file <- file.path(data_dir, "kahekili_16S_taxonomy.csv")
metadata_file <- file.path(data_dir, "kahekili_metadata.csv")

# Read data ####

asv <- read.csv(
  asv_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

tax <- read.csv(
  tax_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  metadata_file,
  row.names = 1,
  header = TRUE,
  check.names = TRUE,
  stringsAsFactors = TRUE
)

# Keep the study sites represented in the original working analysis.
study_sites <- c("control", "S1", "S2", "S3", "S4", "vent", "well")

meta_filtered <- meta %>%
  rownames_to_column(var = "SampleID") %>%
  filter(Site %in% study_sites) %>%
  column_to_rownames(var = "SampleID")

# Create phyloseq object ####

ASV <- otu_table(as.matrix(t(asv)), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
SAM <- sample_data(meta_filtered)

mauiseq16S <- phyloseq(ASV, TAX, SAM)

# Retain bacterial and archaeal ASVs while removing mitochondrial/chloroplast
# assignments. The working script used the same intended filter but without
# parentheses, which could change logical precedence; the intended logic is
# written explicitly here.
mauiseq16S <- subset_taxa(
  mauiseq16S,
  Kingdom %in% c("Bacteria", "Archaea") &
    Family != "Mitochondria" &
    Order != "Chloroplast"
)

# Rarefied objects used in diversity analyses ####

rare_mauiseq <- rarefy_even_depth(
  mauiseq16S,
  sample.size = 5000,
  rngseed = 123,
  verbose = FALSE
)

coral_physeq <- subset_samples(mauiseq16S, Sample_type == "coral")

# Preserved from the working analysis. The original note indicates that 86 was
# selected as the lowest coral-sample depth by manual inspection. Verify this
# depth before using rare_coral in a final manuscript analysis.
rare_coral <- rarefy_even_depth(
  coral_physeq,
  sample.size = 86,
  rngseed = 123,
  verbose = FALSE
)

# Plotting palettes used in the working analysis ####

a_coral_colors <- c(
  "control" = "#FF6F61",
  "S1" = "#40E0D0",
  "S2" = "#2E8B57",
  "S3" = "#4169E1",
  "S4" = "#8E44AD",
  "vent" = "#DAA520",
  "well" = "#FF1493"
)

species_colors <- c(
  "Porites compressa" = "#FF69B4",
  "Porites lobata" = "#F4A460",
  "Pocillopora meandrina" = "#2E8B57"
)

sampletype_colors <- c(
  "water" = "#0072B2",
  "coral" = "#D9534F",
  "mucus" = "#4CAF50"
)

coral_colors3_20 <- c(
  "Acinetobacter" = "#20B2AA",
  "Alteromonas" = "#1E90FF",
  "Ascidiaceihabitans" = "#48D1CC",
  "Clade_Ia" = "#5F9EA0",
  "Cobetia" = "#FFDAB9",
  "Endozoicomonas" = "#C71585",
  "Erythrobacter" = "#87CEEB",
  "Halomonas" = "#BC8F8F",
  "Mesonia" = "#FFEFD5",
  "Muricauda" = "#FFF5EE",
  "NS4_marine_group" = "#4169E1",
  "NS5_marine_group" = "#00BFFF",
  "Prochlorococcus_MIT9313" = "#0047AB",
  "Pseudoalteromonas" = "#E0B0FF",
  "Pseudomonas" = "#004040",
  "Psychrobacter" = "#E15100",
  "Ruegeria" = "#800020",
  "Tenacibaculum" = "#D6A4A4",
  "Thalassotalea" = "#CD9575",
  "Vibrio" = "#2E8B57",
  "Others" = "#9D9D9D"
)
