# Kahekili 16S Tax4Fun2 analysis
# Author: Makeda I. Mills

# Packages ####
require(phyloseq)
require(tidyverse)
require(seqinr)
require(Tax4Fun2)

# Load 16S data ####
# Run this script from the repository root.
source("scripts/analyses/02_16S_setup.R")

# Tax4Fun2 paths ####
tax4fun_reference <- "Tax4Fun2_ReferenceData_v2"
tax4fun_temp <- "kahekili_tax4fun2_ref99nr"

# Filter phyloseq object ####

mauiseq16S <- subset_taxa(
  mauiseq16S,
  Kingdom %in% c("Bacteria", "Archaea") &
    (is.na(Family) | Family != "Mitochondria") &
    (is.na(Order) | Order != "Chloroplast")
)

# Extract filtered ASV count table ####

asv_tax4fun <- as(otu_table(mauiseq16S), "matrix")

# Tax4Fun2 wants ASVs as rows and samples as columns
if (!taxa_are_rows(mauiseq16S)) {
  asv_tax4fun <- t(asv_tax4fun)
}

dim(asv_tax4fun)

head(rownames(asv_tax4fun))
head(colnames(asv_tax4fun))

# Check
sample_names(mauiseq16S)

identical(
  colnames(asv_tax4fun),
  sample_names(mauiseq16S)
)

# Remove absent ASVs ####

sum(rowSums(asv_tax4fun) == 0)

asv_tax4fun <- asv_tax4fun[
  rowSums(asv_tax4fun) > 0,
  ,
  drop = FALSE
]

dim(asv_tax4fun)

# Load ASV sequence reference ####

seq_ref <- read.csv(
  "input_data/kahekili_16S_asv_sequences.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(seq_ref)

head(rownames(seq_ref))
head(seq_ref)

# Match ASVs to sequence reference ####

missing_sequences <- setdiff(
  rownames(asv_tax4fun),
  rownames(seq_ref)
)

length(missing_sequences)

missing_sequences


# Reorder sequences to match the ASV table
seq_ref_tax4fun <- seq_ref[
  rownames(asv_tax4fun),
  ,
  drop = FALSE
]

# Check
nrow(asv_tax4fun)
nrow(seq_ref_tax4fun)

identical(
  rownames(asv_tax4fun),
  rownames(seq_ref_tax4fun)
)

# Write filtered FASTA ####

seq_names <- rownames(seq_ref_tax4fun)

seq_list <- as.list(
  as.character(seq_ref_tax4fun[, 1])
)

write.fasta(
  sequences = seq_list,
  names = seq_names,
  file.out = "kahekili_16S_tax4fun2_filtered.fasta"
)

# Verify
kahekili_fasta_tax4fun <- read.fasta(
  "kahekili_16S_tax4fun2_filtered.fasta",
  seqtype = "DNA",
  as.string = TRUE
)

length(kahekili_fasta_tax4fun)

identical(
  names(kahekili_fasta_tax4fun),
  rownames(asv_tax4fun)
)

# Write filtered ASV abundance table ####

write.table(
  asv_tax4fun,
  file = "kahekili_16S_tax4fun2_asv_table.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Tax4Fun2 input QC ####

cat("Samples:", ncol(asv_tax4fun), "\n")
cat("ASVs:", nrow(asv_tax4fun), "\n")
cat("Zero-sum ASVs:", sum(rowSums(asv_tax4fun) == 0), "\n")
cat(
  "ASV IDs match FASTA:",
  identical(rownames(asv_tax4fun), names(kahekili_fasta_tax4fun)),
  "\n"
)

# Show samples included in final analysis
colnames(asv_tax4fun)

# Configure NCBI BLAST+ ####

# Tax4Fun2 requires NCBI BLAST+.
# If BLAST+ is not already on your PATH, set BLAST_BIN to the directory
# containing blastn and makeblastdb before running this script.
blast_bin <- Sys.getenv("BLAST_BIN")

if (nzchar(blast_bin)) {
  Sys.setenv(
    PATH = paste(blast_bin, Sys.getenv("PATH"), sep = .Platform$path.sep)
  )
}

if (Sys.which("blastn") == "" || Sys.which("makeblastdb") == "") {
  stop(
    "NCBI BLAST+ was not found. Add BLAST+ to PATH or set BLAST_BIN before running Tax4Fun2."
  )
}

system("blastn -version")

# Run RefBlast

runRefBlast(
  path_to_otus = "kahekili_16S_tax4fun2_filtered.fasta",
  path_to_reference_data = tax4fun_reference,
  path_to_temp_folder = tax4fun_temp,
  database_mode = "Ref99NR",
  use_force = TRUE,
  num_threads = 6)

# Functional prediction

makeFunctionalPrediction(
  path_to_otu_table = "kahekili_16S_tax4fun2_asv_table.txt",
  path_to_reference_data = tax4fun_reference,
  path_to_temp_folder = tax4fun_temp,
  database_mode = "Ref99NR",
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = 0.97,
  normalize_pathways = FALSE
)

# Check Tax4Fun2 samples against metadata ####

# Samples in Tax4Fun2 table but NOT in filtered metadata
setdiff(
  colnames(asv_tax4fun),
  rownames(meta3)
)

# Samples in filtered metadata but NOT in Tax4Fun2 table
setdiff(
  rownames(meta3),
  colnames(asv_tax4fun)
)

# Match metadata to Tax4Fun2 samples ####

meta_tax4fun <- meta3[
  colnames(asv_tax4fun),
  ,
  drop = FALSE
]

# Check
dim(meta_tax4fun)

identical(
  rownames(meta_tax4fun),
  colnames(asv_tax4fun)
)

# Extract Tax4Fun2 prediction coverage ####

tax4fun_log <- readLines(file.path(tax4fun_temp, "logfile2.txt"))

# Find the two sections
otu_start <- grep(
  "amount of otus unused",
  tax4fun_log,
  ignore.case = TRUE
)

seq_start <- grep(
  "amount of sequences unused",
  tax4fun_log,
  ignore.case = TRUE
)

# Function to parse "Sample: value"
parse_unknown_fraction <- function(x) {

  x <- x[grepl(": ", x)]

  split_x <- strsplit(x, ": ")

  data.frame(
    Sample = sapply(split_x, `[`, 1),
    Fraction = as.numeric(sapply(split_x, `[`, 2)),
    stringsAsFactors = FALSE
  )
}

# Fraction of ASVs unused
FTU <- parse_unknown_fraction(
  tax4fun_log[(otu_start + 1):(seq_start - 1)]
)

# Fraction of sequence abundance unused
FSU <- parse_unknown_fraction(
  tax4fun_log[(seq_start + 1):length(tax4fun_log)]
)

# Combine
tax4fun_qc <- merge(
  FTU,
  FSU,
  by = "Sample"
)

colnames(tax4fun_qc) <- c(
  "Sample",
  "FTU",
  "FSU"
)

# Fraction of sequence abundance represented
tax4fun_qc$Sequence_Coverage <- 1 - tax4fun_qc$FSU

# Checking
nrow(tax4fun_qc)

setdiff(
  colnames(asv_tax4fun),
  tax4fun_qc$Sample
)

setdiff(
  tax4fun_qc$Sample,
  colnames(asv_tax4fun)
)

# Check Tax4Fun2 sample names ####

# Current ASV table
colnames(asv_tax4fun)[grepl("^V1|^V2", colnames(asv_tax4fun))]

# Written Tax4Fun2 abundance table
tax4fun_input_check <- read.table(
  "kahekili_16S_tax4fun2_asv_table.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

colnames(tax4fun_input_check)[
  grepl("^V1|^V2", colnames(tax4fun_input_check))
]

func_check <- read.table(
  file.path(tax4fun_temp, "functional_prediction.txt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1
)

colnames(func_check)[
  grepl("^V1|^V2", colnames(func_check))
]

# Correct V1 and V2 naming collision ####

tax4fun_qc$Sample[tax4fun_qc$Sample == "V1.y"] <- "V1"
tax4fun_qc$Sample[tax4fun_qc$Sample == "V2.y"] <- "V2"

# Merge Tax4Fun2 QC with metadata ####

meta_tax4fun_df <- meta_tax4fun %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample")

tax4fun_qc_meta <- tax4fun_qc %>%
  left_join(
    meta_tax4fun_df,
    by = "Sample"
  )

nrow(tax4fun_qc_meta)
# Should be 106

sum(is.na(tax4fun_qc_meta$Site))
# Ideally 0

# Create sequence coverage classes ####

tax4fun_qc_meta$Coverage_Class <- cut(
  tax4fun_qc_meta$Sequence_Coverage,
  breaks = c(-Inf, 0.70, 0.80, 0.90, Inf),
  labels = c(
    "<70%",
    "70-80%",
    "80-90%",
    ">=90%"
  ),
  right = FALSE
)

# Check
table(tax4fun_qc_meta$Coverage_Class)

# Sequence coverage class by site ####

table(
  tax4fun_qc_meta$Coverage_Class,
  tax4fun_qc_meta$Site
)

# Percentages

round(
  prop.table(
    table(
      tax4fun_qc_meta$Coverage_Class,
      tax4fun_qc_meta$Site
    ),
    margin = 2
  ) * 100,
  1
)

# Summarize sequence coverage by site ####

tax4fun_qc_meta %>%
  group_by(Site) %>%
  summarise(
    n = n(),
    mean_coverage = mean(Sequence_Coverage, na.rm = TRUE),
    median_coverage = median(Sequence_Coverage, na.rm = TRUE),
    min_coverage = min(Sequence_Coverage, na.rm = TRUE),
    max_coverage = max(Sequence_Coverage, na.rm = TRUE),
    n_below_70 = sum(Sequence_Coverage < 0.70, na.rm = TRUE)
  )

# Test sequence coverage among sites ####

qc_site <- tax4fun_qc_meta %>%
  filter(
    Site %in% c("control", "S1", "S2", "S3", "S4", "vent")
  ) %>%
  droplevels()
# Exclude the well sample because only one well sample is available for comparison.

kruskal.test(
  Sequence_Coverage ~ Site,
  data = qc_site
)
#this test yielded a p-value of 0.4369 ( not significant) and that's a good thing because that means that
# the reference-dtatabase coverage doesn't differ systematically across sites we want to compare

# Sequence coverage by sample type ####

tax4fun_qc_meta %>%
  group_by(Sample_type) %>%
  summarise(
    n = n(),
    mean_coverage = mean(Sequence_Coverage, na.rm = TRUE),
    median_coverage = median(Sequence_Coverage, na.rm = TRUE),
    min_coverage = min(Sequence_Coverage, na.rm = TRUE),
    max_coverage = max(Sequence_Coverage, na.rm = TRUE),
    n_below_70 = sum(Sequence_Coverage < 0.70, na.rm = TRUE)
  )

# Test
kruskal.test(
  Sequence_Coverage ~ Sample_type,
  data = tax4fun_qc_meta
)

#Coverage class
table(
  tax4fun_qc_meta$Coverage_Class,
  tax4fun_qc_meta$Sample_type
)

# Percentages
round(
  prop.table(
    table(
      tax4fun_qc_meta$Coverage_Class,
      tax4fun_qc_meta$Sample_type
    ),
    margin = 2
  ) * 100,
  1
)

# Sequence coverage among coral species ####

coral_qc <- tax4fun_qc_meta %>%
  filter(Sample_type == "coral") %>%
  droplevels()

coral_qc %>%
  group_by(Sample_species) %>%
  summarise(
    n = n(),
    mean_coverage = mean(Sequence_Coverage, na.rm = TRUE),
    median_coverage = median(Sequence_Coverage, na.rm = TRUE),
    min_coverage = min(Sequence_Coverage, na.rm = TRUE),
    max_coverage = max(Sequence_Coverage, na.rm = TRUE)
  )

kruskal.test(
  Sequence_Coverage ~ Sample_species,
  data = coral_qc
)

# Sequence coverage among mucus species ####

mucus_qc <- tax4fun_qc_meta %>%
  filter(Sample_type == "mucus") %>%
  droplevels()

mucus_qc %>%
  group_by(Sample_species) %>%
  summarise(
    n = n(),
    mean_coverage = mean(Sequence_Coverage, na.rm = TRUE),
    median_coverage = median(Sequence_Coverage, na.rm = TRUE),
    min_coverage = min(Sequence_Coverage, na.rm = TRUE),
    max_coverage = max(Sequence_Coverage, na.rm = TRUE)
  )

kruskal.test(
  Sequence_Coverage ~ Sample_species,
  data = mucus_qc
)

# Pairwise coverage comparisons among sample types ####

pairwise.wilcox.test(
  tax4fun_qc_meta$Sequence_Coverage,
  tax4fun_qc_meta$Sample_type,
  p.adjust.method = "BH",
  exact = FALSE
)

# Pairwise species coverage comparisons within mucus ####

pairwise.wilcox.test(
  mucus_qc$Sequence_Coverage,
  mucus_qc$Sample_species,
  p.adjust.method = "BH",
  exact = FALSE
)

# Test site effects on coverage within sample types ####

# Water
water_qc <- tax4fun_qc_meta %>%
  filter(
    Sample_type == "water",
    Site %in% c("control", "S1", "S2", "S3", "S4", "vent")
  ) %>%
  droplevels()

kruskal.test(
  Sequence_Coverage ~ Site,
  data = water_qc
)


# Mucus
mucus_site_qc <- tax4fun_qc_meta %>%
  filter(
    Sample_type == "mucus",
    Site %in% c("control", "S1", "S2", "S3", "S4", "vent")
  ) %>%
  droplevels()

kruskal.test(
  Sequence_Coverage ~ Site,
  data = mucus_site_qc
)


# Coral
coral_site_qc <- tax4fun_qc_meta %>%
  filter(
    Sample_type == "coral",
    Site %in% c("control", "S1", "S2", "S3", "S4", "vent")
  ) %>%
  droplevels()

kruskal.test(
  Sequence_Coverage ~ Site,
  data = coral_site_qc
)

# Check sample distribution
table(
  tax4fun_qc_meta$Site,
  tax4fun_qc_meta$Sample_type
)

# Load Tax4Fun2 predicted functions ####

functional_pred <- read.table(
  file.path(tax4fun_temp, "functional_prediction.txt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

pathway_pred <- read.table(
  file.path(tax4fun_temp, "pathway_prediction.txt"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

# Examine dimensions
dim(functional_pred)
dim(pathway_pred)

# Column names
colnames(functional_pred)
colnames(pathway_pred)

# First few functions/pathways
head(functional_pred[, c(1, ncol(functional_pred))])
head(pathway_pred[, c(1, ncol(pathway_pred))])

# Check functional prediction sample names ####

# Functional prediction should have:
# KO + 106 samples + description

functional_samples <- colnames(functional_pred)[
  !colnames(functional_pred) %in% c("KO", "description")
]

length(functional_samples)

setdiff(
  functional_samples,
  rownames(meta_tax4fun)
)

setdiff(
  rownames(meta_tax4fun),
  functional_samples
)

summary(
  colSums(functional_pred[, functional_samples])
)

# Inspect Tax4Fun2 pathway hierarchy ####

pathway_pred %>%
  select(pathway, level1, level2, level3) %>%
  head(20)

# How many unique categories are in each level?
sapply(
  pathway_pred[, c("level1", "level2", "level3")],
  function(x) length(unique(x))
)

# Look at examples from each hierarchy level
unique(pathway_pred$level1)[1:20]
unique(pathway_pred$level2)[1:20]
unique(pathway_pred$level3)[1:20]

# Build curated functional dictionary ####

functional_dictionary <- tribble(

  ~KO,      ~Function,

# Nitrogen cycling ####

  # Ammonia oxidation potential
  # K10945/K10946 = methane/ammonia monooxygenase components
  # K10535 = hydroxylamine dehydrogenase
  "K10945", "Ammonia oxidation",
  "K10946", "Ammonia oxidation",
  "K10535", "Ammonia oxidation",

  # DNRA - dissimilatory nitrate reduction to ammonium
  "K00370", "DNRA",
  "K00371", "DNRA",
  "K00374", "DNRA",
  "K02567", "DNRA",
  "K00362", "DNRA",
  "K00363", "DNRA",
  "K03385", "DNRA",
  "K15876", "DNRA",

  # Denitrification-related N reduction
  # Note: nosZ / K00376 was not detected,
  # so we are NOT calling this complete denitrification
  "K00368", "Denitrification-related N reduction",
  "K15864", "Denitrification-related N reduction",
  "K04561", "Denitrification-related N reduction",
  "K02305", "Denitrification-related N reduction",

  # Nitrogen fixation - canonical nifHDK
  "K02588", "Nitrogen fixation",
  "K02586", "Nitrogen fixation",
  "K02591", "Nitrogen fixation",

  # Assimilatory nitrate/nitrite reduction
  "K00360", "Assimilatory nitrate reduction",
  "K00367", "Assimilatory nitrate reduction",
  "K00372", "Assimilatory nitrate reduction",
  "K00366", "Assimilatory nitrate reduction",

# Phosphorus cycling ####

  # Inorganic phosphate uptake
  "K02036", "Inorganic phosphate acquisition",
  "K02037", "Inorganic phosphate acquisition",
  "K02038", "Inorganic phosphate acquisition",
  "K02039", "Inorganic phosphate acquisition",
  "K02040", "Inorganic phosphate acquisition",
  "K03306", "Inorganic phosphate acquisition",
  "K03324", "Inorganic phosphate acquisition",
  "K07221", "Inorganic phosphate acquisition",
  "K08176", "Inorganic phosphate acquisition",
  "K16322", "Inorganic phosphate acquisition",

  # P limitation / scavenging machinery
  "K01077", "P scavenging / starvation response",
  "K01113", "P scavenging / starvation response",
  "K06217", "P scavenging / starvation response",
  "K07636", "P scavenging / starvation response",
  "K07657", "P scavenging / starvation response",
  "K07658", "P scavenging / starvation response",

  # Polyphosphate storage / turnover
  "K00937", "Polyphosphate metabolism",
  "K01524", "Polyphosphate metabolism",
  "K22468", "Polyphosphate metabolism",

  # Phosphonate acquisition / utilization
  "K02041", "Phosphonate acquisition / utilization",
  "K02042", "Phosphonate acquisition / utilization",
  "K02043", "Phosphonate acquisition / utilization",
  "K02044", "Phosphonate acquisition / utilization",
  "K03430", "Phosphonate acquisition / utilization",
  "K05781", "Phosphonate acquisition / utilization",
  "K06163", "Phosphonate acquisition / utilization",
  "K09469", "Phosphonate acquisition / utilization",
  "K11081", "Phosphonate acquisition / utilization",
  "K11082", "Phosphonate acquisition / utilization",
  "K11083", "Phosphonate acquisition / utilization",
  "K11084", "Phosphonate acquisition / utilization",
  "K18916", "Phosphonate acquisition / utilization",
  "K21195", "Phosphonate acquisition / utilization",
  "K21196", "Phosphonate acquisition / utilization",


# Carbon cycling ####

  # Inorganic carbon acquisition / fixation
  "K01601", "Inorganic carbon acquisition / fixation",
  "K01602", "Inorganic carbon acquisition / fixation",
  "K01673", "Inorganic carbon acquisition / fixation",
  "K01674", "Inorganic carbon acquisition / fixation",
  "K01726", "Inorganic carbon acquisition / fixation",
  "K11950", "Inorganic carbon acquisition / fixation",
  "K11951", "Inorganic carbon acquisition / fixation",
  "K11952", "Inorganic carbon acquisition / fixation",
  "K11953", "Inorganic carbon acquisition / fixation",

  # Complex carbohydrate degradation / utilization
  "K00688", "Complex carbohydrate utilization",
  "K01183", "Complex carbohydrate utilization",
  "K01187", "Complex carbohydrate utilization",
  "K01188", "Complex carbohydrate utilization",
  "K01210", "Complex carbohydrate utilization",
  "K01215", "Complex carbohydrate utilization",
  "K01222", "Complex carbohydrate utilization",
  "K01223", "Complex carbohydrate utilization",
  "K02438", "Complex carbohydrate utilization",
  "K03791", "Complex carbohydrate utilization",
  "K05349", "Complex carbohydrate utilization",
  "K05350", "Complex carbohydrate utilization",
  "K13381", "Complex carbohydrate utilization",
  "K19668", "Complex carbohydrate utilization",
  "K20547", "Complex carbohydrate utilization",
  "K20829", "Complex carbohydrate utilization",
  "K21571", "Complex carbohydrate utilization",
  "K21572", "Complex carbohydrate utilization",
  "K21573", "Complex carbohydrate utilization",
  "K21574", "Complex carbohydrate utilization",


# Oxidative stress ####

  "K00428", "Oxidative stress defense",
  "K00432", "Oxidative stress defense",
  "K00518", "Oxidative stress defense",
  "K03781", "Oxidative stress defense",
  "K03782", "Oxidative stress defense",
  "K04564", "Oxidative stress defense",
  "K04565", "Oxidative stress defense",
  "K05910", "Oxidative stress defense",
  "K05919", "Oxidative stress defense",
  "K07217", "Oxidative stress defense",
  "K11065", "Oxidative stress defense",
  "K19511", "Oxidative stress defense",


# Sulfur metabolism ####

  "K00380", "Sulfur metabolism",
  "K00381", "Sulfur metabolism",
  "K00385", "Sulfur metabolism",
  "K00392", "Sulfur metabolism",
  "K00394", "Sulfur metabolism",
  "K00395", "Sulfur metabolism",
  "K00956", "Sulfur metabolism",
  "K00957", "Sulfur metabolism",
  "K00958", "Sulfur metabolism",
  "K01011", "Sulfur metabolism",
  "K02439", "Sulfur metabolism",
  "K08352", "Sulfur metabolism",
  "K08353", "Sulfur metabolism",
  "K08354", "Sulfur metabolism",
  "K11180", "Sulfur metabolism",
  "K11181", "Sulfur metabolism",
  "K16936", "Sulfur metabolism",
  "K16937", "Sulfur metabolism",
  "K16950", "Sulfur metabolism",
  "K16951", "Sulfur metabolism",
  "K17218", "Sulfur metabolism",
  "K19713", "Sulfur metabolism"
)

# Refine DNRA category ####

functional_dictionary <- functional_dictionary %>%
  filter(
    !(
      Function == "DNRA" &
        KO %in% c(
          "K00370",
          "K00371",
          "K00374",
          "K02567"
        )
    )
  ) %>%
  mutate(
    Function = ifelse(
      Function == "DNRA",
      "DNRA (nitrite to ammonium)",
      Function
    )
  )

functional_dictionary %>%
  count(Function, sort = TRUE)

nrow(functional_dictionary)

functional_pred %>%
  filter(
    KO %in% c(
      "K28504",
      "K10944",
      "K10945",
      "K10946",
      "K10535"
    )
  ) %>%
  select(KO, description)

# Rename ammonia oxidation category ####

functional_dictionary <- functional_dictionary %>%
  mutate(
    Function = ifelse(
      Function == "Ammonia oxidation",
      "Potential ammonia oxidation",
      Function
    )
  )

functional_dictionary %>%
  filter(Function == "Potential ammonia oxidation")

# Validate curated functional dictionary ####

dictionary_check <- functional_dictionary %>%
  mutate(
    Present = KO %in% functional_pred$KO
  )

table(dictionary_check$Present)

dictionary_check %>%
  filter(!Present)

functional_dictionary %>%
  count(Function, sort = TRUE)

# Add KO descriptions and save dictionary ####

functional_dictionary_full <- functional_dictionary %>%
  left_join(
    functional_pred %>%
      select(KO, description),
    by = "KO"
  )

functional_dictionary_full %>%
  arrange(Function, KO)

write.csv(
  functional_dictionary_full,
  "kahekili_Tax4Fun2_function_dictionary.csv",
  row.names = FALSE
)

# Extract curated KOs ####

functional_long <- functional_pred %>%
  select(
    KO,
    all_of(functional_samples)
  ) %>%
  inner_join(
    functional_dictionary,
    by = "KO"
  ) %>%
  pivot_longer(
    cols = all_of(functional_samples),
    names_to = "Sample",
    values_to = "KO_relative_abundance"
  )

# Sum KOs by functional category ####

functional_by_sample <- functional_long %>%
  group_by(
    Sample,
    Function
  ) %>%
  summarise(
    Predicted_relative_abundance = sum(
      KO_relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

n_distinct(functional_by_sample$Sample)
n_distinct(functional_by_sample$Function)
nrow(functional_by_sample)

# Merge functional predictions with metadata ####

functional_by_sample_meta <- functional_by_sample %>%
  left_join(
    tax4fun_qc_meta %>%
      select(
        Sample,
        Site,
        Sample_type,
        Sample_species,
        Sequence_Coverage
      ),
    by = "Sample"
  )

# Check dimensions
dim(functional_by_sample_meta)

# Check for missing metadata
sum(is.na(functional_by_sample_meta$Site))
sum(is.na(functional_by_sample_meta$Sample_type))

# Preview
head(functional_by_sample_meta)

# Examine functional abundance distributions ####

function_distribution <- functional_by_sample_meta %>%
  group_by(Function) %>%
  summarise(
    n = n(),
    min = min(Predicted_relative_abundance, na.rm = TRUE),
    median = median(Predicted_relative_abundance, na.rm = TRUE),
    mean = mean(Predicted_relative_abundance, na.rm = TRUE),
    max = max(Predicted_relative_abundance, na.rm = TRUE),
    n_zero = sum(Predicted_relative_abundance == 0),
    percent_zero = round(
      100 * n_zero / n,
      1
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(percent_zero))

function_distribution

# Examine zeros by sample type ####

function_distribution_type <- functional_by_sample_meta %>%
  group_by(
    Sample_type,
    Function
  ) %>%
  summarise(
    n = n(),
    median = median(
      Predicted_relative_abundance,
      na.rm = TRUE
    ),
    mean = mean(
      Predicted_relative_abundance,
      na.rm = TRUE
    ),
    n_zero = sum(
      Predicted_relative_abundance == 0
    ),
    percent_zero = round(
      100 * n_zero / n,
      1
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Sample_type,
    desc(percent_zero)
  )

function_distribution_type

# Check very low-abundance functional categories
functional_by_sample_meta %>%
  group_by(Function) %>%
  summarise(
    median = median(Predicted_relative_abundance),
    mean = mean(Predicted_relative_abundance),
    max = max(Predicted_relative_abundance)
  ) %>%
  arrange(mean)

# Calculate Cliff's delta ####

cliffs_delta <- function(x, y) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }

  mean(
    outer(
      x,
      y,
      FUN = function(a, b) sign(a - b)
    )
  )
}

# Calculate control-relative functional effect sizes ####

comparison_sites <- c(
  "S1",
  "S2",
  "S3",
  "S4",
  "vent"
)

functional_effects <- functional_by_sample_meta %>%

  filter(
    Site %in% c("control", comparison_sites)
  ) %>%

  group_by(
    Sample_type,
    Function
  ) %>%

  group_modify(
    ~ {

      control_values <- .x %>%
        filter(Site == "control") %>%
        pull(Predicted_relative_abundance)

      sites_present <- intersect(
        comparison_sites,
        unique(as.character(.x$Site))
      )

      bind_rows(
        lapply(
          sites_present,
          function(current_site) {

            site_values <- .x %>%
              filter(Site == current_site) %>%
              pull(Predicted_relative_abundance)

            tibble(
              Site = current_site,

              n_control = length(control_values),
              n_site = length(site_values),

              Control_mean = mean(control_values),
              Site_mean = mean(site_values),

              Control_median = median(control_values),
              Site_median = median(site_values),

              Control_prevalence =
                mean(control_values > 0),

              Site_prevalence =
                mean(site_values > 0),

              Cliff_delta =
                cliffs_delta(
                  site_values,
                  control_values
                ),

              p_value =
                wilcox.test(
                  site_values,
                  control_values,
                  exact = FALSE
                )$p.value
            )
          }
        )
      )
    }
  ) %>%

  ungroup()

# Adjust p-values ####

functional_effects <- functional_effects %>%
  group_by(Sample_type) %>%
  mutate(
    p_adj = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup()

functional_effects <- functional_effects %>%
  mutate(
    Direction = case_when(
      Cliff_delta > 0 ~ "Higher than Control",
      Cliff_delta < 0 ~ "Lower than Control",
      TRUE ~ "No difference"
    )
  )

functional_effects %>%
  arrange(
    Sample_type,
    Site,
    desc(abs(Cliff_delta))
  )

functional_effects <- functional_effects %>%
  mutate(
    Sparse = case_when(
      Control_prevalence < 0.50 |
        Site_prevalence < 0.50 ~ TRUE,
      TRUE ~ FALSE
    )
  )

functional_effects %>%
  filter(Sparse) %>%
  select(
    Sample_type,
    Site,
    Function,
    n_control,
    n_site,
    Control_prevalence,
    Site_prevalence,
    Cliff_delta,
    p_adj
  )

# Check coral species distribution among sites ####

tax4fun_qc_meta %>%
  filter(
    Sample_type == "mucus",
    Site %in% c("control", "S1", "S2", "S3", "S4")
  ) %>%
  count(
    Site,
    Sample_species
  ) %>%
  tidyr::pivot_wider(
    names_from = Sample_species,
    values_from = n,
    values_fill = 0
  )

tax4fun_qc_meta %>%
  filter(
    Sample_type == "coral",
    Site %in% c("control", "S1", "S2", "S3", "S4")
  ) %>%
  count(
    Site,
    Sample_species
  ) %>%
  tidyr::pivot_wider(
    names_from = Sample_species,
    values_from = n,
    values_fill = 0
  )

# Determine compartment-specific function prevalence ####

function_eligibility <- bind_rows(

  # Water: includes vent
  functional_by_sample_meta %>%
    filter(
      Sample_type == "water",
      Site %in% c("control", "S1", "S2", "S3", "S4", "vent")
    ) %>%
    group_by(Function, Site) %>%
    summarise(
      prevalence = mean(Predicted_relative_abundance > 0),
      .groups = "drop"
    ) %>%
    group_by(Function) %>%
    summarise(
      min_prevalence = min(prevalence),
      .groups = "drop"
    ) %>%
    mutate(
      Sample_type = "water",
      Eligible = min_prevalence >= 0.50
    ),

  # Mucus: no vent
  functional_by_sample_meta %>%
    filter(
      Sample_type == "mucus",
      Site %in% c("control", "S1", "S2", "S3", "S4")
    ) %>%
    group_by(Function, Site) %>%
    summarise(
      prevalence = mean(Predicted_relative_abundance > 0),
      .groups = "drop"
    ) %>%
    group_by(Function) %>%
    summarise(
      min_prevalence = min(prevalence),
      .groups = "drop"
    ) %>%
    mutate(
      Sample_type = "mucus",
      Eligible = min_prevalence >= 0.50
    ),

  # Coral: no vent
  functional_by_sample_meta %>%
    filter(
      Sample_type == "coral",
      Site %in% c("control", "S1", "S2", "S3", "S4")
    ) %>%
    group_by(Function, Site) %>%
    summarise(
      prevalence = mean(Predicted_relative_abundance > 0),
      .groups = "drop"
    ) %>%
    group_by(Function) %>%
    summarise(
      min_prevalence = min(prevalence),
      .groups = "drop"
    ) %>%
    mutate(
      Sample_type = "coral",
      Eligible = min_prevalence >= 0.50
    )
)

function_eligibility %>%
  arrange(
    Sample_type,
    Eligible,
    min_prevalence
  )

functional_effects_plot <- functional_effects %>%
  left_join(
    function_eligibility %>%
      select(
        Sample_type,
        Function,
        Eligible
      ),
    by = c(
      "Sample_type",
      "Function"
    )
  )

functional_effects_plot %>%
  filter(!Eligible) %>%
  distinct(
    Sample_type,
    Function
  ) %>%
  arrange(
    Sample_type,
    Function
  )

# Prepare water functional effect plot ####

water_effects_plot <- functional_effects_plot %>%
  filter(
    Sample_type == "water",
    Eligible,
    Site %in% c("vent", "S1", "S2", "S3", "S4")
  ) %>%
  filter(Direction != "No difference") %>%
  mutate(
    Site = factor(
      Site,
      levels = c("vent", "S1", "S2", "S3", "S4")
    ),
    Direction = factor(
      Direction,
      levels = c("Higher than Control", "Lower than Control")
    )
  )

# Plot water functional response relative to control ####

ggplot(
  water_effects_plot,
  aes(
    x = Cliff_delta,
    y = Function,
    fill = Direction
  )
) +
  geom_col(width = 0.75) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  facet_wrap(
    ~ Site,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Higher than Control" = "#2fbf3a",
      "Lower than Control" = "#f3776b"
    )
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  labs(
    x = "Cliff's delta relative to Control",
    y = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

