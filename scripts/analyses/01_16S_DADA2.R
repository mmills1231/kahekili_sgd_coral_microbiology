# Kahekili 16S analysis: DADA2 sequence processing and taxonomy assignment
# Author: Makeda I. Mills, adapted from https://benjjneb.github.io/dada2/tutorial.html 

require(dada2)

# File paths ####

# Directory containing the paired-end raw FASTQ files.
# Update this path if raw reads are stored elsewhere.
raw_data_dir <- "raw_data/16S"

# SILVA training set used for taxonomy assignment in the original analysis.
# Update this path to the local location of silva_nr_v132_train_set.fa.
silva_train_set <- "reference_data/silva_nr_v132_train_set.fa"

output_dir <- "input_data"
filtered_dir <- file.path(raw_data_dir, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)

# Locate forward and reverse reads ####

fnFs <- sort(list.files(
  raw_data_dir,
  pattern = "_R1_001.fastq",
  full.names = TRUE
))

fnRs <- sort(list.files(
  raw_data_dir,
  pattern = "_R2_001.fastq",
  full.names = TRUE
))

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

# Inspect representative read-quality profiles before filtering.
# Adjust the index range if fewer than 60 samples are present.
plotQualityProfile(fnFs[seq_len(min(20, length(fnFs)))])
plotQualityProfile(fnRs[seq_len(min(20, length(fnRs)))])

# Filter and trim reads ####

filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(250, 225),
  maxN = 0,
  maxEE = c(2, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

head(out)

# Learn error rates and infer ASVs ####

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

derepFs <- lapply(filtFs, derepFastq, verbose = TRUE)
derepRs <- lapply(filtRs, derepFastq, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired reads and construct sequence table ####

mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = TRUE
)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

# Retain sequence lengths used in the original analysis.
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% seq(250, 445)]

# Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(
  seqtab2,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab2)

# -----------------------------------------------------------------------------
# Track reads through the DADA2 workflow
# -----------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)
rownames(track) <- sample.names

write.csv(
  track,
  file.path(output_dir, "kahekili_16S_read_tracking.csv")
)

# DADA2 produces samples as rows and ASVs as columns. The downstream phyloseq
# script transposes this table when constructing the OTU table.
write.csv(
  seqtab.nochim,
  file.path(output_dir, "kahekili_16S_asv_table.csv")
)

# Assign taxonomy ####

taxa <- assignTaxonomy(
  seqtab.nochim,
  silva_train_set,
  multithread = TRUE
)

write.csv(
  taxa,
  file.path(output_dir, "kahekili_16S_taxonomy.csv")
)
