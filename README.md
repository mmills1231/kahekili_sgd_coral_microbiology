# Kahekili Coral Microbiome

This repository contains the R scripts and analysis workflow supporting the manuscript: **Relative Taxonomic Stability of Coral-Associated and Reef-Water Microbial Communities Does Not Preclude Predicted Functional Variation Near Submarine Groundwater Discharge in the Kahekili Herbivore Fisheries Management Area, Maui, Hawaiʻi**, which is currently in preparation for publication.

---

## Overview

This study assessed coral mucus-associated and reef-water prokaryotic communities across a submarine groundwater discharge (SGD) gradient at Kahekili Beach, Maui, Hawaiʻi, USA. Analyses examined microbial diversity, community composition, environmental relationships, and predicted functional potential across coral, mucus, and seawater samples.

---

## Data availability

Raw 16S rRNA gene amplicon sequencing data are publicly available through the **NCBI Sequence Read Archive (SRA)** under **BioProject PRJNA1534070**.

This repository contains the R scripts and processed input data required to reproduce the analyses presented in the manuscript, including ASV count tables, taxonomy tables, sample metadata, and functional prediction outputs. Scripts are numbered and intended to be run sequentially within each analysis folder.

---

## Software

### R using the following packages:

- phyloseq
- tidyverse
- ggplot2
- vegan
- plotly
- patchwork
- cowplot
- gt

### Tax4Fun2

Predicted prokaryotic functional potential was assessed using **Tax4Fun2**.
