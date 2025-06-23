# Cervicovaginal Microbiome and Virome During ART and HIV Discordant Shedding  
This repository contains code for the manuscript **Longitudinal cervicovaginal microbiome and virome alterations during ART and discordant shedding in women living with HIV**, under consideration at _Nature Communications_. Link to preprint: [https://www.researchsquare.com/article/rs-4078561/v1](https://www.researchsquare.com/article/rs-4078561/v1)
***

### Contents:  
**1. Bacterial microbiome sequencing analysis:**
- Bacterial_QC_classification.sh: Commands and instructions for QC of raw bacterial metagenomic sequencing data and classification of clean reads.
- parse_krakenuniqreports_mpa.py: Python script to parse combined KrakenUniq outputs to desired taxonomic level.
- Bacterial_diversity_taxonomic.R: R script for ecological and statistical analysis of bacterial taxonomic profiles.
- Bacterial_diversity_functional.R: R script for ecological and statistical analysis of bacterial functional profiles.

**2. DNA virome sequencing analysis:**
- DNA_virome_QC_contigs.sh: Commands and instructions for QC of DNA virome (MDA) metagenomic sequencing data; contig assembly; viral contig identification and classification; mapping clean reads back to viral contigs.
- DNA_virome_analysis.R: R script for ecological and statistical analysis of virome profiles.
- Anellovirus_phylogeny.sh: Commands and instructions for building the anellovirus phylogeny.
- Anellovirus_analysis.R: R script for downstream analyses of anellovirus diversity.
- Papillomavirus_phylogeny.sh: Commmands and instructions for building the papillomavirus phylogeny.
- Papillomavirus_analysis.R: R script for downstream analysis of papillomavirus diversity.

**3. DNA/RNA virome sequencing analysis:**
- DNA-RNA_Virome_QC_blastx.sh: Commands and instructions for QC of DNA/RNA virome (SIA) metagenomic sequencing data and classification of clean reads.
