# COVID-19 respiratory microbiome analysis code
Code accompanying the manuscript "Clinical practices underlie COVID-19 patient respiratory microbiome composition and its interactions with the host"

Contributors: 

* Verónica Lloréns Rico, PhD
* Ann C. Gregory, PhD


Code description:

* script_covid_preprocessing.R: raw data preprocessing
* script_covid_dataexploration.R: initial data exploration and plots 
* script_covid_alphadiv.R: Alpha diversity analyses and glmm modeling
* script_covid_betadiv.R: dbRDA analyses of microbiome composition
* script_covid_diffabundances.R: differential taxon abundance analyses
* script_covid_species_strain_analyses.R: macro-vs-micro diversity analyses

Additional folders:

* data/: use this folder to download the raw data and metadata tables from EGA (with controlled access, accession number EGAS00001004951). Contains two additional files:
  - sample_ids.txt: contains sample IDs to preprocess the data (used in script_covid_dataexploration.R)
  - coding_table_final.txt: contains explanations of the metadata variables (used in script_covid_betadiv.R)
* R/: additional R functions, called by the different scripts

R packages required:

* dada2
* phyloseq2
* ggplot2
* ggpubr
* cowplot
* tidyverse
* reshape2
* compositions
* vegan
* rstatix
* GLDEX
* ggrepel
* tibble
* DECIPHER
* Biostrings
* biomod2
* wesanderson
* colortools
* ggigraph
* ggíigraphextra
* glmulti
* sjPlot
* lme4
* CoDaSeq
* DESeq2
* mixOmics
* lubridate
