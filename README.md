# COVID-19 respiratory microbiome analysis code

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/vllorens/covid19_respiratory_microbiome/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/409586241.svg)](https://zenodo.org/badge/latestdoi/409586241)

Code accompanying the manuscript "Clinical practices underlie COVID-19 patient respiratory microbiome composition and its interactions with the host". Link to preprint: https://www.medrxiv.org/content/10.1101/2020.12.23.20248425v3

## Contributors: 

* Verónica Lloréns Rico, PhD
* Ann C. Gregory, PhD


## File description:

### Code files:

* script_covid_preprocessing.R: raw data preprocessing
* script_covid_dataexploration.R: initial data exploration and plots 
* script_covid_alphadiv.R: Alpha diversity analyses and glmm modeling
* script_covid_betadiv.R: dbRDA analyses of microbiome composition
* script_covid_diffabundances.R: differential taxon abundance analyses
* script_covid_species_strain_analyses.R: macro-vs-micro diversity analyses
* script_covid_sc_analyses.R: analyses on the scRNA-seq data of the lower respiratory tract cohort

### Additional folders:

* data/: use this folder to download the raw data and metadata tables from EGA (with controlled access, accession number EGAS00001004951). Contains two additional files:
  - sample_ids.txt: contains sample IDs to preprocess the data (used in script_covid_dataexploration.R)
  - coding_table_final.txt: contains explanations of the metadata variables (used in script_covid_betadiv.R)
* R/: additional R functions, called by the different scripts


## Installation instructions

### Hardware required

All code was run on laptop/desktop computers with 8 cores with 16GB RAM

### Software required

Required: R (https://cran.r-project.org/). The code in this repository was run on R v.4.0.5. 
Recommended: Rstudio (https://www.rstudio.com/)

R packages required:

* dada2 (v1.18.0)
* phyloseq (v1.34.0)
* ggplot2 (v3.3.3)
* ggpubr (v0.4.0)
* cowplot (v1.1.1)
* tidyverse (v1.3.1)
* reshape2 (v1.4.4)
* compositions (v2.0.1)
* vegan (v2.5.7)
* rstatix (v0.7.0)
* ALDEX2 (v1.22.0)
* ggrepel (v0.9.1)
* tibble (v3.1.2)
* DECIPHER (v2.18.1)
* Biostrings (v2.58.0)
* biomod2 (v3.5.1)
* wesanderson (v0.3.6)
* colortools (v0.1.5)
* ggiraph (v0.7.10)
* ggiiraphExtra (v0.3.0)
* glmulti (v1.0.8)
* sjPlot (v2.8.7)
* lme4 (v1.1.27)
* CoDaSeq (v0.99.6)
* DESeq2 (v1.30.1)
* mixOmics (v6.14.1)
* lubridate (v1.7.10)
* Seurat (v4.0.4)
* chisq.posthoc.test (v0.1.2)
* ggmosaic (v0.3.3)
* patchwork (v1.1.1)
* gtable(v0.3.0)


## Running instructions

1. Download raw sequencing files and metadata from EGA: [link to study](https://ega-archive.org/studies/EGAS00001004951)
2. Store them in the `data/` folder, or any other directory of your choice.
3. Run `script_covid_preprocessing.R`. This will do all the raw data preprocessing: quality control, trimming, denoising, ASV assignation and taxonomic annotation, as well as decontamination. 
4. Run `script_covid_dataexploration.R`. This will do the initial data exploration, and generate some of the plots available in Figure 1. 
5. Run `script_covid_alphadiv.R`. This will perform the alpha diversity analyses and modeling section. 
6. Run `script_covid_betadiv.R`. This will perform the beta diversity analyses and modeling section. 
7. Run `script_covid_diffabundances.R`. This will perform the differential taxon abundances reported in the manuscript.
8. Run `script_covid_species_strain_analyses.R`. This will perform the analyses on macro- and micro-diversity, shown in Figure 2. 
9. Download and process raw sequencing files from the scRNA-seq lower respiratory tract cohort. [link to study](https://ega-archive.org/studies/EGAS00001004717) and [link to study website](https://lambrechtslab.sites.vib.be/en/immune-atlas)
10. Run `script_covid_sc_analyses.R`. This will perform the analyses on the lower respiratory tract cohort, shown in Figure 3. 

