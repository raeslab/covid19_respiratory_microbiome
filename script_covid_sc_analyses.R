# scRNA-seq analyses
# Lower respiratory microbiome in COVID-19 patients

# load packages and libraries
library(ggpubr)
library(tidyverse)
library(gtable)
library(grid)
library(ggmosaic)
library(rstatix)
library(RColorBrewer)
library(chisq.posthoc.test)
library(wesanderson)
library(patchwork)
source("R/corrmod.R")
source("R/DoMultiBarHeatmap.R")


#### 1. Load data ####

# The processed data was first publish in a separate publication (https://www.nature.com/articles/s41422-020-00455-9)
# Raw data can be accessed with controlled access from the EGA repository. 
# Processed data (i.e. the Seurat object) shall be requested to the lead authors of that study.

# read microbial data [microbial species and associated barcodes obtained after analyzing the raw reads]
microbes <- read_tsv("data/barcode_species_table.tsv")

# read host data [host cell barcodes and annotations, extracted from the Seurat object from the original publication]
host <- read_tsv("data/host_metadata_barcodes.txt")

# rearrange host data
barcode_only <- sapply(host$Barcode, function(X){strsplit(X, split="[[:punct:]]")[[1]][2]})
host <- host %>% mutate(barcode_only)

# join both datasets
microbes_host <- microbes %>% 
  left_join(host, by=c("barcode" = "barcode_only", "sample" = "SampleID")) %>% 
  filter(!is.na(Barcode))

# distribution of host cells with no microbes [background]
host_nomicrobes <- host %>% 
  dplyr::filter(!Barcode %in% microbes_host$Barcode)

# distribution of host cells with microbes (for some calculations the number of bacteria (2108) is used;
# for others the number of unique cells (1703) is used)
host_microbes <- host %>% 
  dplyr::filter(Barcode %in% microbes_host$Barcode)

# put all together for analyses and make a new column indicating whether the barcodes have bacteria or not
full_dataset <- bind_rows(microbes_host %>% mutate(Class="Yes"),
                          host_nomicrobes %>% mutate(Class="No"))
full_dataset$Class <- factor(full_dataset$Class, levels=c("No", "Yes"))

# add virus data [also obtained from the original publication]
virus <- read_csv("data/virus_data.csv", col_names=T) %>% 
  gather(key="barcode", value="counts", -X1) %>% 
  spread(key="X1", value="counts")

virus$barcode <- gsub(virus$barcode, pattern="\\.", replacement="-")

# make a new column indicating whether the cell barcodes have viral reads or not
full_dataset$Virus <- "No"
full_dataset[full_dataset$Barcode %in% virus$barcode, "Virus"] <- "Yes"


#### 2. Table, chi-squared test and barplot from figure 3a ####
table_microbes <- full_dataset %>% dplyr::select(Class, Disease) %>% table()
chisq_test(table_microbes, simulate.p.value = F)

bplot <- table_microbes %>% 
  prop.table(., margin=2) %>% 
  as_tibble()

ggbarplot(bplot, x="Disease", y="n", fill="Class", color="Class",
          palette=wes_palette("Zissou1",5)[c(2:3)])+
  theme_bw() +
  labs(title="Proportion of cells associated with bacteria", 
       fill="Bacteria", x="Disease group", y="Relative proportion") +
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
ggsave(filename = "output/bact_disease_association.pdf", device="pdf", height=4, width=4, useDingbats=F)


#### 3. Table, chi-squared tests and barplot from figure 3b ####

# plot bacteria per cell type
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

# for covid patients
table_distrib_COVID <- full_dataset %>% 
  dplyr::filter(Disease=="COVID19") %>% dplyr::select(Class, CellType) %>% table
chisq_test(table_distrib_COVID, simulate.p.value = F)
std_residuals(chisq_test(table_distrib_COVID, simulate.p.value = F)) %>% as_tibble()
chisq.posthoc.test(table_distrib_COVID,method = "BH",round = 96)

# for controls
table_distrib_contr <-  full_dataset %>% 
  dplyr::filter(Disease!="COVID19") %>% dplyr::select(Class, CellType) %>% table
chisq_test(table_distrib_contr, simulate.p.value = F)
std_residuals(chisq_test(table_distrib_contr, simulate.p.value = F)) %>% as_tibble()
chisq.posthoc.test(table_distrib_contr,method = "BH",round = 96)

# note that only enrichments in bacteria-associated cells [Dimension = "Yes" and Residuals > 0] are reported to facilitate interpretation

bplot1 <- table_distrib_COVID %>% 
  prop.table(., margin=1) %>% 
  as_tibble() %>% 
  mutate(Disease="COVID19")
bplot2 <- table_distrib_contr %>% 
  prop.table(., margin=1) %>% 
  as_tibble() %>% 
  mutate(Disease="Control")

bplot <- bind_rows(bplot1, bplot2)

ggbarplot(bplot, x="Class", y="n", fill="CellType") +#, palette=wes_palette("Darjeeling1", 18, "continuous")) +
  scale_fill_manual(values=mycolors) +
  theme_bw() +
  facet_grid(~ Disease) +
  labs(title="", 
       fill="Bacteria", x="Cells associated with bacteria", y="Relative proportion") +
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
ggsave(filename = "output/bact_cells_association.pdf", device="pdf", height=6, width=6, useDingbats=F)

# supplementary plots showing the standardized residual scores

stdres <- std_residuals(chisq_test(table_distrib_COVID, simulate.p.value = F))

pdf("output/std_residuals_COVID_cells.pdf", height=3, width=7)
ggcorrmod(t(stdres),title = "Std. residuals for cell types associated to bacteria (COVID-19)",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()

stdres <- std_residuals(chisq_test(table_distrib_contr, simulate.p.value = F))

pdf("output/std_residuals_control_cells.pdf", height=3, width=7)
ggcorrmod(t(stdres),title = "Std. residuals for cell types associated to bacteria (controls)",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()


#### 4. Table, chi-squared tests and plots for figure 3c (specific bacteria associations) ####

# filter cell types enriched in bacteria
whichbact <- full_dataset %>% 
  dplyr::filter(CellType %in% c("Md_macrophage", "Monocyte", "Neutrophil")) %>% 
  dplyr::filter(!is.na(species))

# do it at genus level so first load the taxonomy table
taxontable <- read_tsv("data/db_taxonomy.txt")

# join taxonomy to group at the genus level
whichbact <- whichbact %>% 
  left_join(taxontable, by="species") %>% 
  dplyr::select(-c(taxon_oid, genome_name, kingdom, phylum, class, order, family)) %>% 
  distinct

# chi-squre test and plot controls
tablebact_type <- whichbact %>% 
  dplyr::filter(Disease=="control") %>% 
  dplyr::select(genus, CellType) %>% table
chisq_test(tablebact_type, simulate.p.value = F)
std_residuals(chisq_test(tablebact_type, simulate.p.value = F)) %>% as_tibble()
chisq.posthoc.test(tablebact_type, method = "BH", round = 96)
stdres <- std_residuals(chisq_test(tablebact_type, simulate.p.value = T))
stdres <- stdres[apply(stdres,1,function(X){max(abs(X))>3}),]

pdf("output/specbact_cells_association_control.pdf", height=4, width=4)
ggcorrmod(stdres,title = "Controls",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()

# chi-squre test and plot COVID19
tablebact_type <- whichbact %>% 
  dplyr::filter(Disease=="COVID19") %>% 
  dplyr::select(genus, CellType) %>% table
chisq_test(tablebact_type, simulate.p.value = F)
std_residuals(chisq_test(tablebact_type, simulate.p.value = F)) %>% as_tibble()
chisq.posthoc.test(tablebact_type, method = "BH", round = 96)
stdres <- std_residuals(chisq_test(tablebact_type, simulate.p.value = T))
stdres <- stdres[apply(stdres,1,function(X){max(abs(X))>3}),]

pdf("output/specbact_cells_association_covid.pdf", height=4, width=4)
ggcorrmod(stdres,title = "COVID-19",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()


#### 5. Analysis on neutrophils (figures 3d, 3e) ####

# load neutrophils object [Seurat original object, subsetting cells with CellType=="neutrophil"]
neutrophils <- readRDS("data/neutrophils.rds")

# add metadata variable for the samples that have bacteria
neutrophils[["Bacteria"]] <- "No"
neutrophils@meta.data[intersect(rownames(neutrophils@meta.data), microbes_host$Barcode), "Bacteria"] <- "Yes"

# make contingency table of neutrophil subtypes with and without bacteria
tab_neutrophils <- table(neutrophils$Bacteria, neutrophils$CellType1)

# chisq test and plot
stdres <- rstatix::std_residuals(rstatix::chisq_test(tab_neutrophils, simulate.p.value = F))
chisq.posthoc.test(tab_neutrophils, method="BH", round = 16)
stdres <- stdres["Yes",,drop=F]
pdf("output/bact_cells_association_neutrophils.pdf", height=4, width=6)
ggcorrmod(stdres,title = "Neutrophils",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()


# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(neutrophils) <- "integrated"
neutrophils <- SetIdent(neutrophils, cells = colnames(neutrophils), value="CellType1")
neutrophil.markers <- FindAllMarkers(neutrophils, assay = "integrated", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- neutrophil.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# calculate average expression of the markers per each cluster, but separate cells with and without bacteria
cluster.averages <- AverageExpression(neutrophils, add.ident = "Bacteria", return.seurat = TRUE)
cluster.averages[["Bacteria"]] <- c("No", "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No")

pdf("output/heatmap_clusters_neutrophils_bact.pdf", height=7, width=5)
DoMultiBarHeatmap(cluster.averages, features=top10$gene, 
                  additional.group.by = "Bacteria", draw.lines = F, raster=F) + 
  scale_fill_gradientn(colors=c("#053061", "#F7F7F7", "#B2182B"))
dev.off()

neutrophil_markers <- neutrophil.markers %>% dplyr::filter(cluster=="C1_active_S100A12")
write.table(neutrophil_markers, "neutrophil_markers_bacteria.txt", sep="\t", quote=F, row.names=T, col.names=T)


#### 6. Analysis on md-macrophages and monocytes (figures 3d, 3f) ####

# load monocytes/macrophages object [Seurat original object, subsetting cells with CellType== Macrophages or monocytes]
monocytes_macrophages <- readRDS("data/monocytes_macrophages.rds")

## add metadata variable for the samples that have bacteria
monocytes_macrophages[["Bacteria"]] <- "No"
monocytes_macrophages@meta.data[intersect(rownames(monocytes_macrophages@meta.data), microbes_host$Barcode), "Bacteria"] <- "Yes"

# make contingency table of monocyte/macrophages subtypes with and without bacteria
tab_monocytes_macrophages <- table(monocytes_macrophages$Bacteria, monocytes_macrophages$CellType1)

# chisq test and plot (here, remove alveolar macrophages, as they were not significant in the previous analyses)
stdres <- rstatix::std_residuals(rstatix::chisq_test(tab_monocytes_macrophages, simulate.p.value = T))
chisq.posthoc.test(tab_monocytes_macrophages, method="BH", round = 16)

stdres <- stdres[,setdiff(colnames(stdres), c("Alveolar_Macrophage", "Alveolar_Macrophage_CSF1R"))]
stdres <- stdres["Yes",,drop=F]

pdf("output/bact_cells_association_monocytes_macrophages.pdf", height=4, width=6)
ggcorrmod(stdres,title = "Monocytes_macrophages",
          legend.title = "Std. residuals") + 
  theme(plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
dev.off()


# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(monocytes_macrophages) <- "integrated"
monocytes_macrophages <- SetIdent(monocytes_macrophages, cells = colnames(monocytes_macrophages), value="CellType1")
monocytes_macrophages.markers <- FindAllMarkers(monocytes_macrophages, assay = "integrated", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- monocytes_macrophages.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# calculate average expression of the markers per each cluster, but separate cells with and without bacteria
# here, for figure 3f, filter only enriched subtypes
cluster.averages <- AverageExpression(monocytes_macrophages, add.ident = "Bacteria", return.seurat = TRUE)
cluster.averages <- subset(cluster.averages, cells=c("Macrophage_CCL2_No", "Monocyte_IL1B_No", 
                                                     "Macrophage_CCL2_Yes",  "Monocyte_IL1B_Yes"))
cluster.averages[["Bacteria"]] <- c("No", "No", "Yes", "Yes")

pdf("output/heatmap_clusters_monomac_bact.pdf", height=7, width=5)
DoMultiBarHeatmap(cluster.averages, features=top10 %>% dplyr::filter(!grepl("Alveolar",cluster)) %>% pull(gene), 
                  additional.group.by = "Bacteria", draw.lines = F, raster=F) + 
  scale_fill_gradientn(colors=c("#053061", "#F7F7F7", "#B2182B"))
dev.off()

# besides plotting the differences in marker genes (not shown in the paper), 
# plot also the functional myeloid gene set (shown in the original publication)

genes <- c("IL1B", "IL6", "CCL2", "CCL3", "CCL4", "CCL7", "TNF", "CXCL2", "CXCL9", "CXCL10", "SOCS3", 
           "CD163", "MRC1", "MSR1", "CCL13", "CCL18", "IL10", "FOLR2", "STAB1", "CD163L1", "PLD4", "F13A1", "MERTK",
           "AXL", "CCL22", "GPNMB", "CHI3L1","B2M", 
           "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "HLA-DMB")

pdf("output/heatmap_clusters_monomac_bact2.pdf", height=7, width=5)
DoMultiBarHeatmap(cluster.averages, features=genes , raster=F,
                  additional.group.by = "Bacteria", draw.lines = F) + 
  scale_fill_gradientn(colors=c("#053061", "#F7F7F7", "#B2182B"))
dev.off()

monocyte_macrophage_markers <- monocytes_macrophages.markers %>% dplyr::filter(cluster %in% c("Macrophage_CCL2", "Monocyte_IL1B"))
write.table(monocyte_macrophage_markers, "monomac_markers_bacteria.txt", sep="\t", quote=F, row.names=T, col.names=T)


