# Alpha diversity analyses
# Upper respiratory microbiome in COVID-19 patients

# load packages and libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(wesanderson)
library(colortools)
library(ggiraphExtra)
library(ggiraph)
library(DESeq2)
library(CoDaSeq)


#### 1. load data ####
load("data/covid19_genus.rda")

genus_data <- otu_table(covid19_genus)
genus_metadata <- sample_data(covid19_genus)


#### 2. Prepare data ####

min.reads <- 10000 # read number threshold
min.prop <- 0.001 # taxa abundance threshold (minimum abundance in any sample)
cutoff <- 0.1  # taxa prevalence threshold

genus_filt <- codaSeq.filter(genus_data, min.reads=min.reads, min.occurrence=cutoff, 
                             min.prop=min.prop, samples.by.row=T) 

genus_metadata <- as_tibble(genus_metadata) %>% 
  mutate(samples=rownames(genus_metadata)) %>% 
  as.data.frame
rownames(genus_metadata) <- genus_metadata$samples

#### 3. Oxygen support, simplified (ventilation vs no ventilation, current) ####
# generate design matrix
data_deseq2 <- genus_filt
design <- genus_metadata[colnames(data_deseq2),c("patient", "moment", "group", "deval_sup_specoxyg", "specab1_deval_1", "meropenem_treatment", "ventilated_previously")]
design$oxygen <- "Non_invasive"
design[design$deval_sup_specoxyg %in% c("4","5","6", "7"),"oxygen"] <- "Invasive"
design$oxygen <- factor(design$oxygen, levels=c("Non_invasive", "Invasive"))
design <- na.omit(design)

# retrieve data
data_deseq2 <- data_deseq2[,rownames(design)]

# calculate differential expression using deseq2
dds <- DESeqDataSetFromMatrix(countData = data_deseq2,
                              colData = design,
                              design= ~ oxygen )
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
result_deseq <- results(dds)

# filter significant results
result_signif_table <- result_deseq[!is.na(result_deseq$padj) & result_deseq$padj<0.05,]
result_signif_table <- result_signif_table[order(result_signif_table$padj),]
result_signif <- rownames(result_signif_table)

# format and plot hits (p-values<0.05)
counts_mat <- varianceStabilizingTransformation(as.matrix(data_deseq2))
counts_mat <- cbind(ASV=rownames(counts_mat), counts_mat) %>% as.tibble %>%
  dplyr::filter(ASV %in% result_signif)

counts_matlong <- gather(counts_mat, key = "sample", value="counts", -ASV) %>%
  left_join(metadata_table_16S_samples, by="sample") %>%
  left_join(genus_tax, by="ASV") %>% 
  dplyr::select(patient, ASV, Genus, sample, counts, group, moment,deval_sup_specoxyg) %>% 
  mutate(oxygen=ifelse(deval_sup_specoxyg %in% c("4", "5", "6", "7"),"Invasive", "Non-invasive"))

counts_matlong$ASV <- factor(counts_matlong$ASV, levels=result_signif)
counts_matlong$Genus <- factor(counts_matlong$Genus, levels = counts_matlong %>% arrange(ASV) %>% pull(Genus) %>% unique)

counts_matlong$oxygen <- factor(counts_matlong$oxygen, levels=c("Non-invasive", "Invasive"))
counts_matlong$counts <- as.numeric(counts_matlong$counts)
boxes <- ggboxplot(counts_matlong, x="oxygen", y="counts", fill="oxygen", add="jitter",add.params = list(size=0.8),
                   title="Significant differences among oxygen support types", palette=wes_palette("Darjeeling1", 2, "continuous")[c(2,1)]) +
  theme_bw() +
  geom_line(aes(group=patient),alpha = 0.5, colour = "darkgrey", size=0.5) +
  labs(subtitle = "Test: Likelihood ratio test; Design: ~ oxygen_support_type ; Reduced model: ~ 1") +
  rotate_x_text()
facet(boxes, facet.by = "Genus", scales="free_y", nrow=5)
ggsave(filename = "output_plots_2020/taxa_oxygen_support_invasive.pdf", device="pdf", height=6, width=10, useDingbats=F)



#### 4. Oxygen support, simplified (ventilation vs no ventilation, current) - controlling for antibiotics significant in the dbRDA ####
# generate design matrix
data_deseq2 <- genus_filt
design <- genus_metadata[colnames(data_deseq2),c("patient", "moment", "group", "deval_sup_specoxyg", "specab1_deval_1", "meropenem_treatment", "ventilated_previously")]
design$oxygen <- "Non_invasive"
design[design$deval_sup_specoxyg %in% c("4","5","6", "7"),"oxygen"] <- "Invasive"
design$oxygen <- factor(design$oxygen, levels=c("Non_invasive", "Invasive"))
design <- na.omit(design)

# retrieve data
data_deseq2 <- data_deseq2[,rownames(design)]

# calculate differential expression using deseq2
dds <- DESeqDataSetFromMatrix(countData = data_deseq2,
                              colData = design,
                              design= ~ meropenem_treatment + specab1_deval_1 + oxygen )
dds <- DESeq(dds, test = "LRT", reduced = ~ meropenem_treatment + specab1_deval_1)
result_deseq <- results(dds)

# filter significant results
result_signif_table <- result_deseq[!is.na(result_deseq$padj) & result_deseq$padj<0.05,]
result_signif_table <- result_signif_table[order(result_signif_table$padj),]
result_signif <- rownames(result_signif_table)

# format and plot hits (p-values<0.05)
counts_mat <- varianceStabilizingTransformation(as.matrix(data_deseq2))#getVarianceStabilizedData(dds)
counts_mat <- cbind(ASV=rownames(counts_mat), counts_mat) %>% as.tibble %>%
  dplyr::filter(ASV %in% result_signif)

counts_matlong <- gather(counts_mat, key = "sample", value="counts", -ASV) %>%
  left_join(metadata_table_16S_samples, by="sample") %>%
  left_join(genus_tax, by="ASV") %>% 
  dplyr::select(patient, ASV, Genus, sample, counts, group, moment,deval_sup_specoxyg) %>% 
  mutate(oxygen=ifelse(deval_sup_specoxyg %in% c("4", "5", "6", "7"),"Invasive", "Non-invasive"))

counts_matlong$ASV <- factor(counts_matlong$ASV, levels=result_signif)
counts_matlong$Genus <- factor(counts_matlong$Genus, levels = counts_matlong %>% arrange(ASV) %>% pull(Genus) %>% unique)

counts_matlong$oxygen <- factor(counts_matlong$oxygen, levels=c("Non-invasive", "Invasive"))
counts_matlong$counts <- as.numeric(counts_matlong$counts)
boxes <- ggboxplot(counts_matlong, x="oxygen", y="counts", fill="oxygen", add="jitter",add.params = list(size=0.8),
                   title="Significant differences among oxygen support types", palette=wes_palette("Darjeeling1", 2, "continuous")[c(2,1)]) +
  theme_bw() +
  geom_line(aes(group=patient),alpha = 0.5, colour = "darkgrey", size=0.5) +
  labs(subtitle = "Test: Likelihood ratio test; Design: ~ meropenem_treatment + specab1_deval_1 + oxygen_support_type;\n Reduced model: ~ meropenem_treatment + specab1_deval_1") +
  rotate_x_text()
facet(boxes, facet.by = "Genus", scales="free_y", nrow=5)
ggsave(filename = "output_plots_2020/taxa_oxygen_support_invasive_ANTIBIOTICS.pdf", device="pdf", height=6, width=8, useDingbats=F)




#### 5. Viral load ####
# generate design matrix
data_deseq2 <- genus_filt
design <- genus_metadata[colnames(data_deseq2),c("patient", "moment", "group", "total_sarscov2_viral_load")]
design <- na.omit(design)
design$viral_load_scaled <- log2(design$total_sarscov2_viral_load+1)

# retrieve data
data_deseq2 <- data_deseq2[,rownames(design)]

# calculate differential expression using deseq2
dds <- DESeqDataSetFromMatrix(countData = data_deseq2,
                              colData = design,
                              design= ~ viral_load_scaled)
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
result_deseq <- results(dds)

# filter significant results
result_signif <- rownames(result_deseq[!is.na(result_deseq$padj) & result_deseq$padj<0.05,])


# plot hits (p-values<0.05)
counts_mat <- getVarianceStabilizedData(dds)
counts_mat <- cbind(ASV=rownames(counts_mat), counts_mat) %>% as.tibble %>%
  dplyr::filter(ASV %in% result_signif)

counts_matlong <- gather(counts_mat, key = "sample", value="counts", -ASV) %>%
  left_join(metadata_table_16S_samples, by="sample") %>%
  left_join(genus_tax, by="ASV") %>% 
  dplyr::select(patient, ASV, Genus, sample, counts, group, moment,total_sarscov2_viral_load) %>% 
  mutate(viral_load_scaled=log2(total_sarscov2_viral_load)+1)

counts_matlong$ASV <- factor(counts_matlong$ASV, levels=result_signif)
counts_matlong$Genus <- factor(counts_matlong$Genus, levels = counts_matlong %>% arrange(ASV) %>% pull(Genus) %>% unique)

counts_matlong$counts <- as.numeric(counts_matlong$counts)
points <- ggscatter(counts_matlong, x="viral_load_scaled", y="counts", add="reg.line",
                    conf.int = T,
                    title="Significant associations with scaled viral load") +
  theme_bw() +
  labs(subtitle = "Test: Likelihood ratio test; Design: ~ viral load (scaled) ; Reduced model: ~ 1")
facet(points, facet.by = "Genus", scales="free_y", nrow=1)
ggsave(filename = "output_plots_2020/taxa_viral_load.pdf", device="pdf", height=6, width=6)


#### 5. Viral load ####
# generate design matrix
data_deseq2 <- genus_filt
design <- genus_metadata[colnames(data_deseq2),c("patient", "moment", "group", "deval_sup_specoxyg", "specab1_deval_1", "meropenem_treatment", "ventilated_previously", "total_sarscov2_viral_load")]
design$oxygen <- "Non_invasive"
design[design$deval_sup_specoxyg %in% c("4","5","6", "7"),"oxygen"] <- "Invasive"
design$oxygen <- factor(design$oxygen, levels=c("Non_invasive", "Invasive"))
design <- na.omit(design)
design$viral_load_scaled <- log2(design$total_sarscov2_viral_load+1)

design <- design[!is.na(design$viral_load_scaled),]

# retrieve data
data_deseq2 <- data_deseq2[,rownames(design)]

# calculate differential expression using deseq2
dds <- DESeqDataSetFromMatrix(countData = data_deseq2,
                              colData = design,
                              design= ~ oxygen +  viral_load_scaled)
dds <- DESeq(dds, test = "LRT", reduced = ~ oxygen )
result_deseq <- results(dds)

# filter significant resutls
result_signif <- rownames(result_deseq[!is.na(result_deseq$padj) & result_deseq$padj<0.05,])

