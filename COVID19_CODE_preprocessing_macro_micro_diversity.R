######################################
#############PRE-PROCESSING###########
######################################

library(dada2)
library(phyloseq)

path <- "/Users/u0126662/Desktop/CONTAGIOUS_16S" # CHANGE ME to the directory containing your demultiplexed R1 and R2 fastq files
fns <- list.files(path)
fns

fastqs <- fns[grepl(".fq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
#### make sure that R1 is for forward read and R2 for reverse!!!!

fnFs <- fastqs[grepl(".1.fq", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl(".2.fq", fastqs)] # Just the reverse read files
# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, ".1.fq"), `[`, 1) ## check if is 1 or 2!


# Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
###############
pdf("plotQualityProfile.pdf", onefile=T)
plotQualityProfile(fnFs[1:10]) ## remove 20 plus 10
plotQualityProfile(fnRs[1:10])  ## remove 20 plus 10
dev.off()
#### remove primers
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter #### important remove primers!!! and remove low quality regions

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(130,200), ##150
                     trimLeft=c(30, 30),
                     maxN=0, maxEE=c(2,2), truncQ=11, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) #
head(out)

###############
pdf("plotQualityProfile.filt.pdf", onefile=T)
plotQualityProfile(filtFs[1:15])
plotQualityProfile(filtRs[1:15])
dev.off()

set.seed(12345)

# Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

pdf("plotErrors_F.pdf", onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf("plotErrors_R.pdf", onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


derepRs <- derepFastq(filtRs, verbose=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct sequence table and remove chimeras

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)


# Assign taxonomy
taxHS <- assignTaxonomy(seqtab.nochim, "rdp_train_set_18.fa.gz", multithread=TRUE)
taxHS <- addSpecies(taxHS, "rdp_species_assignment_18.fa.gz")
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(head(taxHS))
unname(tail(taxHS))

#ASV_names <- paste0("SVs",1:ncol(seqtab.nochim))
#colnames(seqtab.nochim) <- ASV_names
#rownames(taxHS) <- ASV_names


# Write to disk
write.table(track, file = "APR_track.tsv", quote=FALSE)
write.table(seqtab.nochim, file = "APR_sequence_table_SV.tsv", quote=FALSE)
write.table(taxHS, file = "APR_taxa_SV.tsv", quote=FALSE)

number_individuals <- as.data.frame(rowSums(seqtab.nochim))
second_min_value <- sort(number_individuals[,1])[2]
Samples_IDs <- rownames(seqtab.nochim)
Sample_data_table <- as.data.frame(read.table("Sample_data_table.txt", sep="\t", header=TRUE))
rownames(Sample_data_table) <- Samples_IDs

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(Sample_data_table), 
               tax_table(taxHS))

######################################
############DECONTAMINATION###########
######################################

#remove contamination
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove contamination from phyloseq object
contamination <- as.integer(as.logical(contamdf.prev$contaminant))
contaminants <- c("1")
replacements <- c("3")
contamination <- replace(contamination, contamination %in% contaminants, replacements)
noncontaminants <- c("0")
replacements <- c("1")
contamination <- replace(contamination, contamination %in% noncontaminants, replacements)
contaminants <- c("3")
replacements <- c("0")
contamination <- replace(contamination, contamination %in% contaminants, replacements)
contamination <- as.logical(as.integer(contamination))
ps.nocontam = prune_taxa(contamination, ps)

Sample_data_table2 <- as.data.frame(read.table("Sample_data_table2.txt", sep="\t", header=TRUE))
rownames(Sample_data_table2) <- Samples_IDs
Nonblanks<-as.character(Sample_data_table2$Logic)
Nonblanks <- as.logical(as.integer(Nonblanks))
ps.nocontam.noblanks = prune_samples(Nonblanks, ps.nocontam)

saveRDS(ps.nocontam.noblanks, "COVID_contamremoved_Dec2020.rdata")
COVID_contamremoved <- readRDS("COVID_contamremoved_Dec2020.rdata")

##############################################
############STRAIN/SPECIES ANALYSES###########
##############################################

library("phyloseq")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("tidyverse")
library("reshape2")
library("compositions")
library("vegan")
library("rstatix")
library("GLDEX")
library("ggrepel")
library("tibble")
library("DECIPHER")
library("Biostrings")
library("biomod2")

setwd("/Users/u0126662/Desktop/CONTAGIOUS_16S")
load("metadata_table_16S_samples.RData")
load("coding_table_final.RData")
COVID_contamremoved <-readRDS("COVID_contamremoved_Dec2020.rdata")
metadata_table_16S_samples <- as.data.frame(metadata_table_16S_samples)
rownames(metadata_table_16S_samples)<-metadata_table_16S_samples$sample
TAXA <- as.data.frame(tax_table(COVID_contamremoved))
META <- as.data.frame(read.table("Sample_data_table3.txt", sep="\t", header=TRUE))
IDs <- META$Official_ID
GENUS_ps_object <- tax_glom(COVID_contamremoved, taxrank="Genus") 
GENUS <- otu_table(GENUS_ps_object) 
TAXA_GENUS <- as.data.frame(tax_table(GENUS_ps_object))
colnames(GENUS) <- TAXA_GENUS$Genus 
rownames(GENUS) <- IDs
GENUS <- GENUS[, which(colSums(GENUS) != 0)]
OTU <- otu_table(COVID_contamremoved) 
rownames(OTU) <- IDs
OTU <- as.data.frame(OTU)
Repeats <- c("COVID_ICU_002_D.1","COVID_ICU_002_D.2","COVID_ICU_010_C.1","COVID_ICU_010_C.2","COVID_ICU_010_C.3","COVID_ICU_012_A.1","COVID_ICU_012_A.2","COVID_ICU_014_A.1","COVID_ICU_014_A.2","COVID_ICU_014_A.3","COVID_ICU_014_Y.1","COVID_ICU_014_Y.2","COVID_ICU_014_Y.3","COVID_ICU_023_A.1","COVID_ICU_023_A.2","COVID_ICU_023_A.3","COVID_ICU_037_Z.1","COVID_ICU_037_Z.2","COVID_ICU_042_Z.1","COVID_ICU_042_Z.2","COVID_ICU_042_Z.3","COVID_ICU_044_A.1","COVID_ICU_044_A.2","COVID_ICU_046_Z.1","COVID_ICU_046_Z.2","COVID_ICU_046_Z.3","COVID_W_049_C.1","COVID_W_049_C.2","COVID_W_049_C.3","COVID_ICU_050_D.1","COVID_ICU_050_D.2","COVID_ICU_050_D.3","COVID_ICU_055_E1.1","COVID_ICU_055_E1.2","COVID_ICU_056_E1.1","COVID_ICU_056_E1.2","COVID_ICU_056_E1.3","COVID_ICU_060_A.1","COVID_ICU_060_A.2","COVID_ICU_060_A.3")

OTU_COVID_ICU_002_D <- OTU[c("COVID_ICU_002_D.1","COVID_ICU_002_D.2"),]
OTU_COVID_ICU_010_C <- OTU[c("COVID_ICU_010_C.1","COVID_ICU_010_C.3","COVID_ICU_010_C.3"),]
OTU_COVID_ICU_012_A <- OTU[c("COVID_ICU_012_A.1","COVID_ICU_012_A.2"),]
OTU_COVID_ICU_014_A <- OTU[c("COVID_ICU_014_A.1","COVID_ICU_014_A.2","COVID_ICU_014_A.3"),]
OTU_COVID_ICU_014_Y <- OTU[c("COVID_ICU_014_Y.1","COVID_ICU_014_Y.2","COVID_ICU_014_Y.3"),]
OTU_COVID_ICU_023_A <- OTU[c("COVID_ICU_023_A.1","COVID_ICU_023_A.2","COVID_ICU_023_A.3"),]
OTU_COVID_ICU_037_Z <- OTU[c("COVID_ICU_037_Z.1","COVID_ICU_037_Z.2"),]
OTU_COVID_ICU_042_Z <- OTU[c("COVID_ICU_042_Z.1","COVID_ICU_042_Z.2","COVID_ICU_042_Z.3"),]
OTU_COVID_ICU_044_A <- OTU[c("COVID_ICU_044_A.1","COVID_ICU_044_A.2"),]
OTU_COVID_ICU_046_Z <- OTU[c("COVID_ICU_046_Z.1","COVID_ICU_046_Z.2","COVID_ICU_046_Z.3"),]
OTU_COVID_W_049_C <- OTU[c("COVID_W_049_C.1","COVID_W_049_C.2","COVID_W_049_C.3"),]
OTU_COVID_ICU_050_D <- OTU[c("COVID_ICU_050_D.1","COVID_ICU_050_D.2","COVID_ICU_050_D.3"),]
OTU_COVID_ICU_055_E1 <- OTU[c("COVID_ICU_055_E1.1","COVID_ICU_055_E1.2"),]
OTU_COVID_ICU_056_E1 <- OTU[c("COVID_ICU_056_E1.1","COVID_ICU_056_E1.2","COVID_ICU_056_E1.3"),]
OTU_COVID_ICU_060_A <- OTU[c("COVID_ICU_060_A.1","COVID_ICU_060_A.2","COVID_ICU_060_A.3"),]

COVID_ICU_002_D <- as.data.frame(colMeans(OTU_COVID_ICU_002_D))
colnames(COVID_ICU_002_D) <- "COVID_ICU_002_D"
COVID_ICU_002_D <- t(COVID_ICU_002_D)
COVID_ICU_010_C <- as.data.frame(colMeans(OTU_COVID_ICU_010_C))
colnames(COVID_ICU_010_C) <- "COVID_ICU_010_C"
COVID_ICU_010_C <- t(COVID_ICU_010_C)
COVID_ICU_012_A <- as.data.frame(colMeans(OTU_COVID_ICU_012_A))
colnames(COVID_ICU_012_A) <- "COVID_ICU_012_A"
COVID_ICU_012_A <- t(COVID_ICU_012_A)
COVID_ICU_014_A <- as.data.frame(colMeans(OTU_COVID_ICU_014_A))
colnames(COVID_ICU_014_A) <- "COVID_ICU_014_A"
COVID_ICU_014_A <- t(COVID_ICU_014_A)
COVID_ICU_014_Y <- as.data.frame(colMeans(OTU_COVID_ICU_014_Y))
colnames(COVID_ICU_014_Y) <- "COVID_ICU_014_Y"
COVID_ICU_014_Y <- t(COVID_ICU_014_Y)
COVID_ICU_023_A <- as.data.frame(colMeans(OTU_COVID_ICU_023_A))
colnames(COVID_ICU_023_A) <- "COVID_ICU_023_A"
COVID_ICU_023_A <- t(COVID_ICU_023_A)
COVID_ICU_037_Z <- as.data.frame(colMeans(OTU_COVID_ICU_037_Z))
colnames(COVID_ICU_037_Z) <- "COVID_ICU_037_Z"
COVID_ICU_037_Z <- t(COVID_ICU_037_Z)
COVID_ICU_042_Z <- as.data.frame(colMeans(OTU_COVID_ICU_042_Z))
colnames(COVID_ICU_042_Z) <- "COVID_ICU_042_Z"
COVID_ICU_042_Z <- t(COVID_ICU_042_Z)
COVID_ICU_044_A <- as.data.frame(colMeans(OTU_COVID_ICU_044_A))
colnames(COVID_ICU_044_A) <- "COVID_ICU_044_A"
COVID_ICU_044_A <- t(COVID_ICU_044_A)
COVID_ICU_046_Z <- as.data.frame(colMeans(OTU_COVID_ICU_046_Z))
colnames(COVID_ICU_046_Z) <- "COVID_ICU_046_Z"
COVID_ICU_046_Z <- t(COVID_ICU_046_Z)
COVID_W_049_C <- as.data.frame(colMeans(OTU_COVID_W_049_C))
colnames(COVID_W_049_C) <- "COVID_W_049_C"
COVID_W_049_C <- t(COVID_W_049_C)
COVID_ICU_050_D <- as.data.frame(colMeans(OTU_COVID_ICU_050_D))
colnames(COVID_ICU_050_D) <- "COVID_ICU_050_D"
COVID_ICU_050_D <- t(COVID_ICU_050_D)
COVID_ICU_055_E1 <- as.data.frame(colMeans(OTU_COVID_ICU_055_E1))
colnames(COVID_ICU_055_E1) <- "COVID_ICU_055_E1"
COVID_ICU_055_E1 <- t(COVID_ICU_055_E1)
COVID_ICU_056_E1 <- as.data.frame(colMeans(OTU_COVID_ICU_056_E1))
colnames(COVID_ICU_056_E1) <- "COVID_ICU_056_E1"
COVID_ICU_056_E1 <- t(COVID_ICU_056_E1)
COVID_ICU_060_A <- as.data.frame(colMeans(OTU_COVID_ICU_060_A))
colnames(COVID_ICU_060_A) <- "COVID_ICU_060_A"
COVID_ICU_060_A <- t(COVID_ICU_060_A)

OTU_no_repeats = OTU[!row.names(OTU)%in%Repeats,]
OTU_fixed <- as.data.frame(rbind(OTU_no_repeats,COVID_ICU_002_D,COVID_ICU_010_C,COVID_ICU_012_A,COVID_ICU_014_A,COVID_ICU_014_Y,COVID_ICU_023_A,COVID_ICU_037_Z,COVID_ICU_042_Z,COVID_ICU_044_A,COVID_ICU_046_Z,COVID_W_049_C,COVID_ICU_050_D,COVID_ICU_055_E1,COVID_ICU_056_E1,COVID_ICU_060_A))

target <- rownames(OTU_fixed)
Metadata <- metadata_table_16S_samples[match(target, metadata_table_16S_samples$sample),]
coding_table_final <- as.data.frame(coding_table_final)


#making OTUs otu of ASV data
OTU_fixed <- as.data.frame(OTU_fixed)
OTU_fixed_binary <- BinaryTransformation(OTU_fixed, 1)
OTU_fixed_binary <- as.data.frame(OTU_fixed_binary)
rownames(OTU_fixed_binary) <- rownames(OTU_fixed)

nproc <- 4
seqtab <- OTU_fixed_binary
seqtab2 <- OTU_fixed
asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)


## Find clusters of ASVs to form the new OTUs
set.seed(12)
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03, # use `cutoff = 0.03` for a 97% OTU 
  processors = nproc)

numbers <- c(1:1194)
TAXA$species_full <- paste(TAXA$Genus,TAXA$Species,numbers)


## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  add_column(sequence = asv_sequences)
clusters2 <- clusters %>%
  add_column(species_name = TAXA$species_full)

write.csv(clusters2,file = "clusters.csv")

count_number_ASVs_per_OTU <- seqtab %>%
  # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
  t %>%
  as_tibble(rownames = "sequence") %>%
  # add the cluster information
  left_join(clusters, by = "sequence") %>%
  # merge ASVs in the same cluster, summing abundances within samples
  group_by(cluster) %>%
  summarize_at(vars(-sequence), sum) %>%
  # Set new taxa names to OTU<cluster #> 
  mutate(cluster = paste0("OTU", cluster)) %>%
  # convert back to a matrix in the original orientation
  column_to_rownames("cluster") %>%
  as("matrix") %>%
  t

OTU_real_table <- seqtab2 %>%
  # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
  t %>%
  as_tibble(rownames = "sequence") %>%
  # add the cluster information
  left_join(clusters, by = "sequence") %>%
  # merge ASVs in the same cluster, summing abundances within samples
  group_by(cluster) %>%
  summarize_at(vars(-sequence), sum) %>%
  # Set new taxa names to OTU<cluster #> 
  mutate(cluster = paste0("OTU", cluster)) %>%
  # convert back to a matrix in the original orientation
  column_to_rownames("cluster") %>%
  as("matrix") %>%
  t

write.csv(count_number_ASVs_per_OTU,file = "97_ASV_counts_per_OTU_per_sample.csv")
write.csv(OTU_real_table,file = "97_OTU_real_table.csv")


Microdiversity_dataframe <-as.data.frame(read.table("97_ASV_counts_per_OTU_per_sample.csv", sep=",", row.names = 1, header=TRUE))
Macrodiversity_dataframe <-as.data.frame(read.table("97_OTU_real_table.csv", sep=",", row.names = 1, header=TRUE))

column_numbers <- c() # 5 is minimum 
for (i in 1:nrow(Microdiversity_dataframe)){
	row <- Microdiversity_dataframe[i,]
	row_no_zeros <- row[,colSums(row) > 0]
	cols <- ncol(row_no_zeros)
	column_numbers <- rbind(column_numbers, cols)
}


microdiversity_averages <- c()
for (i in 1:nrow(Microdiversity_dataframe)){
	row <- Microdiversity_dataframe[i,]
	row_no_zeros <- row[,colSums(row) > 0]
	id <- rownames(row)
	if (rowSums(row)[[1]]== 1) {
		microdiversity_averages <- rbind(microdiversity_averages,cbind(id,"NA"))
	}
	else {
		averages <- c()
		for (j in 1:1000) {
			number <- 3+j
			set.seed(number)
			sampled_5 <- row_no_zeros[, sample(ncol(row_no_zeros), 5, replace = TRUE)]
			average <- rowSums(sampled_5)/5
			averages <- rbind(averages, average)
		}
		mean_of_mean <- mean(averages)
		microdiversity_averages <- rbind(microdiversity_averages,cbind(id,mean_of_mean))
	}
}

microdiversity_averages_w_metadata <- as.data.frame(cbind(microdiversity_averages, Metadata))
microdiversity_averages_w_metadata <- as.data.frame(microdiversity_averages_w_metadata)
microdiversity_averages_w_metadata$mean_of_mean <- as.numeric(as.character(microdiversity_averages_w_metadata$mean_of_mean))
write.table(microdiversity_averages_w_metadata, file='97_ASV_microdiversity_averages.tsv', quote=FALSE,col.names=NA, sep='\t')



my_comparisons <- list(c("1","2"),c("1","3"),c("1","4"),c("1","5"),c("1","6"),c("1","7"), c("2","3"),c("2","4"),c("2","5"),c("2","6"),c("2","7"),c("3","4"),c("3","5"),c("3","6"),c("3","7"),c("4","5"),c("4","6"),c("4","7"),c("5","6"),c("5","7"),c("6","7"))
Microdiversity_averages_by_patient_status_plot <- ggboxplot(microdiversity_averages_w_metadata, x="deval_sup_specoxyg", y = "mean_of_mean", alpha = 0.7,lwd=0.25, add = "jitter", shape = 21, fill = "deval_sup_specoxyg")+stat_compare_means(comparisons = my_comparisons)+ylab("Average Microdiversity")+xlab("")+theme(legend.position="right")

Microdiversity_by_oxysupport <- c()
Supports <- c("1","2","3","4","5","6","7")
for (i in 1:length(Supports)) {
	Support <- Supports[i]
	Microdiversity <- microdiversity_averages_w_metadata[microdiversity_averages_w_metadata$deval_sup_specoxyg==Support,2]
	Microdiversity <- Microdiversity[!is.na(Microdiversity)]
	for (j in 1:100) {
		sampled_5 <- sample(Microdiversity, 5, replace = TRUE)
		average <- mean(sampled_5)
		Microdiversity_by_oxysupport <- rbind(Microdiversity_by_oxysupport,cbind(Support,average))
	}
}
Microdiversity_by_oxysupport <- as.data.frame(Microdiversity_by_oxysupport)
Microdiversity_by_oxysupport$average <- as.numeric(as.character(Microdiversity_by_oxysupport$average))
my_comparisons <- list(c("1","2"),c("1","3"),c("1","4"),c("1","5"),c("1","6"),c("1","7"), c("2","3"),c("2","4"),c("2","5"),c("2","6"),c("2","7"),c("3","4"),c("3","5"),c("3","6"),c("3","7"),c("4","5"),c("4","6"),c("4","7"),c("5","6"),c("5","7"),c("6","7"))
my_comparisons2 <- list(c("1","2"),c("2","3"),c("2","4"),c("2","5"),c("2","6"),c("2","7"))
Microdiversity_averages_by_oxysupport_sub5i <- ggboxplot(Microdiversity_by_oxysupport, x="Support", y = "average", alpha = 0.7,lwd=0.25, add = "jitter", shape = 21, fill = "Support")+stat_compare_means(comparisons = my_comparisons2)+ylab("Average Strain-Level Diversity")+xlab("")+theme(legend.position="none")
Microdiversity_averages_by_oxysupport_sub5 <- ggboxplot(Microdiversity_by_oxysupport, x="Support", y = "average", alpha = 0.7,lwd=0.25, fill = "Support")+ scale_fill_manual(values =c("#EE332E","#50675A","#51A967","#F1AB28","#F79226","#BD914D","#51BDD4"))+ylab("Average Strain-Level Diversity")+xlab("")+theme(legend.position="right")

Richness<-specnumber(Macrodiversity_dataframe) #calculates richness
Shannons_H<- diversity(Macrodiversity_dataframe, index = "shannon") #calculates shannon index
Simpson<- diversity(Macrodiversity_dataframe, index = "simpson") #calculates simpson index
InvSimpson<- diversity(Macrodiversity_dataframe, index = "invsimpson") #calculates inverse simpson index
Peilous_J <- Shannons_H/log(Richness) #calculates peilou's J
alpha_diversity <- cbind(Richness, Shannons_H, Simpson, InvSimpson, Peilous_J)
write.table(alpha_diversity, file = "OTU_Alpha_diveristy_stats.tsv", sep = "\t", row.names = TRUE, col.names=NA)
Macrodiversity_w_metadata <- as.data.frame(cbind(Richness,Shannons_H,Simpson, Metadata))
write.table(Macrodiversity_w_metadata, file='97_OTU_macrodiversity_averages.tsv', quote=FALSE,col.names=NA, sep='\t')

Macrodiversity_by_oxygen_support_plot <- ggboxplot(Macrodiversity_w_metadata, x="deval_sup_specoxyg", y = "Richness", alpha = 0.7,lwd=0.25, add = "jitter", shape = 21,fill = "deval_sup_specoxyg")+ylab("Richness")+xlab("")

Macrodiversity_by_oxysupport <- c()
Supports <- c("1","2","3","4","5","6","7")
for (i in 1:length(Supports)) {
	Support <- Supports[i]
	Macrodiversity <- Macrodiversity_w_metadata[Macrodiversity_w_metadata$deval_sup_specoxyg==Support,1]
	Macrodiversity <-Macrodiversity[!is.na(Macrodiversity)]
	for (j in 1:100) {
		sampled_5 <- sample(Macrodiversity, 5, replace = TRUE)
		average <- mean(sampled_5)
		Macrodiversity_by_oxysupport <- rbind(Macrodiversity_by_oxysupport,cbind(Support,average))
	}
}
Macrodiversity_by_oxysupport <- as.data.frame(Macrodiversity_by_oxysupport)
Macrodiversity_by_oxysupport$average <- as.numeric(as.character(Macrodiversity_by_oxysupport$average))
my_comparisons <- list(c("1","2"),c("1","3"),c("1","4"),c("1","5"),c("1","6"),c("1","7"), c("2","3"),c("2","4"),c("2","5"),c("2","6"),c("2","7"),c("3","4"),c("3","5"),c("3","6"),c("3","7"),c("4","5"),c("4","6"),c("4","7"),c("5","6"),c("5","7"),c("6","7"))
Macrodiversity_averages_by_oxysupport_sub5i <- ggboxplot(Macrodiversity_by_oxysupport, x="Support", y = "average", alpha = 0.7,lwd=0.25, add = "jitter", shape = 21, fill = "Support")+stat_compare_means(comparisons = my_comparisons)+ylab("Average Macrodiversity")+xlab("")+theme(legend.position="none")
Macrodiversity_averages_by_oxysupport_sub5 <- ggboxplot(Macrodiversity_by_oxysupport, x="Support", y = "average", alpha = 0.7,lwd=0.25, fill = "Support")+ scale_fill_manual(values =c("#EE332E","#50675A","#51A967","#F1AB28","#F79226","#BD914D","#51BDD4"))+ylab("Average Species-Level Diversity")+xlab("")+theme(legend.position="none")

Macrodiversity_averages_of_averages <- aggregate(Macrodiversity_by_oxysupport[, 2], list(Macrodiversity_by_oxysupport$Support), mean)
Microdiversity_averages_of_averages <- aggregate(Microdiversity_by_oxysupport[, 2], list(Microdiversity_by_oxysupport$Support), mean)
macro_micro_bind <- cbind(Macrodiversity_averages_of_averages, Microdiversity_averages_of_averages[,2])
colnames(macro_micro_bind) <- c("Support","Macro","Micro")
macro_micro_bind <- as.data.frame(macro_micro_bind)
macro_micro_bind$Support <- as.numeric(as.character(macro_micro_bind$Support))
#macro_micro_bind_no_7 <- head(macro_micro_bind,-1)
#macro_micro_bind_only_7 <- tail(macro_micro_bind,1)


scatter_plot <- ggplot() + 
  geom_point(data=macro_micro_bind, aes(x=Macro, y=Micro))+
  geom_smooth(data=macro_micro_bind, aes(x=Macro, y=Micro),method=lm, se=FALSE, color = "#525456")+
  geom_point(data=macro_micro_bind, aes(x=Macro, y=Micro))+
  geom_label_repel(data=macro_micro_bind,aes(x=Macro, y=Micro, label=Support))+
  ylab("Average Strain-Level Diversity")+xlab("Average Species-Level Diversity")+
  theme_classic()+
  theme(aspect.ratio=1)


ggscatter(macro_micro_bind, x = "Macro", y = "Micro",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)
  
  
grid::current.viewport()

plot_grid(Macrodiversity_averages_by_oxysupport_sub5, Microdiversity_averages_by_oxysupport_sub5, scatter_plot, ncol=3, rel_widths = c(0.75,1,1), align = 'hv')

