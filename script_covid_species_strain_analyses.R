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

COVID_contamremoved <-readRDS("data/COVID_contamremoved.rdata")mple
TAXA <- as.data.frame(tax_table(COVID_contamremoved))
META <- as.data.frame(sample_data(Sample_data_table))
OTU <- otu_table(COVID_contamremoved) 

#making OTUs otu of ASV data
OTU <- as.data.frame(OTU)
OTU_binary <- BinaryTransformation(OTU, 1)
OTU_binary <- as.data.frame(OTU)
rownames(OTU_binary) <- rownames(OTU)

nproc <- 4
seqtab <- OTU_binary
seqtab2 <- OTU
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

