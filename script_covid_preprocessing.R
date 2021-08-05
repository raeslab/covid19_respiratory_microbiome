######################################
#############PRE-PROCESSING###########
######################################

library(dada2)
library(phyloseq)

path <- "/path/to/reads" 
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
Sample_data_table <- as.data.frame(read.table("Metadata_table.txt", sep="\t", header=TRUE))
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

Metadata <- as.data.frame(read.table("Metadata_table.txt", sep="\t", header=TRUE))
rownames(Metadata) <- Samples_IDs
Nonblanks<-as.character(Metadata$Logic)
Nonblanks <- as.logical(as.integer(Nonblanks))
ps.nocontam.noblanks = prune_samples(Nonblanks, ps.nocontam)

saveRDS(ps.nocontam.noblanks, "data/COVID_contamremoved.rdata")

