# Beta diversity analyses
# Upper respiratory microbiome in COVID-19 patients

# load packages and libraries
library(tidyverse)
library(phyloseq)
library(CoDaSeq)
library(vegan)
library(ggpubr)
library(DESeq2)
library(mixOmics)
library(wesanderson)
library(colortools)
library(ggiraphExtra)
library(ggiraph)
library(rstatix)


#### 1. load data ####
load("data/covid19_genus.rda")

genus_data <- otu_table(covid19_genus)
genus_metadata <- sample_data(covid19_genus)

#### 2. beta diversity -- preprocess data ####

min.reads <- 10000 # read number threshold
min.prop <- 0.001 # taxa abundance threshold (minimum abundance in any sample)
cutoff <- 0.1  # taxa prevalence threshold

genus_filt <- codaSeq.filter(genus_data, min.reads=min.reads, min.occurrence=cutoff, 
                             min.prop=min.prop, samples.by.row=T) 

# function to estimate zeros based on the minimum relative abundance
estimate0.min <- function(matrix.test) { 
  print(paste("You have",ncol(matrix.test),"samples. If this is not correct, transpose matrix!"))
  matrix.test.p <- t(t(matrix.test)/rowSums(t(matrix.test)))
  matrix.test.p[is.nan(matrix.test.p)] <- 0
  samplesums <- colSums(matrix.test)
  
  matrix.f.n0 <- matrix.test
  for (i in 1:nrow(matrix.f.n0)) {
    min <- min(matrix.test.p[i,][matrix.test.p[i,] > 0])
    for (j in 1:ncol(matrix.f.n0 )) {
      if (matrix.f.n0 [i,j] == 0)
        matrix.f.n0 [i,j] <- min*samplesums[j]
    }
  }
  return(matrix.f.n0)
}

# replace zeros using the function above
genus_filt_n0 <- estimate0.min(genus_filt)

# perform CLR transformation on the filtered, zero-replaced data
genus_filt_n0_clr <- codaSeq.clr(genus_filt_n0, samples.by.row=F)


#### 5. Univariate RDA ####
set.seed(123)

# prepare metadata matrix 
genus_metadata <- as_tibble(genus_metadata) %>% 
  mutate(samples=rownames(genus_metadata))

metadata_curated <- genus_metadata  %>% 
  as.data.frame
rownames(metadata_curated) <- metadata_curated$sample
metadata_curated <- metadata_curated[,setdiff(colnames(metadata_curated),c("sample", "samples"))] # remove sample names from the metadata variable list
metadata_curated <- metadata_curated[colnames(genus_filt_n0_clr),] # make sure the samples are sorted

# filter metadata matrix
cols_toremove <- names(which(apply(metadata_curated, 2, function(x){length(unique(na.omit(x)))})==1))
metadata_curated <- metadata_curated[,setdiff(colnames(metadata_curated), cols_toremove)]
cols_toremove <- names(which(apply(metadata_curated, 2, function(x){length(which(is.na(x)))/length(x)})>0.25))
metadata_curated <- metadata_curated[,setdiff(colnames(metadata_curated), cols_toremove)]
# remove also dates
cols_toremove <- colnames(metadata_curated)[grep(colnames(metadata_curated), pattern="date")]
metadata_curated <- metadata_curated[,setdiff(colnames(metadata_curated), cols_toremove)]
# remove technical variables that refer only to whether a test was performed or not (not the result of the test)
cols_toremove <- c("hc_bas","resppanel_swab_bas", "broncho_deval")
metadata_curated <- metadata_curated[,setdiff(colnames(metadata_curated), cols_toremove)]
# days in ICU variable -> for patients in ward count the whole stay, for patients in ICU count current days
metadata_curated[is.na(metadata_curated$days_icu_sampling),"days_icu_sampling"] <- metadata_curated[is.na(metadata_curated$days_icu_sampling),"length_icu"]


# dbRDA, univariate, loop to test all covariates included in the metadata_curated data frame
all <- c()
for(i in 1:ncol(metadata_curated)){
  cap <- capscale(t(genus_filt_n0_clr) ~ metadata_curated[,i], distance="euclidean", na.action=na.omit) 
  av <- cap %>% 
    anova.cca %>% 
    as.data.frame
  pval <- av$`Pr(>F)`[1]
  Fstat <- av$`F`[1]
  r2 <- RsquareAdj(cap)[[1]]
  adjr2 <- RsquareAdj(cap)[[2]]
  all <- rbind(all, cbind(Fstat,r2,adjr2,pval))
}

# add corresponding covariate names and select only those significant ones after multiple testing correction
rownames(all) <- colnames(metadata_curated)
all <- as.data.frame(cbind(all, padj=p.adjust(all[,"pval"], method="BH")), stringsAsFactors=F)
significant <- all[all[,"padj"]<0.05,]
significant <- na.omit(significant[order(significant$adjr2, decreasing = T),])
print(paste0(nrow(significant), " variables can individually explain part of the (active) microbiome variation in these samples."))

write.table(all, "output_plots_2020/table_dbRDA.txt", col.names = T,row.names=T, quote=F, sep="\t")

# Stepwise multivariate RDA, with only significant variables
metadata_curated_filtered <- metadata_curated[,rownames(significant), drop=F] # select only significant variables
metadata_curated_filtered <- na.omit(metadata_curated_filtered) # samples with NAs need to be removed from the entire matrix for multivariate analysis

# remove samples with NAs also from the CLR object
reduced_genus_filt_n0_clr <- genus_filt_n0_clr[,rownames(metadata_curated_filtered)]

# set limits for the stepwise analysis, from null model (fit0) to a full model with the 20 selected variables
fit0 <- capscale(t(reduced_genus_filt_n0_clr) ~ 1, data=metadata_curated_filtered, distance="euclidean")
fit1 <- capscale(t(reduced_genus_filt_n0_clr) ~ ., data=metadata_curated_filtered, distance="euclidean")

# run stepwise analysis between the two limits established
fin_dbrda <- ordiR2step(fit0, fit1, data=metadata_curated_filtered, Pin = 0.05, R2scope = TRUE, 
                        pstep = 100, perm.max = 200, trace=F)

# Print call and anova
fin_dbrda$call
fin_dbrda$anova

# test effect of oxygen support w/o confounders:
# individual effect
capscale(t(reduced_genus_filt_n0_clr) ~ deval_sup_specoxyg, 
         data=metadata_curated_filtered, distance="euclidean") %>% anova
# effect independent of patient
capscale(t(reduced_genus_filt_n0_clr) ~ deval_sup_specoxyg + 
           Condition(patient), data=metadata_curated_filtered, distance="euclidean") %>% anova
# effect independent of patient AND antibiotics: meropenem/piperazillin-tazobactam previous treatment (meropenem_treatment) 
# and ongoind ceftriaxone treatment (specab1_deval_1)
capscale(t(reduced_genus_filt_n0_clr) ~ deval_sup_specoxyg + 
           Condition(patient + meropenem_treatment + specab1_deval_1), 
         data=metadata_curated_filtered, distance="euclidean") %>% anova 


# reorganize univariate and multivariate results for plotting
significant <- as.tibble(cbind(significant, covariates=rownames(significant)))
anova_table <- as.tibble(cbind(fin_dbrda$anova, covariates=gsub(rownames(fin_dbrda$anova), 
                                                                pattern="\\+ ",
                                                                replacement="")))

significant <- significant %>% 
  full_join(anova_table, by="covariates")

significant[is.na(significant$R2.adj),"R2.adj"] <- max(significant$R2.adj,na.rm=T)
significant <- significant %>% 
  filter(covariates!="<All variables>")

# gather, define sorting order, and plot
significantL <- gather(significant, key = "Class", value="EffectSize", adjr2, R2.adj)

significantL[significantL$Class=="adjr2", "Class"] <- "Individual"
significantL[significantL$Class=="R2.adj", "Class"] <- "Cumulative, non-redundant"

load("data/coding_table_final.RData") # table including covariate explanations and categories assigned
colnames(coding_table_final)[5] <- "Question_category"
significantL <- significantL %>% 
  left_join(coding_table_final, by=c("covariates"="Variable"))
significantL[significantL$covariates=="days_icu_sampling","Question_category"] <- "COVID"
positions <- rev.default(significant$covariates[order(significant$R2.adj, significant$AIC,
                                                      -significant$adjr2,decreasing=F)])

# plot all effect sizes
ggplot(significantL, aes(x = covariates, y = EffectSize)) +
  geom_bar(aes(color = NULL, fill = Question_category, alpha=Class), 
           stat = "identity", position = position_dodge()) + 
  ggpubr::rotate() + 
  theme_bw() +
  scale_x_discrete(limits = positions) +
  scale_alpha_manual(values=c(0.6,1)) +
  geom_vline(xintercept = nrow(significant)-(nrow(anova_table)-1.5))+
  labs(x="Covariates", y="Effect size") +
  scale_fill_manual(values=wes_palette("Darjeeling1", 6, type = "continuous"))
ggsave(filename = "output_plots_2020/covariates.pdf", device="pdf", width=7, height=6)


# for the significant covariates in the multivariate analyses, plot the RDA and colour by these variables
covs <- anova_table$covariates[1:(length(anova_table$covariates)-1)]
otus <- otu_table(t(reduced_genus_filt_n0_clr), taxa_are_rows = F)
envdata <- sample_data(cbind(metadata_curated_filtered, group=metadata_curated[rownames(metadata_curated_filtered), "group"]))
physeq_object <- phyloseq(otus, envdata)
GP.ord <- ordinate(physeq_object, "RDA", "euclidean", as.formula(paste0("~ ", paste(covs, collapse=" + "))))


for(cov in covs){
  print(plot_ordination(physeq_object, GP.ord, type="samples", color=paste(cov),
                        shape="group", title = paste0("RDA - colored by ", cov)) +
          theme_bw() + geom_point(size=2) +
          theme(legend.text = element_text(size = 6),
                legend.key = element_rect(fill =NA), legend.background = element_rect(fill = NA),
                legend.position="right"))
  ggsave(filename = paste0("output_plots_2020/RDA_", cov, ".pdf"), device="pdf", width=7, height=6, useDingbats=FALSE)
}


## PERMANOVA tests 
genus_metadata <- sample_data(covid19_genus)
genus_metadata <- as_tibble(genus_metadata) %>% 
  mutate(samples=rownames(genus_metadata))


genus_metadata_filt <- genus_metadata %>% 
  dplyr::filter(samples %in% colnames(genus_filt_n0_clr)) %>% 
  mutate(ventilation=ifelse(deval_sup_specoxyg %in% c("4", "5", "6", "7"),"invasive", "noninvasive")) %>% 
  dplyr::select(sample, patient, moment, group, deval_cs_pat, ventilation, deval_sup_specoxyg, 
                samples, meropenem_treatment, ventilated_previously, specab1_deval_1) %>% 
  drop_na()

genus_filt_n0_clr <- genus_filt_n0_clr[,genus_metadata_filt$samples]

# permanova for differences across groups
permanova_oxygen_allcategories <- adonis(t(genus_filt_n0_clr) ~ deval_sup_specoxyg, data = genus_metadata_filt, method="euclidean")
permanova_ventilation_ongoing <- adonis(t(genus_filt_n0_clr) ~ ventilation, data = genus_metadata_filt, method="euclidean")
permanova_ventilated_prev <- adonis(t(genus_filt_n0_clr) ~ ventilated_previously, data = genus_metadata_filt, method="euclidean")
permanova_meropenem_prev <- adonis(t(genus_filt_n0_clr) ~ meropenem_treatment, data = genus_metadata_filt, method="euclidean")
permanova_ceftriaxone_current <- adonis(t(genus_filt_n0_clr) ~ specab1_deval_1, data = genus_metadata_filt, method="euclidean")

