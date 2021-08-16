# Exploratory plots
# Upper respiratory microbiome in COVID-19 patients


# load packages and libraries
library(tidyverse)
library(lubridate)
library(phyloseq)
library(vegan)
library(ggpubr)
library(wesanderson)
library(colortools)
library(ggiraphExtra)
library(ggiraph)

#### 1. Load data ####

# metadata 
metadata_table_16S_samples <- read_tsv("data/metadata_table_16S_samples.txt") # data object per samples (only info from day of sampling)
metadata_table_16S_patients <- read_tsv("data/metadata_table_16S_patients.txt") # data object per patients (info from all days of the study)

metadata_table_16S_patients <- metadata_table_16S_patients %>% 
  filter(!is.na(subjid)) # remove blanks/controls from metadata table

# if date of hospital admission is posterior to date of ICU admission, replace data of hospital admission
# this means the patient has been transferred from another hospital
metadata_table_16S_samples <- metadata_table_16S_samples %>% 
  filter(!is.na(patient)) %>% 
  mutate(date_hospital_admission=if_else(date_hospital_admission>=date_icu_admission, date_icu_admission, date_hospital_admission)) 



#### 2. summary statistics, for use in table 1 ####
sum_table <- metadata_table_16S_samples %>% 
  mutate(length_hospi=date_hospital_last-date_hospital_admission) %>% 
  dplyr::select(patient, length_hospi,age_bas,sex_bas, height_bas, weight_bas, bmi_bas, dia_bas, smoking_bas, length_icu) %>% 
  distinct() %>% 
  summarize(Number_patients=n(),Age=mean(age_bas, na.rm=T),Age_range=paste(range(age_bas, na.rm=T), collapse="-"),
            BMI=mean(bmi_bas, na.rm=T), BMI_range=paste(range(bmi_bas, na.rm=T), collapse="-"),
            Diabetic=length(which(dia_bas==1)), Diabetic_percentage= 100*length(which(dia_bas==1))/n(),
            length_ICU=mean(length_icu, na.rm=T), length_ICU_range=paste(range(length_icu, na.rm=T), collapse="-"),
            length_hospital=round(as.numeric(mean(length_hospi, na.rm=T)),1), length_hospital_range=paste(round(as.numeric(range(length_hospi, na.rm=T)),1), collapse="-"),
            Female_sex=length(which(sex_bas==2)))
write_tsv(sum_table, "output_plots_2020/table_patients.txt", col_names = T)



#### 3. Plot for figure 1a ####
toplot <- metadata_table_16S_samples %>% 
  dplyr::select(patient,date_hospital_admission, date_icu_admission,date_icu_last, 
                date_hospital_last) %>% 
  distinct() %>% 
  gather(key="Event", value="Date", date_icu_admission,date_icu_last, 
         date_hospital_admission,date_hospital_last)

##add new rows for patients still hospitalized on July 1st (to limit the x axis)
newrows <- tibble(patient=c("COVID_ICU_010", "COVID_ICU_053"),  Event=c("ICU","ICU"), Date=as.Date("2020-07-01"))
toplot <- bind_rows(toplot,newrows)
# and remove the rows corresponding to the last day of data
toplot <- toplot %>% 
  dplyr::filter(!(patient %in% c("COVID_ICU_010", "COVID_ICU_053") & Event %in% c("date_hospital_last", "date_icu_last")))

##add new rows for patients who when to ICU more than once
newrows <- tibble(patient="COVID_W_031", Event=c("Ward","ICU"), Date=as.Date(c("2020-04-20", "2020-04-23")))
toplot <- bind_rows(toplot,newrows)
newrows <- tibble(patient="COVID_W_049", Event=c("Ward","ICU","Ward"), Date=as.Date(c("2020-04-24", "2020-04-25", "2020-04-26")))
toplot <- bind_rows(toplot,newrows)
newrows <- tibble(patient="COVID_ICU_059", Event=c("Ward","ICU"), Date=as.Date(c("2020-05-19", "2020-04-23")))
toplot <- bind_rows(toplot,newrows)
newrows <- tibble(patient="COVID_ICU_041", Event=c("Ward","ICU"), Date=as.Date(c("2020-04-09", "2020-04-16")))
toplot <- bind_rows(toplot,newrows)

toplot$Event <- gsub(toplot$Event, pattern="date_", replacement="")
toplot$Event <- gsub(toplot$Event, pattern="hospital_admission", replacement="Ward")
toplot$Event <- gsub(toplot$Event, pattern="icu_last", replacement="Ward")
toplot$Event <- gsub(toplot$Event, pattern="icu_admission", replacement="ICU")
toplot$Event <- gsub(toplot$Event, pattern="hospital", replacement="Hospital")
toplot$Event <- factor(toplot$Event, levels=c("Ward", "ICU", "Hospital_last"))

toplot <- toplot %>% dplyr::arrange(patient,Date,Event)

toplot <- toplot %>% 
  left_join(metadata_table_16S_samples %>% dplyr::select(patient, date), by="patient")

toplot <- toplot %>% dplyr::filter(patient!="COVID_W_061") # remove for plotting because no data on dates is available

# plot timelines for the patients
ggplot(toplot, aes(x=Date, y=patient, group=patient, color=Event)) + 
  scale_color_manual(values=wes_palette(n=5, name="Zissou1")[c(4,2,5)]) +
  geom_line(size=1) + 
  geom_point(shape=15, size=0.7) +
  geom_point(aes(x=date, y=patient), color="black", size=1, shape=4) +
  theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "16S samples - CONTAGIOUS cohort", y="Patient",
       x = "Date", colour = "Event")
ggsave(filename = "output_plots_2020/samples_16S_full.pdf", device = "pdf", width = 5, height=7)




#### 4. Plot timelines of antibiotic administration (for figure S2) ####

toplot <- metadata_table_16S_patients %>% 
  dplyr::select(subjid, deval_date, deval_icu_ongoing, specab1_deval_1, specab1_deval_meropip, ab_deval, deval_sup_specoxyg) %>% 
  distinct() %>% 
  mutate(oxygen_type=ifelse(deval_sup_specoxyg %in% c(4,5,6,7), yes = "Invasive ventilation", no=NA)) %>% 
  mutate(ICU=ifelse(deval_icu_ongoing=="1", yes="ICU", no=NA)) %>% 
  dplyr::filter(deval_date<="2020-07-01") 

# ANTIBIOTIC (any) + VENTILATION (not used in the figure S2, but for additional reference)
ggplot(toplot, aes(x=deval_date, y=subjid, group=subjid, color=ICU)) +
  geom_line(size=4)+
  geom_line(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=oxygen_type), size=3)+
  geom_point(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=ab_deval), size=1, shape=15)+
  scale_color_manual(values=c(wes_palette(n=5, name="Zissou1")[c(4,2)], "gray", "gray30"), na.translate=F) +
  geom_point(inherit.aes = F, data=metadata_table_16S_samples, aes(x=date, y=patient), shape=4, size=1, color="black") +
  theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Antibiotic administration", y="Patient",
       x = "Date", colour = "Antibiotic")
ggsave(filename = "output_plots_2020/samples_antibiotic_ventilation.pdf", device = "pdf", width = 6, height=7)


# ANTIBIOTIC CEFTRIAXONE + VENTILATION
ggplot(toplot, aes(x=deval_date, y=subjid, group=subjid, color=ICU)) +
  geom_line(size=4)+
  geom_line(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=oxygen_type), size=3)+
  geom_point(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=specab1_deval_1), size=1, shape=15)+
  scale_color_manual(values=c(wes_palette(n=5, name="Zissou1")[c(4,2)], "gray", "gray30"), na.translate=F) +
  geom_point(inherit.aes = F, data=metadata_table_16S_samples, aes(x=date, y=patient), shape=4, size=1) +
  theme_bw()+ 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Ceftriaxone administration", y="Patient",
       x = "Date", colour = "Ceftriaxone")
ggsave(filename = "output_plots_2020/samples_ceftriaxone_ventilation.pdf", device = "pdf", width = 6, height=7)


# ANTIBIOTIC MEROPENEM + VENTILATION
ggplot(toplot, aes(x=deval_date, y=subjid, group=subjid, color=ICU)) +
  geom_line(size=4)+
  geom_line(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=oxygen_type), size=3)+
  geom_point(inherit.aes = F, data=toplot, aes(x=deval_date, y=subjid, group=subjid, color=specab1_deval_meropip), size=1, shape=15)+
  scale_color_manual(values=c(wes_palette(n=5, name="Zissou1")[c(4,2)], "gray", "gray30"), na.translate=F) +
  geom_point(inherit.aes = F, data=metadata_table_16S_samples, aes(x=date, y=patient), shape=4, size=1) +
  theme_bw()+ 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Meropenem/piperacillin-tazobactam administration", y="Patient",
       x = "Date", colour = "Meropenem/Pip-tazo")
ggsave(filename = "output_plots_2020/samples_meropenem_ventilation.pdf", device = "pdf", width = 6, height=7)



#### 5. Abundance table: load and convert to phyloseq OTU table ####
ps <- readRDS("data/COVID_contamremoved.rdata")
abundance_table <- otu_table(ps) 

abundance_table <- abundance_table %>% 
  as.data.frame() %>% 
  as_tibble() %>%
  mutate(ID=rownames(abundance_table)) %>% 
  gather(key="ASV", value="counts", -ID) 

sample_ids <- read_tsv("sample_ids.txt", col_names=F)
abundance_table <- abundance_table %>% left_join(sample_ids, by=c("ID" = "X1"))
colnames(abundance_table)[4] <- "Sample"

abundance_table$Sample <- gsub(abundance_table$Sample, pattern="\\.[[:digit:]]", replacement="")

# pool technical replicates data
abundance_table <- abundance_table %>%
  left_join(metadata_table_16S_samples, by=c("Sample"="sample")) %>%
  dplyr::select(Sample, ASV, counts, patient, moment) %>%
  group_by(Sample, patient, moment, ASV) %>%
  summarize(counts=sum(counts)) %>%
  ungroup %>% 
  dplyr::filter(!is.na(patient)) %>% 
  dplyr::select(Sample, ASV, counts) %>%
  spread(key = ASV, value = counts) %>% 
  as.data.frame() 

rownames(abundance_table) <- abundance_table$Sample
abundance_table <- abundance_table[,-1]

# prepare metadata for phyloseq object
metadata_reduced <- metadata_table_16S_samples %>% 
  dplyr::filter(!is.na(patient))

metadata_reduced_df <- metadata_reduced %>% 
  as.data.frame()
rownames(metadata_reduced_df) <- metadata_reduced_df$sample

# assign sample_data and otu_table
samples <- sample_data(metadata_reduced_df)
ASV <- otu_table(as.matrix(abundance_table), taxa_are_rows = F)

# get taxonomy
taxonomy <- tax_table(ps)

# make phyloseq object
covid19 <- phyloseq(ASV, taxonomy, samples)
save(covid19, file = "data/covid19.rda")

# group at genus level and filter samples with >= 10000 reads
covid19_genus <- tax_glom(covid19, "Genus")
covid19_genus <- prune_samples(sample_sums(covid19_genus)>=10000, covid19_genus)
save(covid19_genus, file = "data/covid19_genus.rda")




#### 6. Plot samples ####

# first extract the genus-level data (perform all manipulations outside of the phyloseq object)
genus_data <- otu_table(covid19_genus)
genus_metadata <- sample_data(covid19_genus)
genus_tax <- tax_table(covid19_genus)
genus_tax <- as.data.frame(covid19_genus@tax_table@.Data) %>% 
  as_tibble() %>% 
  mutate(ASV=rownames(genus_tax))

# gather data and calculate relative proportions
genus_long <- as_tibble(as.data.frame(genus_data)) %>%
  mutate(samples=rownames(genus_data)) %>% 
  gather(key="ASV", value="counts", -samples) %>% 
  left_join(genus_tax, by="ASV") %>% 
  dplyr::select(samples, Genus, counts)

genus_long_relative <- genus_long %>% 
  group_by(samples) %>%
  mutate(countT= sum(counts)) %>%
  group_by(Genus, .add=TRUE) %>%
  mutate(counts=counts/countT) %>% 
  ungroup %>% 
  dplyr::select(-countT) %>% 
  mutate(Genus=as.character(Genus))

# choose top 15 per sample
top15 <- genus_long_relative %>% 
  group_by(Genus) %>%
  summarize(counts = sum(counts, na.rm=T)) %>%
  arrange(desc(counts)) %>%
  pull(Genus) %>%
  head(.,15)

simplifiedData <- genus_long_relative 
simplifiedData[!simplifiedData$Genus %in% top15,"Genus"] <- "Other"
simplifiedData$Genus <- factor(simplifiedData$Genus, levels=c(top15, "Other"))

simplifiedData <- simplifiedData %>%
  group_by(Genus, samples) %>%
  summarize(counts=sum(counts, na.rm=T)) %>% 
  left_join(metadata_reduced, by=c("samples" = "sample")) %>% 
  ungroup

# divide the sampling points in: ICU admission, stay and discharge/ward
simplifiedData$moment_simplified <- "ICU_Stay"
simplifiedData[simplifiedData$moment=="Admission_Date_sampling","moment_simplified"] <- "ICU_Admission"
simplifiedData[simplifiedData$moment=="Discharge_Date_sampling","moment_simplified"] <- "ICU_Discharge/Ward"
simplifiedData[simplifiedData$group=="Ward","moment_simplified"] <- "ICU_Discharge/Ward"

simplifiedData$moment_simplified <- factor(simplifiedData$moment_simplified, levels=c("ICU_Admission", "ICU_Stay", "ICU_Discharge/Ward"))

# plot
ggbarplot(simplifiedData, x = "samples", y="counts", color="black", fill="Genus",
          legend="right", 
          legend.title="Top genera", main="Relative counts per genus",
          font.main = c(14,"bold", "black"), font.x = c(12, "bold"), 
          font.y=c(12,"bold")) + 
  theme_bw() +
  rotate_x_text() + 
  scale_fill_manual(values=c(wheel("skyblue3", num = 15),"gray")) +
  facet_grid(~ moment_simplified, scales = "free_x", space='free') + 
  labs(x = "Sample", y = "Relative proportion") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) 
ggsave(filename = "output_plots_2020/relative_counts.pdf", device="pdf", width=8, height=4)



