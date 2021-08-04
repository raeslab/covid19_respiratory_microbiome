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
library(glmulti)
library(sjPlot)
library(lme4)
source("R/importance_glmulti.R")


#### 1. load data ####
load("data/covid19_genus.rda")

genus_data <- otu_table(covid19_genus)
genus_metadata <- sample_data(covid19_genus)

#### 2. alpha diversity - glmulti model selection ####

# calculate shannon diversity
shannon <- diversity(genus_data, index="shannon")

# prepare data for analyses and plots
plot_alphadiv <- tibble(samples=names(shannon), shannon, 
                        patient=genus_metadata$patient,
                        ICU_days = genus_metadata$days_icu_sampling,
                        Ward_days = genus_metadata$days_ward_sampling,
                        total_ICU = genus_metadata$length_icu,
                        Viral_load_scaled=log2(genus_metadata$total_sarscov2_viral_load+1),
                        CRP=genus_metadata$crp_deval, 
                        calpro=genus_metadata$s100a8,
                        mero_treatment=genus_metadata$meropenem_treatment, # refers to either meropenem or piperazillin/tazobactam, they were merged
                        ceftr_treatment=genus_metadata$ceftriaxone_treatment,
                        ventilated_previously=genus_metadata$ventilated_previously)

plot_alphadiv <- plot_alphadiv %>% 
  mutate(ICU_days = coalesce(ICU_days, total_ICU)) # for patients in ward (recently discharged), include the total days in ICU

## Use glmulti to select best model with selected features
alphadiv_simplified <- plot_alphadiv %>% 
  dplyr::select(shannon, patient, ICU_days, Viral_load_scaled, CRP, ceftr_treatment, mero_treatment,calpro, ventilated_previously) %>% 
  drop_na() %>% 
  mutate(calpro_scaled=log2(calpro)) %>% 
  mutate(crp_scaled=log2(CRP))

# create function for the GLMULTI algorithm
glmer.glmulti <- function (formula, data, random = "", family=gaussian(link = "log"), ...) {
  glmer(paste(deparse(formula), random), data = data, family=family,
        ...)
}

# run glmulti with the random patient effect
glmulti_alphadiv_patient <- glmulti(y=shannon~ICU_days+Viral_load_scaled +
                                      crp_scaled+mero_treatment+ceftr_treatment+ventilated_previously+calpro_scaled, 
                                    data=alphadiv_simplified, 
                                    level = 1,              
                                    method = "h",         
                                    crit="aicc",
                                    confsetsize = 128,       
                                    fitfunc = glmer.glmulti,
                                    random="+ (1|patient)",
                                    intercept=TRUE)

# same, without incuding the random effect
glm.glmulti <- function (formula, data, family=gaussian(link ="log"), ...) {
  glm(paste(deparse(formula)), data = data, family=family,
      ...)
}

glmulti_alphadiv_nop <- glmulti(y=shannon~ICU_days+Viral_load_scaled+
                                  crp_scaled+mero_treatment+ceftr_treatment+ventilated_previously+calpro_scaled,  
                                data=alphadiv_simplified, 
                                level = 1,              
                                method = "h",         
                                crit="aicc",
                                confsetsize = 128,       
                                fitfunc = glm.glmulti,
                                intercept=TRUE)

# lower AIC corresponds to the model including random effects, so we  use the top model from that algorithm
# plot variable importance (for fixed effects only)
plot(glmulti_alphadiv_patient, type="s")

# make nicer plot to save. First get importance values (I extracted the code from the plot.glmulti function)
importance_coefs <- importance_glmulti(glmulti_alphadiv_patient) %>% sort
importance_coefs <- importance_coefs %>% 
  as_tibble() %>% 
  mutate(Terms=names(importance_coefs)) 

importance_coefs <- importance_coefs %>% 
  mutate(Terms=factor(Terms, levels=(importance_coefs$Terms)))

# plot
ggbarplot(importance_coefs, x="Terms", y="value", color="white", 
          fill=wes_palette("Darjeeling1",5)[1]) +
  theme_bw()+
  rotate() +
  labs(title = "Model-averaged importance of terms", y="Importance") +
  # geom_hline(yintercept=0.8) +
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"))
ggsave(filename = "output_plots_2020/alpha_diversity_ALLMODELS_importance_coefficients.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)


# Choose the glmer best model including patient as random effect
glmermodel <- glmulti_alphadiv_patient@objects[[1]]


plot_model(glmermodel, show.values = TRUE, value.offset = .1, show.intercept = T) +
  theme_bw() +
  geom_hline(yintercept = 0, lwd=1, color="gray70")+
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Fixed effects estimates on Shannon diversity",
       subtitle = "Model = shannon ~ 1 + ICU_days + Viral_load_scaled + calpro_scaled + (1 | patient)\nFamily=gaussian(link='log')")
ggsave(filename = "output_plots_2020/alpha_diversity_BESTMODEL_coefficients.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)


p <- plot_model(glmermodel, type="pred", colors=get_palette("Blues",10)[c(4,6,8,10)],
                terms=c("ICU_days", "Viral_load_scaled [5,10,15,20]"), 
                pred.type="re", show.data=T)

p$layers[[1]]$aes_params$size <- 1.5

p + theme_bw() +
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Predicted values of Shannon diversity",
       y="Shannon diversity",
       subtitle = "Model = shannon ~ 1 + ICU_days + Viral_load_scaled + calpro_scaled + (1 | patient)\nFamily=gaussian(link='log')")
ggsave(filename = "output_plots_2020/alpha_diversity_BESTMODEL_predict_days_viralload.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)

p <- plot_model(glmermodel, type="pred", colors=get_palette("Reds",10)[c(4,6,8,10)],
                terms=c("ICU_days", "calpro_scaled [10,13,16,18]"), 
                pred.type="re", show.data=T)

p$layers[[1]]$aes_params$size <- 1.5

p + theme_bw() +
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +
  labs(title = "Predicted values of Shannon diversity",
       y="Shannon diversity",
       subtitle = "Model = shannon ~ 1 + ICU_days + Viral_load_scaled + calpro_scaled + (1 | patient)\nFamily=gaussian(link='log')")
ggsave(filename = "output_plots_2020/alpha_diversity_BESTMODEL_predict_days_calpro.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)



# compare to a null model using ANOVA
fullmodel <- glmer(shannon~ICU_days+Viral_load_scaled+calpro_scaled+(1|patient), 
                   data=alphadiv_simplified, family=gaussian(link="log"))
nullmodel <- glmer(shannon~(1|patient), data=alphadiv_simplified, family=gaussian(link="log"))
anova(fullmodel, nullmodel)


#### 3. alpha diversity - intra patient effect ####

## test further the effect of specific antibiotics, ventilation 
# select patients with samples before/after meropenem/piperacillin-tazobactam (for antibiotics in general or ceftriaxone there are no cases)
samples_meropenem <- c()
for(patientid in unique(genus_metadata$patient)){
  abs <- genus_metadata %>% as_tibble %>% dplyr::filter(patient==patientid) %>% 
    pull(meropenem_treatment)
  if(length(unique(na.omit(abs)))>1){
    samples_meropenem <- c(samples_meropenem, genus_metadata %>% as_tibble %>% 
                             dplyr::filter(patient==patientid) %>% 
                             dplyr::filter(meropenem_treatment=="0") %>% 
                             dplyr::filter(date==max(date)) %>% 
                             pull(sample))
    samples_meropenem <- c(samples_meropenem, genus_metadata %>% as_tibble %>% 
                             dplyr::filter(patient==patientid) %>% 
                             dplyr::filter(meropenem_treatment=="1") %>% 
                             dplyr::filter(date==min(date)) %>% 
                             pull(sample))
  }
}

paired_data <- plot_alphadiv %>% 
  dplyr::filter(samples %in% samples_meropenem) %>% 
  dplyr::select(shannon, patient, ICU_days, mero_treatment) %>% 
  arrange(mero_treatment)
paired_data$mero_treatment <- factor(paired_data$mero_treatment)

ggpaired(paired_data, x="mero_treatment", y="shannon", fill="mero_treatment", id="patient",
         add="jitter",line.color = "gray", line.size = 0.4,
         palette=inlmisc::GetColors(5, scheme = "sunset")[c(2,1)],
         xlab="Meropenem administered", 
         ylab="Shannon diversity", 
         legend.title="Meropenem",
         title="Intra-patient effect of meropenem/piperacillin-tazobactam") + 
  stat_compare_means(method="wilcox", paired = T, comparisons=list(c("0", "1"))) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "gray97"), 
        panel.grid.minor = element_line(colour = "gray97")) + 
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(size = 14, 
                                  face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output_plots_2020/alpha_diversity_INTRAPATIENT_meropenem.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)





# select patients with samples before/after MV 
samples_MV <- c()
for(patientid in unique(genus_metadata$patient)){
  mv <- genus_metadata %>% as_tibble %>% dplyr::filter(patient==patientid) %>% 
    pull(ventilated_previously)
  if(length(unique(na.omit(mv)))>1){
    samples_MV <- c(samples_MV, genus_metadata %>% as_tibble %>% 
                      dplyr::filter(patient==patientid) %>% 
                      dplyr::filter(ventilated_previously==0) %>% 
                      dplyr::filter(date==max(date)) %>% 
                      pull(sample))
    samples_MV <- c(samples_MV, genus_metadata %>% as_tibble %>% 
                      dplyr::filter(patient==patientid) %>% 
                      dplyr::filter(ventilated_previously==1) %>% 
                      dplyr::filter(date==min(date)) %>% 
                      pull(sample))
  }
}

paired_data <- plot_alphadiv %>% 
  dplyr::filter(samples %in% samples_MV) %>% 
  dplyr::select(shannon, patient, ICU_days, ventilated_previously) %>% 
  arrange(ventilated_previously)
paired_data$ventilated_previously <- factor(paired_data$ventilated_previously)

ggpaired(paired_data, x="ventilated_previously", y="shannon", fill="ventilated_previously", id="patient",
         add="jitter",line.color = "gray", line.size = 0.4,
         palette=inlmisc::GetColors(5, scheme = "sunset")[c(2,1)],
         xlab="Mechanical ventilation", 
         ylab="Shannon diversity", title = "Intra-patient effect of mechanical ventilation",
         legend.title="Mechanical ventilation") + 
  stat_compare_means(method="wilcox", paired = T, comparisons=list(c("0", "1"))) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "gray97"), 
        panel.grid.minor = element_line(colour = "gray97")) + 
  theme(axis.title = element_text(face = "bold"), 
        plot.title = element_text(size = 14, 
                                  face = "bold"), legend.title = element_text(face = "bold"))
ggsave(filename = "output_plots_2020/alpha_diversity_INTRAPATIENT_ventilation.pdf", device="pdf", height=6, width=6, useDingbats=FALSE)

