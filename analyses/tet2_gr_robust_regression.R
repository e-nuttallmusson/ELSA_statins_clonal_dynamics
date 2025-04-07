library(dplyr)
library(readxl)
library(olsrr)
library(jtools)
library(broom.mixed)
library(lmtest)
library(Metrics)
library(caret)
library(ggplot2)
library(glmtoolbox)
library(truncnorm)
library(broom)
library(data.table)
library(mice)
library(MASS)
library(broom)
library(repmod)

#setwd

setwd("N:/My Documents/ELSA_data/FBC_results")

#read in data
fbc_data <- read.csv("fbc_clean_20250123.csv")

#convert to numeric
fbc_data$WBC_num <- as.numeric(fbc_data$WBC_10.9.L)
fbc_data$lymph_num <- as.numeric(fbc_data$LYA_10.9.L)

#calculate myeloid count

fbc_data$myeloid <- fbc_data$WBC_num - fbc_data$lymph_num

#calculate myeloid:lymphoid ratio

fbc_data$m_l_ratio <- fbc_data$myeloid/fbc_data$lymph_num


#simplify to just keep idauniq, wave and m_l_ratio

fbc_data_simp <- subset(fbc_data, select = c(ID, wave, m_l_ratio))


#check for duplications in wave 6
library(dplyr)
fbc_data_simp6 <- fbc_data_simp %>% dplyr::filter(wave == 6)
dups <- duplicated(fbc_data_simp6$ID)

table(dups)
#dups
#FALSE  TRUE 
#6456     2 

which(dups == TRUE)
#[1] 1395 5571
#> View(fbc_data_simp6)
#> fbc_data_simp6$ID[1395]
#[1] 108023
#> fbc_data__simp6$ID[5571]
#[1] 106847

#duplicated wave 6 ids with NA values:
fbc_data_simp[14086,]
fbc_data_simp[14302,]

#drop these
fbc_data_simp_dedup <- fbc_data_simp %>% filter(!row_number() %in% c(14086, 14302))

#reshape data long to wide
fbc_data_simp_w <- reshape(fbc_data_simp_dedup, idvar = "ID", v.names = "m_l_ratio", timevar = "wave", direction = "wide")


#calculate growth rate


setwd("N:/My Documents/ELSA_data/upgrade/tet2_vaf_progression_linear_regression/")

#read in repeat chip mutation data with wave 8/9 comorbidities and wave 6 medications

chip_repeats_comorbs_drugs <- read_xlsx("chip_repeats_all_comorbs_wave6_drug.xlsx")

#add variable of number of mutations
chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% group_by(idauniq) %>% mutate(count = n())

#select tet2 repeats
tet2_chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% dplyr::filter(grepl("TET2", Gene.refGene))

#remove those with haem disorder
tet2_chip_repeats_comorbs_drugs <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(blood_dis == 0)
tet2_chip_repeats_comorbs_drugs <- as.data.frame(tet2_chip_repeats_comorbs_drugs)

#create variable indicating which elsa_wave control sample was (labelled with vaf 0.01)
wave_2_ctls <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(tet2_chip_repeats_comorbs_drugs$AF_2 == 0.01)
wave_4_ctls <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(tet2_chip_repeats_comorbs_drugs$AF_4 == 0.01)
wave_6_ctls <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(tet2_chip_repeats_comorbs_drugs$AF_6 == 0.01)
wave_8_ctls <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(tet2_chip_repeats_comorbs_drugs$AF_8 == 0.01)
wave_9_ctls <- tet2_chip_repeats_comorbs_drugs %>% dplyr::filter(tet2_chip_repeats_comorbs_drugs$AF_9 == 0.01)

wave_2_ctls$elsa_wave <- 2
wave_4_ctls$elsa_wave <- 4
wave_6_ctls$elsa_wave <- 6
wave_9_ctls$elsa_wave <- 9

#recombine
tet2_chip_repeats_ctlwave <- do.call("rbind", list(wave_2_ctls, wave_4_ctls, wave_6_ctls, wave_9_ctls))



#import control details with sample number of control
setwd("N://My Documents/ELSA_data")
control_samples <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")
case_samples <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")
#simplify
control_samples <- subset(control_samples, select = c(Sample, idauniq, elsa_wave))

case_samples <- subset(case_samples, select = c(idauniq, Sample))
#merge
tet2_chip_repeats_ctlwave <- merge(tet2_chip_repeats_ctlwave, control_samples, by = c("idauniq", "elsa_wave"), all.x = TRUE)


#import mutation data to get mutation coordinates
setwd("N://My Documents/ELSA_data/elsa_variant_calling/")

variant_data <- read_xlsx("chip_cases_whitelist_elsa.xlsx")
#merge with sample info
variant_data <- merge(variant_data, case_samples, by = "Sample", all.x = TRUE)
#simplify
vars_simp <- subset(variant_data, select = c(idauniq, Gene.refGene, NonsynOI, Chr, Start, End))

#merge mutation data
tet2_chip_repeats_mutation_coords <- merge(tet2_chip_repeats_ctlwave, vars_simp, by = c("idauniq", "Gene.refGene", "NonsynOI"), all.x = TRUE)


#import df with locus specific depths manually added
setwd("N://My Documents/ELSA_data/upgrade/tet2_vaf_progression_linear_regression")
tet2_chip_repeats_depth <- read_xlsx("tet2_chip_repeats_mutation_anon_depth.xlsx")
#simplify
tet2_chip_repeats_depth <- subset(tet2_chip_repeats_depth, select = c(Gene.refGene, NonsynOI, Sample, growth_rate, Depth))

#add idauniq to depth df
tet2_chip_repeats_mutation_coords_simp <- subset(tet2_chip_repeats_mutation_coords, select = c(idauniq, Gene.refGene, NonsynOI, Sample))
tet2_chip_repeats_depth_idauniq <- merge(tet2_chip_repeats_depth, tet2_chip_repeats_mutation_coords_simp, by = "Sample")
tet2_chip_repeats_depth_idauniq <- unique(tet2_chip_repeats_depth_idauniq)
tet2_chip_repeats_depth_idauniq <- tet2_chip_repeats_depth_idauniq %>% dplyr::filter(!is.na(Depth))
tet2_chip_repeats_depth_idauniq <- subset(tet2_chip_repeats_depth_idauniq, select = c(idauniq, Depth))

tet2_chip_repeats_mutations_coords_depth <- merge(tet2_chip_repeats_mutation_coords, tet2_chip_repeats_depth_idauniq, by = "idauniq")

#one individual has technical repeats all confirming wave 6 sample is control - remove duplicates
tet2_chip_repeats_mutations_coords_depth1 <- tet2_chip_repeats_mutations_coords_depth[ !duplicated(tet2_chip_repeats_mutations_coords_depth$idauniq), ]

#make new variable which is maximum allele frequency that could be present at limit of detection. Filtering strategy: 20 alt reads needed
tet2_chip_repeats_mutations_coords_depth1$lod <- 20/tet2_chip_repeats_mutations_coords_depth1$Depth
#simplify
lod <- subset(tet2_chip_repeats_mutations_coords_depth1, select = c(idauniq, lod))

#make new df with lod 
tet2_chip_repeats_lod <- merge(tet2_chip_repeats_comorbs_drugs, lod, by = "idauniq", all.x = TRUE)

#drop those obs with a lod > 0.02 (locus specific depth <1000 when designated control sample)

tet2_chip_repeats_lod <- tet2_chip_repeats_lod %>% dplyr::filter(is.na(lod) | lod < 0.02)

#now replace 0.01 VAFs with lod
tet2_chip_repeats_lod$AF_2 <- ifelse(tet2_chip_repeats_lod$AF_2 == 0.01, tet2_chip_repeats_lod$lod, tet2_chip_repeats_lod$AF_2)
tet2_chip_repeats_lod$AF_4 <- ifelse(tet2_chip_repeats_lod$AF_4 == 0.01, tet2_chip_repeats_lod$lod, tet2_chip_repeats_lod$AF_4)
tet2_chip_repeats_lod$AF_6 <- ifelse(tet2_chip_repeats_lod$AF_6 == 0.01, tet2_chip_repeats_lod$lod, tet2_chip_repeats_lod$AF_6)
tet2_chip_repeats_lod$AF_8 <- ifelse(tet2_chip_repeats_lod$AF_8 == 0.01, tet2_chip_repeats_lod$lod, tet2_chip_repeats_lod$AF_8)
tet2_chip_repeats_lod$AF_9 <- ifelse(tet2_chip_repeats_lod$AF_9 == 0.01, tet2_chip_repeats_lod$lod, tet2_chip_repeats_lod$AF_9)



#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with tet2 chip mutation data
tet2_chip_repeats_comorbs_drugs_lod <- merge(tet2_chip_repeats_lod, w6_data_simp, by = "idauniq", all.x = TRUE)


#two people (row 60 and row 113) don't have drug info and survey data indicates not on statin at wave 6 timepoint
tet2_chip_repeats_comorbs_drugs_lod$lld[60] <- FALSE
tet2_chip_repeats_comorbs_drugs_lod$lld[113] <- FALSE

#two people don't have data from wave 6 on medications
tet2_chip_repeats_comorbs_drugs_lod_comp <- tet2_chip_repeats_comorbs_drugs_lod %>% dplyr::filter(!is.na(lld))

#merge mutation data with myeloid lymphoid ratio data
colnames(fbc_data_simp_w)[1] <- "idauniq"

tet2_chip_repeats_ml <- merge(tet2_chip_repeats_comorbs_drugs_lod_comp, fbc_data_simp_w, by = "idauniq", all.x = TRUE)

#correct allele frequencies for m/l ratio
tet2_chip_repeats_ml$AF_2c <- tet2_chip_repeats_ml$AF_2/tet2_chip_repeats_ml$m_l_ratio.2
tet2_chip_repeats_ml$AF_4c <- tet2_chip_repeats_ml$AF_4/tet2_chip_repeats_ml$m_l_ratio.4
tet2_chip_repeats_ml$AF_6c <- tet2_chip_repeats_ml$AF_6/tet2_chip_repeats_ml$m_l_ratio.6
tet2_chip_repeats_ml$AF_8c <- tet2_chip_repeats_ml$AF_8/tet2_chip_repeats_ml$m_l_ratio.8
tet2_chip_repeats_ml$AF_9c <- tet2_chip_repeats_ml$AF_9/tet2_chip_repeats_ml$m_l_ratio.9

#calculate growth rate- method 1
tet2_chip_repeats_ml$i_1 <- ((tet2_chip_repeats_ml$AF_9c / tet2_chip_repeats_ml$AF_2c)^{1/15} - 1)*100
tet2_chip_repeats_ml$i_2 <- ((tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_2c)^{1/13} - 1)*100
tet2_chip_repeats_ml$i_3 <- ((tet2_chip_repeats_ml$AF_9c / tet2_chip_repeats_ml$AF_4c)^{1/11} - 1)*100
tet2_chip_repeats_ml$i_4 <- ((tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_4c)^{1/9} - 1)*100
tet2_chip_repeats_ml$i_5 <- ((tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_6c)^{1/5} - 1)*100
tet2_chip_repeats_ml$i_6 <- ((tet2_chip_repeats_ml$AF_6c / tet2_chip_repeats_ml$AF_2c)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
tet2_chip_repeats_ml$growth_rate <- tet2_chip_repeats_ml$i_1
tet2_chip_repeats_ml$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate), tet2_chip_repeats_ml$i_2, tet2_chip_repeats_ml$growth_rate)
tet2_chip_repeats_ml$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate), tet2_chip_repeats_ml$i_3, tet2_chip_repeats_ml$growth_rate)
tet2_chip_repeats_ml$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate), tet2_chip_repeats_ml$i_4, tet2_chip_repeats_ml$growth_rate)
tet2_chip_repeats_ml$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate), tet2_chip_repeats_ml$i_5, tet2_chip_repeats_ml$growth_rate)
tet2_chip_repeats_ml$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate), tet2_chip_repeats_ml$i_6, tet2_chip_repeats_ml$growth_rate)


#for independent observations select single variant - keep largest growth rate
tet2_chip_repeats_ordered <- tet2_chip_repeats_ml[order(tet2_chip_repeats_ml$idauniq, abs(tet2_chip_repeats_ml$growth_rate)),]
tet2_chip_repeats_single <- tet2_chip_repeats_ordered[ !duplicated(tet2_chip_repeats_ordered$idauniq), ]


#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
tet2_chip_repeats_single <- merge(tet2_chip_repeats_single, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
tet2_chip_repeats_single$indager <- ifelse(tet2_chip_repeats_single$indager < 0, 90, tet2_chip_repeats_single$indager)

#count how many missing cholesterol measurement
nrow(tet2_chip_repeats_single[tet2_chip_repeats_single$chol <0, ])

#replace chol <0 with NA
tet2_chip_repeats_single$chol <- ifelse(tet2_chip_repeats_single$chol < 0, NA, tet2_chip_repeats_single$chol)


#combine ihd, heart_dise and stroke into cvd
tet2_chip_repeats_single$cvd <- rep(0, nrow(tet2_chip_repeats_single))
tet2_chip_repeats_single <- tet2_chip_repeats_single %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
tet2_chip_repeats_single <- tet2_chip_repeats_single %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))


#simplify to just cvd relevant variables for robust regression
tet2_chip_repeats_cvd <- subset(tet2_chip_repeats_single, select = c(idauniq, growth_rate, indager, indsex, cvd, lld, bb, chol, htn, antiplatelets, RAAS, current_smoker))


#assess missing data
md.pattern(tet2_chip_repeats_cvd)

#drop those without calculable growth rate - drop 2, 107 individuals remain
tet2_chip_repeats_cvd <- tet2_chip_repeats_cvd %>% dplyr::filter(!is.na(growth_rate))

#drop those without serum cholesterol - leaves 93 individuals
tet2_chip_repeats_cvd <- tet2_chip_repeats_cvd %>% dplyr::filter(!is.na(chol))


#run robust regression for growth rate method 1
robust_tet2_gr1 <- rlm(growth_rate ~ indager + indsex + cvd + lld + bb + chol + htn + antiplatelets + RAAS, data = tet2_chip_repeats_cvd)
summary(robust_tet2_gr1)
tidy(robust_tet2_gr1, conf.int = TRUE)
rob.pvals(robust_tet2_gr1)

gvif(robust_tet2_gr1)

#details of growth rate
summary(tet2_chip_repeats_cvd$growth_rate)
sd(tet2_chip_repeats_cvd$growth_rate)

#method 2 for calculating growth rate
#log(VAFfollowup - VAFbaseline)/(years between measurements)

#calculate growth rate
tet2_chip_repeats_ml$j_1 <- log(tet2_chip_repeats_ml$AF_9c/tet2_chip_repeats_ml$AF_2c)/15
tet2_chip_repeats_ml$j_2 <- log(tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_2c)/13
tet2_chip_repeats_ml$j_3 <- log(tet2_chip_repeats_ml$AF_9c / tet2_chip_repeats_ml$AF_4c)/11
tet2_chip_repeats_ml$j_4 <- log(tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_4c)/9
tet2_chip_repeats_ml$j_5 <- log(tet2_chip_repeats_ml$AF_8c / tet2_chip_repeats_ml$AF_6c)/5
tet2_chip_repeats_ml$j_6 <- log(tet2_chip_repeats_ml$AF_6c / tet2_chip_repeats_ml$AF_2c)/8 


#coalesce delta vaf per year into one column
tet2_chip_repeats_ml$growth_rate2 <- tet2_chip_repeats_ml$j_1
tet2_chip_repeats_ml$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate2), tet2_chip_repeats_ml$j_2, tet2_chip_repeats_ml$growth_rate2)
tet2_chip_repeats_ml$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate2), tet2_chip_repeats_ml$j_3, tet2_chip_repeats_ml$growth_rate2)
tet2_chip_repeats_ml$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate2), tet2_chip_repeats_ml$j_4, tet2_chip_repeats_ml$growth_rate2)
tet2_chip_repeats_ml$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate2), tet2_chip_repeats_ml$j_5, tet2_chip_repeats_ml$growth_rate2)
tet2_chip_repeats_ml$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_ml$growth_rate2), tet2_chip_repeats_ml$j_6, tet2_chip_repeats_ml$growth_rate2)


#for independent observations select single variant - keep largest growth rate
tet2_chip_repeats_ordered <- tet2_chip_repeats_ml[order(tet2_chip_repeats_ml$idauniq, abs(tet2_chip_repeats_ml$growth_rate2)),]
tet2_chip_repeats_single2 <- tet2_chip_repeats_ordered[ !duplicated(tet2_chip_repeats_ordered$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data

tet2_chip_repeats_single2 <- merge(tet2_chip_repeats_single2, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
tet2_chip_repeats_single2$indager <- ifelse(tet2_chip_repeats_single2$indager < 0, 90, tet2_chip_repeats_single2$indager)

#count how many missing cholesterol measurement
nrow(tet2_chip_repeats_single2[tet2_chip_repeats_single2$chol <0, ])

#replace chol <0 with NA
tet2_chip_repeats_single2$chol <- ifelse(tet2_chip_repeats_single2$chol < 0, NA, tet2_chip_repeats_single2$chol)

#combine ihd, heart_dise and stroke into cvd
tet2_chip_repeats_single2$cvd <- rep(0, nrow(tet2_chip_repeats_single2))
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))


#drop those without calculable growth rate - drop 2, 107 individuals remain
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% dplyr::filter(!is.na(growth_rate2))

#drop those without serum cholesterol - leaves 93 individuals
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% dplyr::filter(!is.na(chol))


#multiply growth_rate2 by 100 to get percentage (so comparable)
tet2_chip_repeats_single2$growth_rate2_percentage <- tet2_chip_repeats_single2$growth_rate2*100

#run robust regression - method 2

robust_tet2_gr2 <- rlm(growth_rate2_percentage ~ indager + indsex + cvd + lld + bb + chol + htn + antiplatelets + RAAS, data = tet2_chip_repeats_single2 )
summary(robust_tet2_gr2)
tidy(robust_tet2_gr2, conf.int = TRUE)
rob.pvals(robust_tet2_gr2)


