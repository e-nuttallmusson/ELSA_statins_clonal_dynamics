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
library(sfsmisc)

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
table(dups)
#dups
#FALSE  TRUE 
#6456     2 

#check for duplications in wave 6
fbc_data_simp6 <- fbc_data_simp %>% dplyr::filter(wave == 6)
dups <- duplicated(fbc_data_simp6$ID)

#duplicated wave 6 ids with NA values:
fbc_data_simp[14086,]
fbc_data_simp[14302,]

#drop these
fbc_data_simp_dedup <- fbc_data_simp %>% filter(!row_number() %in% c(14086, 14302))

#reshape data long to wide
fbc_data_simp_w <- reshape(fbc_data_simp_dedup, idvar = "ID", v.names = "m_l_ratio", timevar = "wave", direction = "wide")


#calculate growth rate - method 1
setwd("N:/My Documents/ELSA_data/upgrade/tet2_vaf_progression_linear_regression/")

#read in repeat chip mutation data with wave 8/9 comorbidities and wave 6 medications

chip_repeats_comorbs_drugs <- read_xlsx("chip_repeats_all_comorbs_wave6_drug.xlsx")

#add variable of number of mutations
chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% group_by(idauniq) %>% mutate(count = n())

#select dnmt3a repeats
dnmt3a_chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% dplyr::filter(grepl("DNMT3A", Gene.refGene))

#remove haem disorder
dnmt3a_chip_repeats_comorbs_drugs <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(blood_dis == 0)
dnmt3a_chip_repeats_comorbs_drugs <- as.data.frame(dnmt3a_chip_repeats_comorbs_drugs)

#create variable indicating which elsa_wave control sample was - labelled with placeholder VAF 0.01
wave_2_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_2 == 0.01)
wave_4_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_4 == 0.01)
wave_6_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_6 == 0.01)
wave_8_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_8 == 0.01)
wave_9_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_9 == 0.01)

wave_2_ctls$elsa_wave <- 2
wave_4_ctls$elsa_wave <- 4
wave_6_ctls$elsa_wave <- 6
wave_9_ctls$elsa_wave <- 9

#recombine
dnmt3a_chip_repeats_ctlwave <- do.call("rbind", list(wave_2_ctls, wave_4_ctls, wave_6_ctls, wave_9_ctls))


#import control details with sample number of control
setwd("N://My Documents/ELSA_data")
control_samples <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")
case_samples <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")
#simplify
control_samples <- subset(control_samples, select = c(Sample, idauniq, elsa_wave))

case_samples <- subset(case_samples, select = c(idauniq, Sample))
#merge
dnmt3a_chip_repeats_ctlwave <- merge(dnmt3a_chip_repeats_ctlwave, control_samples, by = c("idauniq", "elsa_wave"), all.x = TRUE)


#import mutation data to get mutation coordinates
setwd("N://My Documents/ELSA_data/elsa_variant_calling/")

variant_data <- read_xlsx("chip_cases_whitelist_elsa.xlsx")
#merge with sample info
variant_data <- merge(variant_data, case_samples, by = "Sample", all.x = TRUE)
#simplify
vars_simp <- subset(variant_data, select = c(idauniq, Gene.refGene, NonsynOI, Chr, Start, End))

#merge mutation data
dnmt3a_chip_repeats_mutation_coords <- merge(dnmt3a_chip_repeats_ctlwave, vars_simp, by = c("idauniq", "Gene.refGene", "NonsynOI"), all.x = TRUE)

#import df with locus specific depths manually added
setwd("N://My Documents/ELSA_data/upgrade/statin_focussed_analyses_jan25")
dnmt3a_chip_repeats_depth <- read_xlsx("dnmt3a_lod.xlsx")
#simplify
dnmt3a_chip_repeats_depth <- subset(dnmt3a_chip_repeats_depth, select = c(idauniq, depth))

dnmt3a_chip_repeats_mutation_coords_depth <- merge(dnmt3a_chip_repeats_mutation_coords, dnmt3a_chip_repeats_depth, by = "idauniq", all.x = TRUE)

#make new variable which is maximum allele frequency that could be present at limit of detection. Filtering strategy: 20 alt reads needed
dnmt3a_chip_repeats_mutation_coords_depth$lod <- 20/dnmt3a_chip_repeats_mutation_coords_depth$depth

#drop those obs with a lod > 0.02 (locus specific depth <1000 when designated control sample)

dnmt3a_chip_repeats_lod <- dnmt3a_chip_repeats_mutation_coords_depth %>% dplyr::filter(lod < 0.02)
dnmt3a_lod <- subset(dnmt3a_chip_repeats_lod, select= c(idauniq, lod))

#split dnmt3a chip repeats into two - to replace lod and those where not required
dnmt3a_chip_repeats_comp <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(!idauniq %in% dnmt3a_chip_repeats_ctlwave$idauniq)
dnmt3a_chip_repeats_rep <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(idauniq %in% dnmt3a_chip_repeats_lod$idauniq)

dnmt3a_chip_repeats_rep <- merge(dnmt3a_chip_repeats_rep, dnmt3a_lod, by = "idauniq")

#now replace 0.01 VAFs with lod
dnmt3a_chip_repeats_rep$AF_2 <- ifelse(dnmt3a_chip_repeats_rep$AF_2 == 0.01, dnmt3a_chip_repeats_rep$lod, dnmt3a_chip_repeats_rep$AF_2)
dnmt3a_chip_repeats_rep$AF_4 <- ifelse(dnmt3a_chip_repeats_rep$AF_4 == 0.01, dnmt3a_chip_repeats_rep$lod, dnmt3a_chip_repeats_rep$AF_4)
dnmt3a_chip_repeats_rep$AF_6 <- ifelse(dnmt3a_chip_repeats_rep$AF_6 == 0.01, dnmt3a_chip_repeats_rep$lod, dnmt3a_chip_repeats_rep$AF_6)
dnmt3a_chip_repeats_rep$AF_8 <- ifelse(dnmt3a_chip_repeats_rep$AF_8 == 0.01, dnmt3a_chip_repeats_rep$lod, dnmt3a_chip_repeats_rep$AF_8)
dnmt3a_chip_repeats_rep$AF_9 <- ifelse(dnmt3a_chip_repeats_rep$AF_9 == 0.01, dnmt3a_chip_repeats_rep$lod, dnmt3a_chip_repeats_rep$AF_9)

#drop lod var now incorporated
dnmt3a_chip_repeats_rep <- subset(dnmt3a_chip_repeats_rep, select = -c(lod))

#merge with complete cases
dnmt3a_chip_repeats_lod <- rbind(dnmt3a_chip_repeats_comp, dnmt3a_chip_repeats_rep)

#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with te2 chip mutation data
dnmt3a_chip_repeats_comorbs_drugs_lod <- merge(dnmt3a_chip_repeats_lod, w6_data_simp, by = "idauniq", all.x = TRUE)


#two people don't have data from wave 6 on medications
dnmt3a_chip_repeats_comorbs_drugs_lod_comp <- dnmt3a_chip_repeats_comorbs_drugs_lod %>% dplyr::filter(!is.na(lld))

#merge mutation data with myeloid lymphoid ratio data
colnames(fbc_data_simp_w)[1] <- "idauniq"

dnmt3a_chip_repeats_ml <- merge(dnmt3a_chip_repeats_comorbs_drugs_lod_comp, fbc_data_simp_w, by = "idauniq", all.x = TRUE)

#correct allele frequencies for m/l ratio
dnmt3a_chip_repeats_ml$AF_2c <- dnmt3a_chip_repeats_ml$AF_2/dnmt3a_chip_repeats_ml$m_l_ratio.2
dnmt3a_chip_repeats_ml$AF_4c <- dnmt3a_chip_repeats_ml$AF_4/dnmt3a_chip_repeats_ml$m_l_ratio.4
dnmt3a_chip_repeats_ml$AF_6c <- dnmt3a_chip_repeats_ml$AF_6/dnmt3a_chip_repeats_ml$m_l_ratio.6
dnmt3a_chip_repeats_ml$AF_8c <- dnmt3a_chip_repeats_ml$AF_8/dnmt3a_chip_repeats_ml$m_l_ratio.8
dnmt3a_chip_repeats_ml$AF_9c <- dnmt3a_chip_repeats_ml$AF_9/dnmt3a_chip_repeats_ml$m_l_ratio.9


#calculate growth rate - method 1
dnmt3a_chip_repeats_ml$i_1 <- ((dnmt3a_chip_repeats_ml$AF_9c / dnmt3a_chip_repeats_ml$AF_2c)^{1/15} - 1)*100
dnmt3a_chip_repeats_ml$i_2 <- ((dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_2c)^{1/13} - 1)*100
dnmt3a_chip_repeats_ml$i_3 <- ((dnmt3a_chip_repeats_ml$AF_9c / dnmt3a_chip_repeats_ml$AF_4c)^{1/11} - 1)*100
dnmt3a_chip_repeats_ml$i_4 <- ((dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_4c)^{1/9} - 1)*100
dnmt3a_chip_repeats_ml$i_5 <- ((dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_6c)^{1/5} - 1)*100
dnmt3a_chip_repeats_ml$i_6 <- ((dnmt3a_chip_repeats_ml$AF_6c / dnmt3a_chip_repeats_ml$AF_2c)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
dnmt3a_chip_repeats_ml$growth_rate <- dnmt3a_chip_repeats_ml$i_1
dnmt3a_chip_repeats_ml$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate), dnmt3a_chip_repeats_ml$i_2, dnmt3a_chip_repeats_ml$growth_rate)
dnmt3a_chip_repeats_ml$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate), dnmt3a_chip_repeats_ml$i_3, dnmt3a_chip_repeats_ml$growth_rate)
dnmt3a_chip_repeats_ml$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate), dnmt3a_chip_repeats_ml$i_4, dnmt3a_chip_repeats_ml$growth_rate)
dnmt3a_chip_repeats_ml$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate), dnmt3a_chip_repeats_ml$i_5, dnmt3a_chip_repeats_ml$growth_rate)
dnmt3a_chip_repeats_ml$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate), dnmt3a_chip_repeats_ml$i_6, dnmt3a_chip_repeats_ml$growth_rate)


#for independent observations select single variant - keep largest growth rate
dnmt3a_chip_repeats_ordered <- dnmt3a_chip_repeats_ml[order(dnmt3a_chip_repeats_ml$idauniq, abs(dnmt3a_chip_repeats_ml$growth_rate)),]
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_ordered[ !duplicated(dnmt3a_chip_repeats_ordered$idauniq), ]


#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
dnmt3a_chip_repeats_single <- merge(dnmt3a_chip_repeats_single, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
dnmt3a_chip_repeats_single$indager <- ifelse(dnmt3a_chip_repeats_single$indager < 0, 90, dnmt3a_chip_repeats_single$indager)

#count how many missing cholesterol measurement
nrow(dnmt3a_chip_repeats_single[dnmt3a_chip_repeats_single$chol <0, ])

#replace chol <0 with NA
dnmt3a_chip_repeats_single$chol <- ifelse(dnmt3a_chip_repeats_single$chol < 0, NA, dnmt3a_chip_repeats_single$chol)

#combine ihd, heart_dise and stroke into cvd
dnmt3a_chip_repeats_single$cvd <- rep(0, nrow(dnmt3a_chip_repeats_single))
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_single %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_single %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))


#simplify to just cvd relevant variables that I would use in linear regression 
dnmt3a_chip_repeats_cvd <- subset(dnmt3a_chip_repeats_single, select = c(idauniq, growth_rate, lld, indager, indsex, cvd, htn, chol, bb, current_smoker))


#assess missing data
md.pattern(dnmt3a_chip_repeats_cvd)

#drop those without calculable growth rate 
dnmt3a_chip_repeats_cvd <- dnmt3a_chip_repeats_single %>% dplyr::filter(!is.na(growth_rate))

#drop those without serum cholesterol 
dnmt3a_chip_repeats_cvd <- dnmt3a_chip_repeats_cvd %>% dplyr::filter(!is.na(chol))


#run robust regression

dnmt3a_rqfit1 <- rlm(growth_rate ~ indager + indsex + cvd + lld + bb + chol + htn + antiplatelets + RAAS + current_smoker, data = dnmt3a_chip_repeats_cvd)
summary(dnmt3a_rqfit1)
tidy(dnmt3a_rqfit1, conf.int = TRUE)
rob.pvals(dnmt3a_rqfit1)

gvif(dnmt3a_rqfit1)

#details of growth rate
summary(dnmt3a_chip_repeats_cvd$growth_rate)
sd(dnmt3a_chip_repeats_cvd$growth_rate)


#method 2 for calculating growth rate
#log(VAFfollowup - VAFbaseline)/(years between measurements)

#calculate growth rate
dnmt3a_chip_repeats_ml$j_1 <- log(dnmt3a_chip_repeats_ml$AF_9c/dnmt3a_chip_repeats_ml$AF_2c)/15
dnmt3a_chip_repeats_ml$j_2 <- log(dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_2c)/13
dnmt3a_chip_repeats_ml$j_3 <- log(dnmt3a_chip_repeats_ml$AF_9c / dnmt3a_chip_repeats_ml$AF_4c)/11
dnmt3a_chip_repeats_ml$j_4 <- log(dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_4c)/9
dnmt3a_chip_repeats_ml$j_5 <- log(dnmt3a_chip_repeats_ml$AF_8c / dnmt3a_chip_repeats_ml$AF_6c)/5
dnmt3a_chip_repeats_ml$j_6 <- log(dnmt3a_chip_repeats_ml$AF_6c / dnmt3a_chip_repeats_ml$AF_2c)/8 


#coalesce delta vaf per year into one column
dnmt3a_chip_repeats_ml$growth_rate2 <- dnmt3a_chip_repeats_ml$j_1
dnmt3a_chip_repeats_ml$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate2), dnmt3a_chip_repeats_ml$j_2, dnmt3a_chip_repeats_ml$growth_rate2)
dnmt3a_chip_repeats_ml$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate2), dnmt3a_chip_repeats_ml$j_3, dnmt3a_chip_repeats_ml$growth_rate2)
dnmt3a_chip_repeats_ml$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate2), dnmt3a_chip_repeats_ml$j_4, dnmt3a_chip_repeats_ml$growth_rate2)
dnmt3a_chip_repeats_ml$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate2), dnmt3a_chip_repeats_ml$j_5, dnmt3a_chip_repeats_ml$growth_rate2)
dnmt3a_chip_repeats_ml$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_ml$growth_rate2), dnmt3a_chip_repeats_ml$j_6, dnmt3a_chip_repeats_ml$growth_rate2)


#for independent observations select single variant - keep largest growth rate
dnmt3a_chip_repeats_ordered <- dnmt3a_chip_repeats_ml[order(dnmt3a_chip_repeats_ml$idauniq, abs(dnmt3a_chip_repeats_ml$growth_rate2)),]
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_ordered[ !duplicated(dnmt3a_chip_repeats_ordered$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data

dnmt3a_chip_repeats_single2 <- merge(dnmt3a_chip_repeats_single2, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace neg numbers (missing data) with median for age and cholesterol
dnmt3a_chip_repeats_single2$indager <- ifelse(dnmt3a_chip_repeats_single2$indager < 0, 90, dnmt3a_chip_repeats_single2$indager)

#count how many missing cholesterol measurement
nrow(dnmt3a_chip_repeats_single2[dnmt3a_chip_repeats_single2$chol <0, ])

#replace chol <0 with NA
dnmt3a_chip_repeats_single2$chol <- ifelse(dnmt3a_chip_repeats_single2$chol < 0, NA, dnmt3a_chip_repeats_single2$chol)


#combine ihd, heart_dise and stroke into cvd
dnmt3a_chip_repeats_single2$cvd <- rep(0, nrow(dnmt3a_chip_repeats_single2))
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))


#drop those without calculable growth rate 
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% dplyr::filter(!is.na(growth_rate2))

#drop those without serum cholesterol 
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% dplyr::filter(!is.na(chol))


#multiply growth_rate2 by 100 to get percentage (so comparable)
dnmt3a_chip_repeats_single2$growth_rate2_percentage <- dnmt3a_chip_repeats_single2$growth_rate2*100

#run robust regression

dnmt3a_rqfit2 <- rlm(growth_rate2_percentage ~ indager + indsex + cvd + lld + bb + chol + htn + antiplatelets + RAAS + current_smoker, data = dnmt3a_chip_repeats_single2)
summary(dnmt3a_rqfit2)
tidy(dnmt3a_rqfit2, conf.int = TRUE)
rob.pvals(dnmt3a_rqfit2)


