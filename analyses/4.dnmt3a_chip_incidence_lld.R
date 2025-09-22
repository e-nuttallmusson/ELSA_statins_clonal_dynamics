#script to look at whether statins have impact on incidence of dnmt3a chip
library(dplyr)
library(readxl)
library(tidyr)

#script to calculate growth rate in serially monitored TET2 and DNMT3A clones

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


#read in repeat chip mutation data with wave 8/9 comorbidities and wave 6 medications
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions_2/integrate_variants_with_elsa_data")
#this has wave 9 drugs
chip_cases_w8_9_metadata <- read.table("chip_cases_comorb_drugs_w8_9.txt",header = TRUE, sep = "\t")

#need to get info on one TET2 individual who had last follow up sample in wave 6 rather than 8 or 9
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions")
load("integrated_longitudinal_variants_and_late_wave_comorbs_drugs.RData")
w6_case <- chip_cases_comorb %>% dplyr::filter(idauniq == "[censored]")
w6_case_comorb_drugs_w8_9 <- merge(w6_case, w8_9_drug_classes, by = "idauniq", all.x = TRUE)
w6_case_comorb_drugs_w8_9 <- merge(w6_case_comorb_drugs_w8_9, wave_8_9_bloods, by = "idauniq", all.x = TRUE)

#combine
#need to remove Sample columns from chip_cases_w8_9_metadata
chip_cases_w8_9_metadata <- chip_cases_w8_9_metadata[, !grepl("^Sample", names(chip_cases_w8_9_metadata))]
chip_cases_w8_9_metadata <- rbind(chip_cases_w8_9_metadata, w6_case_comorb_drugs_w8_9)

#read in wave 6 drugs
setwd("N:/My Documents/ELSA_data/drug_regression_analysis")
wave_6_drugs <- read.table("wave_6_drug_categories.txt", header = TRUE, sep = "\t")

#drop wave 8/9 drug data
chip_cases_w8_9_comorbs <- subset(chip_cases_w8_9_metadata, select = -c(statins, ccb, ppi, acei, antiplatelets, metformin, adr_ag, tca, non_op_an, inh_cort, nsaids, bisphos, gout_rx, ur_ret, glauc, opioids, epilep, oestr, loop))

#add in wave 6 data
chip_cases_w8_9_comorbs_w6drugs <- merge(chip_cases_w8_9_comorbs, wave_6_drugs, by = "idauniq")

#keep only those with repeat measurements
chip_repeats_comorbs_drugs <- chip_cases_w8_9_comorbs_w6drugs %>% dplyr::filter(chip_cases_w8_9_comorbs_w6drugs$missing < 8)

#summarise numbers
table(chip_repeats_comorbs_drugs$Gene.refGene)

#add variable of number of mutations
chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% group_by(idauniq) %>% mutate(count = n())

#select dnmt3a repeats
dnmt3a_chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% dplyr::filter(grepl("DNMT3A", Gene.refGene))

#remove those with haem disorder
dnmt3a_chip_repeats_comorbs_drugs <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(blood_dis == 0)
dnmt3a_chip_repeats_comorbs_drugs <- as.data.frame(dnmt3a_chip_repeats_comorbs_drugs)

#create variable indicating which elsa_wave control sample was (labelled with vaf 0.01)
wave_2_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_2 == 0.01)
wave_4_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_4 == 0.01)
wave_6_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_6 == 0.01)
wave_8_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_8 == 0.01)
wave_9_ctls <- dnmt3a_chip_repeats_comorbs_drugs %>% dplyr::filter(dnmt3a_chip_repeats_comorbs_drugs$AF_9 == 0.01)

wave_2_ctls$elsa_wave <- 2
wave_4_ctls$elsa_wave <- 4
wave_6_ctls$elsa_wave <- 6
wave_8_ctls$elsa_wave <- 8
wave_9_ctls$elsa_wave <- 9

#recombine
dnmt3a_chip_repeats_ctlwave <- do.call("rbind", list(wave_2_ctls, wave_4_ctls, wave_6_ctls, wave_8_ctls, wave_9_ctls))


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
setwd("N://My Documents/ELSA_data/leukaemia_letter_statins/revisions/growth_rates")
dnmt3a_chip_repeats_depth <- read_xlsx("dnmt3a_lod_revised.xlsx")
#simplify
dnmt3a_chip_repeats_depth <- subset(dnmt3a_chip_repeats_depth, select = c(idauniq, sample, depth))
colnames(dnmt3a_chip_repeats_depth)[2] <- "Sample"

#add idauniq to depth df
dnmt3a_chip_repeats_mutation_coords_simp <- subset(dnmt3a_chip_repeats_mutation_coords, select = c(idauniq, Gene.refGene, NonsynOI, Sample))
dnmt3a_chip_repeats_depth_idauniq <- merge(dnmt3a_chip_repeats_depth, dnmt3a_chip_repeats_mutation_coords, by = "Sample")
dnmt3a_chip_repeats_depth_idauniq <- unique(dnmt3a_chip_repeats_depth_idauniq)
dnmt3a_chip_repeats_depth_idauniq <- dnmt3a_chip_repeats_depth_idauniq %>% dplyr::filter(!is.na(depth))
dnmt3a_chip_repeats_depth_idauniq <- subset(dnmt3a_chip_repeats_depth_idauniq, select = c(idauniq.x, depth))

colnames(dnmt3a_chip_repeats_depth_idauniq)[1] <- "idauniq"

dnmt3a_chip_repeats_mutations_coords_depth <- merge(dnmt3a_chip_repeats_mutation_coords, dnmt3a_chip_repeats_depth_idauniq, by = "idauniq")

#remove duplicates
dnmt3a_chip_repeats_mutations_coords_depth1 <- unique(dnmt3a_chip_repeats_mutations_coords_depth)


#make new variable which is maximum allele frequency that could be present at limit of detection. Filtering strategy: 20 alt reads needed
dnmt3a_chip_repeats_mutations_coords_depth1$lod <- 20/dnmt3a_chip_repeats_mutations_coords_depth1$depth
#simplify
lod <- subset(dnmt3a_chip_repeats_mutations_coords_depth1, select = c(idauniq, lod))

#make new df with lod 
dnmt3a_chip_repeats_lod <- merge(dnmt3a_chip_repeats_comorbs_drugs, lod, by = "idauniq", all.x = TRUE)

#drop those obs with a lod > 0.02 (locus specific depth <1000 when designated control sample)

dnmt3a_chip_repeats_lod <- dnmt3a_chip_repeats_lod %>% dplyr::filter(is.na(lod) | lod < 0.02)

#now replace 0.01 VAFs with lod
dnmt3a_chip_repeats_lod$AF_2 <- ifelse(dnmt3a_chip_repeats_lod$AF_2 == 0.01, dnmt3a_chip_repeats_lod$lod, dnmt3a_chip_repeats_lod$AF_2)
dnmt3a_chip_repeats_lod$AF_4 <- ifelse(dnmt3a_chip_repeats_lod$AF_4 == 0.01, dnmt3a_chip_repeats_lod$lod, dnmt3a_chip_repeats_lod$AF_4)
dnmt3a_chip_repeats_lod$AF_6 <- ifelse(dnmt3a_chip_repeats_lod$AF_6 == 0.01, dnmt3a_chip_repeats_lod$lod, dnmt3a_chip_repeats_lod$AF_6)
dnmt3a_chip_repeats_lod$AF_8 <- ifelse(dnmt3a_chip_repeats_lod$AF_8 == 0.01, dnmt3a_chip_repeats_lod$lod, dnmt3a_chip_repeats_lod$AF_8)
dnmt3a_chip_repeats_lod$AF_9 <- ifelse(dnmt3a_chip_repeats_lod$AF_9 == 0.01, dnmt3a_chip_repeats_lod$lod, dnmt3a_chip_repeats_lod$AF_9)

#remove duplicates
dnmt3a_chip_repeats_lod <- unique(dnmt3a_chip_repeats_lod)

#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with dnmt3a chip mutation data
dnmt3a_chip_repeats_comorbs_drugs_lod <- merge(dnmt3a_chip_repeats_lod, w6_data_simp, by = "idauniq", all.x = TRUE)

#take those with repeated sampling for dnmt3a mutations where lod < 0.02

incident_dnmt3a <- dnmt3a_chip_repeats_comorbs_drugs_lod %>% dplyr::filter(lod < 0.02)

#simplify

incident_dnmt3a_simp <- subset(incident_dnmt3a, select = c(idauniq, AF_2, AF_4, AF_6, AF_8, AF_9))

lose_clones <- c([censored])
incident_dnmt3a_simp <- incident_dnmt3a_simp %>% dplyr::filter(!idauniq %in% lose_clones)

#add time variable for wave of incident chip

incident_dnmt3a_simp$time <- NA
incident_dnmt3a_simp$time <- ifelse(incident_dnmt3a_simp$AF_4 >= 0.02 & !is.na(incident_dnmt3a_simp$AF_4), 4, incident_dnmt3a_simp$time)
incident_dnmt3a_simp$time <- ifelse(incident_dnmt3a_simp$AF_6 >= 0.02 & !is.na(incident_dnmt3a_simp$AF_6), 6, incident_dnmt3a_simp$time)
#one individual had incident chip at wave 6
incident_dnmt3a_simp$time <- ifelse(incident_dnmt3a_simp$AF_8 >= 0.02 & !is.na(incident_dnmt3a_simp$AF_8) & is.na(incident_dnmt3a_simp$time), 8, incident_dnmt3a_simp$time)
incident_dnmt3a_simp$time <- ifelse(incident_dnmt3a_simp$AF_9 >= 0.02 & !is.na(incident_dnmt3a_simp$AF_9) & is.na(incident_dnmt3a_simp$time), 9, incident_dnmt3a_simp$time)

#add outcome
incident_dnmt3a_simp$outcome <- "dnmt3a_chip"

#simplify just to idauniq, time, outcome
incident_dnmt3a_simp2 <- subset(incident_dnmt3a_simp, select = c(idauniq, time, outcome))

#remove duplicates 
incident_dnmt3a_simp2 <- unique(incident_dnmt3a_simp2)

#now need to get control data
#read in control cases
setwd("N://My Documents/ELSA_data")

control_samples <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")

control_simp <- subset(control_samples, select = c(idauniq, elsa_wave))
#drop if idauniq is na
control_simp <- control_simp %>% dplyr::filter(!is.na(idauniq))
#drop duplicates
control_simp <- unique(control_simp)

#reshape long to wide
control_simp_wide <- pivot_wider(control_simp, names_from = "elsa_wave", values_from = "elsa_wave")

#count nas across rows
control_simp_wide$missing <- rowSums(is.na(control_simp_wide))

#keep those with at least two control samples 
control_repeats <- control_simp_wide %>% dplyr::filter(missing < 4)

colnames(control_repeats) <- c("idauniq", "AF_2", "AF_9", "AF_4", "AF_8", "AF_6", "missing")
#remove missing column
control_repeats <- subset(control_repeats, select = -c(missing))

#bloods from wave 8 and 9 are mutually exclusive so assume censor date for controls is wave 9
control_repeats$time <- 9
control_repeats$outcome <- "control"

control_repeats_simp <- subset(control_repeats, select = c(idauniq, time, outcome))

#combine dnmt3a chip cases and controls
repeats_for_coxph <- rbind(incident_dnmt3a_simp2, control_repeats_simp)

#load cumulative incidence data of main comorbidities and statin treatment
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions_2/cvd_incidence")
load("revisions2_cvd_incidence.RData")

#get all longitudinal data
#combine all longitudinal data
longitudinal_data <- merge(age_longitudinal, ihd_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, stroke_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, htn_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, diabetes_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, hm_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, smoke_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
longitudinal_data <- merge(longitudinal_data, statin_longitudinal, by = "idauniq", all.x = TRUE, all.y = TRUE)
#chol - keep only those with at lease one cholesterol measurement
longitudinal_data <- merge(longitudinal_data, chol_longitudinal_nomissing, by = "idauniq")

#now need to get rid of negative numbers (censored) for age - comes in in wave 5
#wave 5
for (i in (1:nrow(longitudinal_data))) {
  if (!is.na(longitudinal_data$indager5[i]) & longitudinal_data$indager5[i] < 0) {
    longitudinal_data$indager5[i] <- longitudinal_data$indager4[i] + 2
  }
}

#wave 6
for (i in (1:nrow(longitudinal_data))) {
  if (!is.na(longitudinal_data$indager6[i]) & longitudinal_data$indager6[i] < 0) {
    longitudinal_data$indager6[i] <- longitudinal_data$indager5[i] + 2
  }
}

#wave 7
for (i in (1:nrow(longitudinal_data))) {
  if (!is.na(longitudinal_data$indager7[i]) & longitudinal_data$indager7[i] < 0) {
    longitudinal_data$indager7[i] <- longitudinal_data$indager6[i] + 2
  }
}

#wave 8
for (i in (1:nrow(longitudinal_data))) {
  if (!is.na(longitudinal_data$indager8[i]) & longitudinal_data$indager8[i] < 0) {
    longitudinal_data$indager8[i] <- longitudinal_data$indager7[i] + 2
  }
}

#wave 9
for (i in (1:nrow(longitudinal_data))) {
  if (!is.na(longitudinal_data$indager9[i]) & longitudinal_data$indager9[i] < 0) {
    longitudinal_data$indager9[i] <- longitudinal_data$indager8[i] + 2
  }
}



#reshape
#drop 'missing' var from cholesterol
longitudinal_data <- subset(longitudinal_data, select = -c(missing))
colnames(longitudinal_data)
long_data <- reshape(longitudinal_data, 
                     idvar = "idauniq",
                     varying = 2:82,
                     sep = "",
                     timevar = "wave",
                     times = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                     direction = "long")

#order
ordered_long_data <- long_data[order(long_data$idauniq, long_data$wave),]

#add in sex
indsex <- subset(wave_2_summary, select = c(idauniq, indsex))

ordered_long_data <- merge(ordered_long_data, indsex, by = "idauniq", all.x = TRUE)

#once has comorbidity then continues to have it
#ihd
ihd_cont <- as.integer(ave(ordered_long_data$ihd, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$ihd_cont <- ihd_cont

#stroke
stroke_cont <- as.integer(ave(ordered_long_data$stroke, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$stroke_cont <- stroke_cont

#htn
htn_cont <- as.integer(ave(ordered_long_data$htn, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$htn_cont <- htn_cont

#diabetes
diabetes_cont <- as.integer(ave(ordered_long_data$diabetes, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$diabetes_cont <- diabetes_cont

#hm
hm_cont <- as.integer(ave(ordered_long_data$hm, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$hm_cont <- hm_cont

#cvd
ordered_long_data$cvd <- NA
for (i in (1:nrow(ordered_long_data))) {
  ordered_long_data$cvd[i] <- ifelse(ordered_long_data$ihd_cont[i] == 1 | ordered_long_data$stroke_cont[i] == 1, 1, 0)
}

#smoke_cont
smoke_cont <- as.integer(ave(ordered_long_data$smoke, ordered_long_data$idauniq, FUN = cumsum) >= 1)
ordered_long_data$smoke_cont <- smoke_cont

#need to make a first wave variable
ordered_long_data <- ordered_long_data[order(ordered_long_data$idauniq, ordered_long_data$wave),]

ids <- repeats_for_coxph$idauniq
first_wave <- c()



for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- which(!is.na(a$indager))
  c <- min(b)
  first_wave <- append(first_wave, c)
}

firstwave_id <- cbind(ids, first_wave)

#add this information to outcome data
colnames(firstwave_id)[1] <- "idauniq"
repeats_for_coxph_fw <- merge(repeats_for_coxph, firstwave_id, by = "idauniq")

#last wave on study
last_wave <- c()



for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- which(!is.na(a$indager))
  c <- max(b)
  last_wave <- append(last_wave, c)
}

lastwave_id <- cbind(ids, last_wave)

colnames(lastwave_id)[1] <- "idauniq"
repeats_for_coxph_fw <- merge(repeats_for_coxph_fw, lastwave_id, by = "idauniq")


#if on statin in wave 3, 4, 5, 6 then say on statin as repeats largely wave 8 and 9
ids <- repeats_for_coxph$idauniq
statin <- c()



for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- ifelse(a$statin[3] == 1 | a$statin[4] == 1 | a$statin[5] == 1 | a$statin[6] == 1, 1, 0)
  c <- ifelse(a$statin[6] == 0 & is.na(b), 0, b)
  statin <- append(statin, c)
}

statin_id <- cbind(ids, statin)

#2 individuals have NA at wave 6 for statins but have never been prescribed statin
statin_id <- as.data.frame(statin_id)
statin_id[statin_id$ids == "[censored]",][2] <- 0
statin_id[statin_id$ids == "[censored]",][2] <- 0
colnames(statin_id)[1] <- "idauniq"

#combine with outcome info
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw, statin_id, by = "idauniq")

#add relevant comorbidities accrued any time before wave 7
#ihd
ihd <- c()

for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- any(a$ihd[1:6] == 1, na.rm = TRUE)
  ihd <- append(ihd, b)
}

ihd_id <- cbind(ids, ihd)
colnames(ihd_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, ihd_id, by = "idauniq")

#stroke
stroke <- c()

for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- any(a$stroke[1:6] == 1, na.rm = TRUE)
  stroke <- append(stroke, b)
}

stroke_id <- cbind(ids, stroke)
colnames(stroke_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, stroke_id, by = "idauniq")


#smoking
smoke <- c()

for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- any(a$smoke[1:6] == 1, na.rm = TRUE)
  smoke <- append(smoke, b)
}

smoke_id <- cbind(ids, smoke)
colnames(smoke_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, smoke_id, by = "idauniq")


#diabetes
diabetes <- c()

for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- any(a$diabetes[1:6] == 1, na.rm = TRUE)
  diabetes <- append(diabetes, b)
}

diabetes_id <- cbind(ids, diabetes)
colnames(diabetes_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, diabetes_id, by = "idauniq")

#htn
htn <- c()

for (i in ids){
  a <- ordered_long_data %>% dplyr::filter(idauniq == i)
  b <- any(a$htn[1:6] == 1, na.rm = TRUE)
  htn <- append(htn, b)
}

htn_id <- cbind(ids, htn)
colnames(htn_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, htn_id, by = "idauniq")

#chol
#impute cholesterol
chol <- subset(ordered_long_data, select = c(idauniq, wave, indager, chol, statin))
pred <- make.predictorMatrix(chol)


chol_imps <- mice(chol, predictorMatrix = pred)

chol_imputed <- complete(chol_imps)

#now add to df

chol <- c()

for (i in ids){
  a <- chol_imputed %>% dplyr::filter(idauniq == i)
  b <- a$chol[6] 
  chol <- append(chol, b)
}

chol_id <- cbind(ids, chol)
colnames(chol_id)[1] <- "idauniq"
repeats_for_coxph_fw_statins <- merge(repeats_for_coxph_fw_statins, chol_id, by = "idauniq")



#age and sex
wave_6_age_sex <- subset(wave_6, select = c(idauniq, indager, indsex))
repeats_for_coxph_fw_statins_age <- merge(repeats_for_coxph_fw_statins, wave_6_age_sex, by = "idauniq", all.x = TRUE)

#two individuals again missing wave 6 
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "[censored]",]$indager <- 64
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "[censored]",]$indsex <- 2

repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "[censored]",]$indager <- 64
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "[censored]",]$indsex <- 2

#make time observed variable
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 1, 2002, NA)
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 2, 2004, repeats_for_coxph_fw_statins_age$year_1)
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 3, 2006, repeats_for_coxph_fw_statins_age$year_1)
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 4, 2008, repeats_for_coxph_fw_statins_age$year_1)
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 5, 2010, repeats_for_coxph_fw_statins_age$year_1)
repeats_for_coxph_fw_statins_age$year_1 <- ifelse(repeats_for_coxph_fw_statins_age$first_wave == 6, 2012, repeats_for_coxph_fw_statins_age$year_1)


repeats_for_coxph_fw_statins_age$year_f <- ifelse(repeats_for_coxph_fw_statins_age$last_wave == 8, 2017, NA)
repeats_for_coxph_fw_statins_age$year_f <- ifelse(repeats_for_coxph_fw_statins_age$last_wave == 9, 2019, repeats_for_coxph_fw_statins_age$year_f)

#time_year
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 4, 2008, NA)
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 6, 2012, repeats_for_coxph_fw_statins_age$time_year)
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 8, 2017, repeats_for_coxph_fw_statins_age$time_year)
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 9, 2019, repeats_for_coxph_fw_statins_age$time_year)

#time in years between start and censoring
repeats_for_coxph_fw_statins_age$time_in_years <- repeats_for_coxph_fw_statins_age$time_year - repeats_for_coxph_fw_statins_age$year_1

#check for baseline differences - age, sex, cholesterol, smoking, diabetes, htn
#age 
dnmt3a_incidence_data <- repeats_for_coxph_fw_statins_age

dnmt3a_incidence_cases <- dnmt3a_incidence_data %>% dplyr::filter(outcome == "dnmt3a_chip")
dnmt3a_incidence_controls <- dnmt3a_incidence_data %>% dplyr::filter(outcome == "control")

t.test(dnmt3a_incidence_cases$indager, dnmt3a_incidence_controls$indager)
sd(dnmt3a_incidence_cases$indager)
sd(dnmt3a_incidence_controls$indager)


#sex 
sex.m <- as.matrix(table(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$indsex))
fisher.test(sex.m)

#cholesterol 
t.test(dnmt3a_incidence_cases$chol, dnmt3a_incidence_controls$chol)
sd(dnmt3a_incidence_cases$chol)
sd(dnmt3a_incidence_controls$chol)

#smoking 
smoke.m <- as.matrix(table(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$smoke))
fisher.test(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$smoke)

#diabetes 
dm.m <- as.matrix(table(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$diabetes))
fisher.test(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$diabetes)

#htn 
htn.m <- as.matrix(table(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$htn))
fisher.test(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$htn)

#cvd
dnmt3a_incidence_data$cvd <- ifelse(dnmt3a_incidence_data$ihd == 1 | dnmt3a_incidence_data$stroke == 1, 1, 0)
cvd.m <- as.matrix(table(dnmt3a_incidence_data$outcome, dnmt3a_incidence_data$cvd))
fisher.test(cvd.m)

#follow up
dnmt3a_incidence_cases$follow_up <- dnmt3a_incidence_cases$year_f - dnmt3a_incidence_cases$year_1
dnmt3a_incidence_controls$follow_up <- dnmt3a_incidence_controls$year_f - dnmt3a_incidence_controls$year_1

t.test(dnmt3a_incidence_cases$follow_up, dnmt3a_incidence_controls$follow_up)
sd(dnmt3a_incidence_cases$follow_up, na.rm = TRUE)
sd(dnmt3a_incidence_controls$follow_up, na.rm = TRUE)


#do coxph 
dnmt3a_incidence_data$status <- ifelse(dnmt3a_incidence_data$outcome == "dnmt3a_chip", TRUE, FALSE)

mod_dnmt3a_incidence <- Surv(dnmt3a_incidence_data$time_in_years, dnmt3a_incidence_data$status)

dnmt3a_coxph_model <- coxph(Surv(time_in_years, status) ~ cvd + statin, data = dnmt3a_incidence_data)
#cox.zph(dnmt3a_coxph_model)
#summary(dnmt3a_coxph_model)

dnmt3a_coxph_model2 <- coxph(Surv(time_in_years, status) ~ cvd + strata(statin), data = dnmt3a_incidence_data)
summary(dnmt3a_coxph_model2)

fit2 <- survfit(dnmt3a_coxph_model2)

ggsurvplot(fit2, data = dnmt3a_incidence_data, conf.int = TRUE, risk.table = TRUE, risk.table.col = "strata", fun = "event",
           break.time.by = 1, 
           ylab = "DNMT3A CHIP cumulative events",
           xlab = "Time (years)",
           legend.labs = c("No statin", "Statin"))




