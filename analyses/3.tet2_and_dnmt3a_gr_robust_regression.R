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
library(boot)   
library(purrr)
library(tibble)


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

fbc_data_simp <- subset(fbc_data, select = c(ID, wave, WBC_num, lymph_num, m_l_ratio))


#check for duplications in wave 6
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
fbc_data_simp_w <- reshape(fbc_data_simp_dedup, idvar = "ID", v.names = c("WBC_num", "lymph_num","m_l_ratio"), timevar = "wave", direction = "wide")


#calculate growth rate

#read in repeat chip mutation data with wave 8/9 comorbidities and wave 6 medications
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions_2/integrate_variants_with_elsa_data")
#this has wave 9 drugs
chip_cases_w8_9_metadata <- read.table("chip_cases_comorb_drugs_w8_9.txt",header = TRUE, sep = "\t")

#need to get info on one TET2 individual who had last follow up sample in wave 6 rather than 8 or 9
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions")
load("integrated_longitudinal_variants_and_late_wave_comorbs_drugs.RData")
w6_case <- chip_cases_comorb %>% dplyr::filter(idauniq == [censored])
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


#add variable of number of mutations
chip_repeats_comorbs_drugs <- chip_repeats_comorbs_drugs %>% group_by(idauniq) %>% mutate(count = n())


####################TET2
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

#add idauniq to depth df
tet2_chip_repeats_mutation_coords_simp <- subset(tet2_chip_repeats_mutation_coords, select = c(idauniq, Gene.refGene, NonsynOI, Sample))
tet2_chip_repeats_depth_idauniq <- merge(tet2_chip_repeats_depth, tet2_chip_repeats_mutation_coords_simp, by = "Sample")
tet2_chip_repeats_depth_idauniq <- unique(tet2_chip_repeats_depth_idauniq)
tet2_chip_repeats_depth_idauniq <- tet2_chip_repeats_depth_idauniq %>% dplyr::filter(!is.na(Depth))
tet2_chip_repeats_depth_idauniq <- subset(tet2_chip_repeats_depth_idauniq, select = c(idauniq, Depth))

tet2_chip_repeats_mutations_coords_depth <- merge(tet2_chip_repeats_mutation_coords, tet2_chip_repeats_depth_idauniq, by = "idauniq", all.x = TRUE)

#remove duplicates
tet2_chip_repeats_mutations_coords_depth1 <- unique(tet2_chip_repeats_mutations_coords_depth)

#one individual has technical repeats all confirming wave 6 sample is control - remove duplicates
tet2_chip_repeats_mutations_coords_depth1 <- subset(tet2_chip_repeats_mutations_coords_depth1, select = -c(Sample))


tet2_chip_repeats_mutations_coords_depth1 <- unique(tet2_chip_repeats_mutations_coords_depth1)

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

#remove duplicates
tet2_chip_repeats_lod <- unique(tet2_chip_repeats_lod)

#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with tet2 chip mutation data
tet2_chip_repeats_comorbs_drugs_lod <- merge(tet2_chip_repeats_lod, w6_data_simp, by = "idauniq", all.x = TRUE)

#ensure all cases complete with medication info
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
tet2_chip_repeats_ordered <- tet2_chip_repeats_ml[order(tet2_chip_repeats_ml$idauniq, abs(tet2_chip_repeats_ml$growth_rate), decreasing = TRUE),]
tet2_chip_repeats_single <- tet2_chip_repeats_ordered[ !duplicated(tet2_chip_repeats_ordered$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
tet2_chip_repeats_single <- merge(tet2_chip_repeats_single, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
tet2_chip_repeats_single$indager <- ifelse(tet2_chip_repeats_single$indager < 0, 90, tet2_chip_repeats_single$indager)

#replate missing chol data with na
tet2_chip_repeats_single$chol <- ifelse(tet2_chip_repeats_single$chol.y < 0, NA, tet2_chip_repeats_single$chol.y)


#combine ihd, heart_dise and stroke into cvd
tet2_chip_repeats_single$cvd <- rep(0, nrow(tet2_chip_repeats_single))
tet2_chip_repeats_single <- tet2_chip_repeats_single %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
tet2_chip_repeats_single <- tet2_chip_repeats_single %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
tet2_chip_repeats_single$initial_vaf <- tet2_chip_repeats_single$AF_2c
tet2_chip_repeats_single$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single$initial_vaf), tet2_chip_repeats_single$AF_4c, tet2_chip_repeats_single$initial_vaf)
tet2_chip_repeats_single$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single$initial_vaf), tet2_chip_repeats_single$AF_6c, tet2_chip_repeats_single$initial_vaf)


#simplify initial vaf as will have outliers - could take log
tet2_chip_repeats_single$log_initial_vaf <- log(tet2_chip_repeats_single$initial_vaf)



#simplify to just cvd relevant variables for robust regression
tet2_chip_repeats_cvd <- subset(tet2_chip_repeats_single, select = c("idauniq", "growth_rate", "indager", "indsex", "cvd", "lld", "bb.y", "chol", "htn", "antiplatelets", "RAAS", "current_smoker", "log_initial_vaf", "count", "diabetes", "hscrp.y"))


#assess missing data
md.pattern(tet2_chip_repeats_cvd)

#drop those without cholesterol

tet2_chip_repeats_cvd <- tet2_chip_repeats_cvd %>% dplyr::filter(!is.na(chol))


#drop those without calculable growth rate 
tet2_chip_repeats_cvd <- tet2_chip_repeats_cvd %>% dplyr::filter(!is.na(growth_rate))


#run robust regression for growth rate method 1
robust_tet2_gr1 <- rlm(growth_rate ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count, data = tet2_chip_repeats_cvd)
summary(robust_tet2_gr1)
tidy(robust_tet2_gr1, conf.int = TRUE)
rob.pvals(robust_tet2_gr1)

gvif(robust_tet2_gr1)


# DATA:
dat <- tet2_chip_repeats_cvd  


form <- growth_rate ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count

#Fit once to get baseline estimates & term names
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))   
p <- length(coef_names)

#Bootstrap function: resample rows, refit rlm, return aligned coefs
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_results_perc <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_results_perc

#gvif
gvif(fit0)






#details of growth rate
summary(tet2_chip_repeats_cvd$growth_rate)
sd(tet2_chip_repeats_cvd$growth_rate)


#now look for additional associations outlined by reviewers adjusted for age and initial vaf
#current smoker - now underpowered for this analysis as only 4 individuals
table(tet2_chip_repeats_cvd$bb) #underpowered 9/84
table(tet2_chip_repeats_cvd$RAAS) #26/84

#raas
#DATA:
dat <- tet2_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + RAAS

#Fit once to get baseline estimates & term names
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))        
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_raas <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_raas


#antipletlets given association in logistic regression

table(tet2_chip_repeats_cvd$antiplatelets)

#DATA:
dat <- tet2_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + antiplatelets

#Fit once to get baseline estimates & term names
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))      
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_ap <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_ap



#hscrp
#DATA:
dat <- tet2_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + hscrp.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_crp <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_crp



#diabetes
#DATA:
dat <- tet2_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + diabetes

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_dm <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_dm



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
tet2_chip_repeats_ordered2 <- tet2_chip_repeats_ml[order(tet2_chip_repeats_ml$idauniq, abs(tet2_chip_repeats_ml$growth_rate2), decreasing = TRUE),]

tet2_chip_repeats_ordered2 <- tet2_chip_repeats_ordered2 %>% group_by(idauniq) %>% mutate(tet2_count = n())

tet2_chip_repeats_single2 <- tet2_chip_repeats_ordered2[ !duplicated(tet2_chip_repeats_ordered2$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data

tet2_chip_repeats_single2 <- merge(tet2_chip_repeats_single2, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
tet2_chip_repeats_single2$indager <- ifelse(tet2_chip_repeats_single2$indager < 0, 90, tet2_chip_repeats_single2$indager)

#count how many missing cholesterol measurement
nrow(tet2_chip_repeats_single2[tet2_chip_repeats_single2$chol.y <0, ])

#replace chol <0 with NA
tet2_chip_repeats_single2$chol.y <- ifelse(tet2_chip_repeats_single2$chol.y < 0, NA, tet2_chip_repeats_single2$chol.y)

#combine ihd, heart_dise and stroke into cvd
tet2_chip_repeats_single2$cvd <- rep(0, nrow(tet2_chip_repeats_single2))
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
tet2_chip_repeats_single2$initial_vaf <- tet2_chip_repeats_single2$AF_2c
tet2_chip_repeats_single2$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single2$initial_vaf), tet2_chip_repeats_single2$AF_4c, tet2_chip_repeats_single2$initial_vaf)
tet2_chip_repeats_single2$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single2$initial_vaf), tet2_chip_repeats_single2$AF_6c, tet2_chip_repeats_single2$initial_vaf)

#take log of initial VAF so more interpretable
tet2_chip_repeats_single2$log_initial_vaf <- log(tet2_chip_repeats_single2$initial_vaf)

#drop those without calculable growth rate 
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% dplyr::filter(!is.na(growth_rate2))

#drop those without serum cholesterol 
tet2_chip_repeats_single2 <- tet2_chip_repeats_single2 %>% dplyr::filter(!is.na(chol.y))


#multiply growth_rate2 by 100 to get percentage (so comparable)
tet2_chip_repeats_single2$growth_rate2_percentage <- tet2_chip_repeats_single2$growth_rate2*100

#run robust regression - method 2

#DATA:
dat <- tet2_chip_repeats_single2  


form <- growth_rate2_percentage ~ indager + indsex + cvd + lld + chol.y + htn + log_initial_vaf + count

#Fit once to get baseline estimates & term names
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_results_perc <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_results_perc

 
#look at antiplatelets given association in logistic regression
#DATA:
dat <- tet2_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + antiplatelets

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_ap <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_ap

#RAAS
#DATA:
dat <- tet2_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + RAAS

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_raas <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_raas


#now look for additional associations outlined by reviewers adjusted for age and initial vaf

#hscrp

#DATA:
dat <- tet2_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + hscrp.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_crp <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_crp





#diabetes
form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + diabetes

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_dm <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_dm


##########################DNMT3A

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

#ensure all cases complete with medication info
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

#calculate growth rate- method 1
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
dnmt3a_chip_repeats_ordered <- dnmt3a_chip_repeats_ml[order(dnmt3a_chip_repeats_ml$idauniq, abs(dnmt3a_chip_repeats_ml$growth_rate), decreasing = TRUE),]
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_ordered[ !duplicated(dnmt3a_chip_repeats_ordered$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
dnmt3a_chip_repeats_single <- merge(dnmt3a_chip_repeats_single, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
dnmt3a_chip_repeats_single$indager <- ifelse(dnmt3a_chip_repeats_single$indager < 0, 90, dnmt3a_chip_repeats_single$indager)


#combine ihd, heart_dise and stroke into cvd
dnmt3a_chip_repeats_single$cvd <- rep(0, nrow(dnmt3a_chip_repeats_single))
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_single %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
dnmt3a_chip_repeats_single <- dnmt3a_chip_repeats_single %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
dnmt3a_chip_repeats_single$initial_vaf <- dnmt3a_chip_repeats_single$AF_2c
dnmt3a_chip_repeats_single$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single$initial_vaf), dnmt3a_chip_repeats_single$AF_4c, dnmt3a_chip_repeats_single$initial_vaf)
dnmt3a_chip_repeats_single$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single$initial_vaf), dnmt3a_chip_repeats_single$AF_6c, dnmt3a_chip_repeats_single$initial_vaf)

#take log of initial vaf
dnmt3a_chip_repeats_single$log_initial_vaf <- log(dnmt3a_chip_repeats_single$initial_vaf)


#simplify to just cvd relevant variables for robust regression
dnmt3a_chip_repeats_cvd <- subset(dnmt3a_chip_repeats_single, select = c("idauniq", "growth_rate", "indager", "indsex", "cvd", "lld", "bb.y", "chol.y", "htn", "antiplatelets", "RAAS", "current_smoker", "log_initial_vaf", "count", "diabetes", "hscrp.y"))


#assess missing data
md.pattern(dnmt3a_chip_repeats_cvd)


#drop those without calculable growth rate 
dnmt3a_chip_repeats_cvd <- dnmt3a_chip_repeats_cvd %>% dplyr::filter(!is.na(growth_rate))

#drop those without cholesterol
dnmt3a_chip_repeats_cvd <- dnmt3a_chip_repeats_cvd %>% dplyr::filter(chol.y > 0)


#run robust regression for growth rate method 1
#DATA:
dat <- dnmt3a_chip_repeats_cvd  


form <- growth_rate ~ indager + indsex + cvd + lld + chol.y + htn + log_initial_vaf + count

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_results_perc <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_results_perc

#gvif
gvif(fit0)


#details of growth rate
summary(dnmt3a_chip_repeats_cvd$growth_rate)
sd(dnmt3a_chip_repeats_cvd$growth_rate)

#antiplatelets
table(dnmt3a_chip_repeats_cvd$antiplatelets) #21/113

#DATA:
dat <- dnmt3a_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + antiplatelets

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_ap <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_ap

#raas

table(dnmt3a_chip_repeats_cvd$RAAS) #30/113

#DATA:
dat <- dnmt3a_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + RAAS

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_raas <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_raas

#bb
#DATA:
dat <- dnmt3a_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + bb.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_bb <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_bb

#now look for additional associations outlined by reviewers adjusted for age and initial vaf
#current smoker
#DATA:
dat <- dnmt3a_chip_repeats_cvd


form <- growth_rate ~ indager + indsex + log_initial_vaf + count + current_smoker

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_smoke <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_smoke


#hscrp
form <- growth_rate ~ indager + indsex + log_initial_vaf + count + hscrp.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_crp <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_crp


#diabetes
form <- growth_rate ~ indager + indsex + log_initial_vaf + count + diabetes

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_dm <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_dm

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
dnmt3a_chip_repeats_ordered2 <- dnmt3a_chip_repeats_ml[order(dnmt3a_chip_repeats_ml$idauniq, abs(dnmt3a_chip_repeats_ml$growth_rate2), decreasing = TRUE),]

dnmt3a_chip_repeats_ordered2 <- dnmt3a_chip_repeats_ordered2 %>% group_by(idauniq) %>% mutate(dnmt3a_count = n())

dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_ordered2[ !duplicated(dnmt3a_chip_repeats_ordered2$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data

dnmt3a_chip_repeats_single2 <- merge(dnmt3a_chip_repeats_single2, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
dnmt3a_chip_repeats_single2$indager <- ifelse(dnmt3a_chip_repeats_single2$indager < 0, 90, dnmt3a_chip_repeats_single2$indager)

#count how many missing cholesterol measurement
nrow(dnmt3a_chip_repeats_single2[dnmt3a_chip_repeats_single2$chol.y <0, ])

#replace chol <0 with NA
dnmt3a_chip_repeats_single2$chol.y <- ifelse(dnmt3a_chip_repeats_single2$chol.y < 0, NA, dnmt3a_chip_repeats_single2$chol.y)

#combine ihd, heart_dise and stroke into cvd
dnmt3a_chip_repeats_single2$cvd <- rep(0, nrow(dnmt3a_chip_repeats_single2))
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
dnmt3a_chip_repeats_single2$initial_vaf <- dnmt3a_chip_repeats_single2$AF_2c
dnmt3a_chip_repeats_single2$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single2$initial_vaf), dnmt3a_chip_repeats_single2$AF_4c, dnmt3a_chip_repeats_single2$initial_vaf)
dnmt3a_chip_repeats_single2$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single2$initial_vaf), dnmt3a_chip_repeats_single2$AF_6c, dnmt3a_chip_repeats_single2$initial_vaf)

#take log of initial vaf
dnmt3a_chip_repeats_single2$log_initial_vaf <- log(dnmt3a_chip_repeats_single2$initial_vaf)

#drop those without calculable growth rate 
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% dplyr::filter(!is.na(growth_rate2))

#drop those without serum cholesterol - 
dnmt3a_chip_repeats_single2 <- dnmt3a_chip_repeats_single2 %>% dplyr::filter(!is.na(chol.y))

#review how many have > 1 mutations
table(dnmt3a_chip_repeats_single2$dnmt3a_count)

#multiply growth_rate2 by 100 to get percentage (so comparable)
dnmt3a_chip_repeats_single2$growth_rate2_percentage <- dnmt3a_chip_repeats_single2$growth_rate2*100

#run robust regression - method 2

#DATA:
dat <- dnmt3a_chip_repeats_single2  


form <- growth_rate2_percentage ~ indager + indsex + cvd + lld + chol.y + htn + log_initial_vaf + count

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_results_perc <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_results_perc



#other cvd meds
#antiplatelets
#DATA:
dat <- dnmt3a_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + antiplatelets

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_ap <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_ap


#bb
#DATA:
dat <- dnmt3a_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + bb.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_bb <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_bb

#raas
#DATA:
dat <- dnmt3a_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + RAAS

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_raas <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_raas


#now look for additional associations outlined by reviewers adjusted for age and initial vaf
#current smoker
#DATA:
dat <- dnmt3a_chip_repeats_single2


form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + current_smoker

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_smoke <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_smoke



#hscrp
form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + hscrp.y

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_crp <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_crp






#diabetes
form <- growth_rate2_percentage ~ indager + indsex + log_initial_vaf + count + diabetes

#Fit once to get baseline estimates
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_diabetes <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_diabetes


##########sensitivity analysis to replace any lod vaf <0.01 with 0.01

#################################
#drop those obs with a lod > 0.02 (locus specific depth <1000 when designated control sample)
#make new df with lod 
#redefine lod
lod <- subset(tet2_chip_repeats_mutations_coords_depth1, select = c(idauniq, lod))

tet2_chip_repeats_lod2 <- merge(tet2_chip_repeats_comorbs_drugs, lod, by = "idauniq", all.x = TRUE)

tet2_chip_repeats_lod2 <- tet2_chip_repeats_lod2 %>% dplyr::filter(is.na(lod) | lod < 0.02)

#now replace 0.01 VAFs with lod
tet2_chip_repeats_lod2$AF_2 <- ifelse(tet2_chip_repeats_lod2$AF_2 == 0.01 & !is.na(tet2_chip_repeats_lod2$lod) & tet2_chip_repeats_lod2$lod > 0.01, tet2_chip_repeats_lod2$lod, tet2_chip_repeats_lod2$AF_2)
tet2_chip_repeats_lod2$AF_4 <- ifelse(tet2_chip_repeats_lod2$AF_4 == 0.01 & !is.na(tet2_chip_repeats_lod2$lod) & tet2_chip_repeats_lod2$lod > 0.01, tet2_chip_repeats_lod2$lod, tet2_chip_repeats_lod2$AF_4)
tet2_chip_repeats_lod2$AF_6 <- ifelse(tet2_chip_repeats_lod2$AF_6 == 0.01 & !is.na(tet2_chip_repeats_lod2$lod) & tet2_chip_repeats_lod2$lod > 0.01, tet2_chip_repeats_lod2$lod, tet2_chip_repeats_lod2$AF_6)
tet2_chip_repeats_lod2$AF_8 <- ifelse(tet2_chip_repeats_lod2$AF_8 == 0.01 & !is.na(tet2_chip_repeats_lod2$lod) & tet2_chip_repeats_lod2$lod > 0.01, tet2_chip_repeats_lod2$lod, tet2_chip_repeats_lod2$AF_8)
tet2_chip_repeats_lod2$AF_9 <- ifelse(tet2_chip_repeats_lod2$AF_9 == 0.01 & !is.na(tet2_chip_repeats_lod2$lod) & tet2_chip_repeats_lod2$lod > 0.01, tet2_chip_repeats_lod2$lod, tet2_chip_repeats_lod2$AF_9)

#remove duplicates
tet2_chip_repeats_lod2 <- unique(tet2_chip_repeats_lod2)

#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with tet2 chip mutation data
tet2_chip_repeats_comorbs_drugs_lod2 <- merge(tet2_chip_repeats_lod2, w6_data_simp, by = "idauniq", all.x = TRUE)

#ensure all cases complete with medication info
tet2_chip_repeats_comorbs_drugs_lod_comp2 <- tet2_chip_repeats_comorbs_drugs_lod2 %>% dplyr::filter(!is.na(lld))

#merge mutation data with myeloid lymphoid ratio data
colnames(fbc_data_simp_w)[1] <- "idauniq"

tet2_chip_repeats_ml2 <- merge(tet2_chip_repeats_comorbs_drugs_lod_comp2, fbc_data_simp_w, by = "idauniq", all.x = TRUE)

#correct allele frequencies for m/l ratio
tet2_chip_repeats_ml2$AF_2c <- tet2_chip_repeats_ml2$AF_2/tet2_chip_repeats_ml2$m_l_ratio.2
tet2_chip_repeats_ml2$AF_4c <- tet2_chip_repeats_ml2$AF_4/tet2_chip_repeats_ml2$m_l_ratio.4
tet2_chip_repeats_ml2$AF_6c <- tet2_chip_repeats_ml2$AF_6/tet2_chip_repeats_ml2$m_l_ratio.6
tet2_chip_repeats_ml2$AF_8c <- tet2_chip_repeats_ml2$AF_8/tet2_chip_repeats_ml2$m_l_ratio.8
tet2_chip_repeats_ml2$AF_9c <- tet2_chip_repeats_ml2$AF_9/tet2_chip_repeats_ml2$m_l_ratio.9

#calculate growth rate- method 1
tet2_chip_repeats_ml2$i_1 <- ((tet2_chip_repeats_ml2$AF_9c / tet2_chip_repeats_ml2$AF_2c)^{1/15} - 1)*100
tet2_chip_repeats_ml2$i_2 <- ((tet2_chip_repeats_ml2$AF_8c / tet2_chip_repeats_ml2$AF_2c)^{1/13} - 1)*100
tet2_chip_repeats_ml2$i_3 <- ((tet2_chip_repeats_ml2$AF_9c / tet2_chip_repeats_ml2$AF_4c)^{1/11} - 1)*100
tet2_chip_repeats_ml2$i_4 <- ((tet2_chip_repeats_ml2$AF_8c / tet2_chip_repeats_ml2$AF_4c)^{1/9} - 1)*100
tet2_chip_repeats_ml2$i_5 <- ((tet2_chip_repeats_ml2$AF_8c / tet2_chip_repeats_ml2$AF_6c)^{1/5} - 1)*100
tet2_chip_repeats_ml2$i_6 <- ((tet2_chip_repeats_ml2$AF_6c / tet2_chip_repeats_ml2$AF_2c)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
tet2_chip_repeats_ml2$growth_rate <- tet2_chip_repeats_ml2$i_1
tet2_chip_repeats_ml2$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml2$growth_rate), tet2_chip_repeats_ml2$i_2, tet2_chip_repeats_ml2$growth_rate)
tet2_chip_repeats_ml2$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml2$growth_rate), tet2_chip_repeats_ml2$i_3, tet2_chip_repeats_ml2$growth_rate)
tet2_chip_repeats_ml2$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml2$growth_rate), tet2_chip_repeats_ml2$i_4, tet2_chip_repeats_ml2$growth_rate)
tet2_chip_repeats_ml2$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml2$growth_rate), tet2_chip_repeats_ml2$i_5, tet2_chip_repeats_ml2$growth_rate)
tet2_chip_repeats_ml2$growth_rate <- ifelse(is.na(tet2_chip_repeats_ml2$growth_rate), tet2_chip_repeats_ml2$i_6, tet2_chip_repeats_ml2$growth_rate)


#for independent observations select single variant - keep largest growth rate
tet2_chip_repeats_ordered2 <- tet2_chip_repeats_ml2[order(tet2_chip_repeats_ml2$idauniq, abs(tet2_chip_repeats_ml2$growth_rate), decreasing = TRUE),]
tet2_chip_repeats_single_0.01 <- tet2_chip_repeats_ordered2[ !duplicated(tet2_chip_repeats_ordered2$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
tet2_chip_repeats_single_0.01 <- merge(tet2_chip_repeats_single_0.01, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
tet2_chip_repeats_single_0.01$indager <- ifelse(tet2_chip_repeats_single_0.01$indager < 0, 90, tet2_chip_repeats_single_0.01$indager)

#count how many missing cholesterol measurement - 13/105 (12%)
nrow(tet2_chip_repeats_single_0.01[tet2_chip_repeats_single_0.01$chol.y <0, ])

#replate missing chol data with na
tet2_chip_repeats_single_0.01$chol <- ifelse(tet2_chip_repeats_single_0.01$chol.y < 0, NA, tet2_chip_repeats_single_0.01$chol.y)


#combine ihd, heart_dise and stroke into cvd
tet2_chip_repeats_single_0.01$cvd <- rep(0, nrow(tet2_chip_repeats_single_0.01))
tet2_chip_repeats_single_0.01 <- tet2_chip_repeats_single_0.01 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
tet2_chip_repeats_single_0.01 <- tet2_chip_repeats_single_0.01 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
tet2_chip_repeats_single_0.01$initial_vaf <- tet2_chip_repeats_single_0.01$AF_2c
tet2_chip_repeats_single_0.01$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single_0.01$initial_vaf), tet2_chip_repeats_single_0.01$AF_4c, tet2_chip_repeats_single_0.01$initial_vaf)
tet2_chip_repeats_single_0.01$initial_vaf <- ifelse(is.na(tet2_chip_repeats_single_0.01$initial_vaf), tet2_chip_repeats_single_0.01$AF_6c, tet2_chip_repeats_single_0.01$initial_vaf)


#simplify initial vaf as will have outliers - could take log
tet2_chip_repeats_single_0.01$log_initial_vaf <- log(tet2_chip_repeats_single_0.01$initial_vaf)


#simplify to just cvd relevant variables for robust regression
tet2_chip_repeats_cvd_0.01 <- subset(tet2_chip_repeats_single_0.01, select = c("idauniq", "growth_rate", "indager", "indsex", "cvd", "lld", "bb.y", "chol", "htn", "antiplatelets", "RAAS", "current_smoker", "log_initial_vaf", "count"))


#assess missing data
md.pattern(tet2_chip_repeats_cvd_0.01)


#drop those without calculable growth rate 
tet2_chip_repeats_cvd_0.01 <- tet2_chip_repeats_cvd_0.01 %>% dplyr::filter(!is.na(growth_rate))

#drop those without chol 
tet2_chip_repeats_cvd_0.01 <- tet2_chip_repeats_cvd_0.01 %>% dplyr::filter(!is.na(chol))


#run robust regression for growth rate method 1
dat <- tet2_chip_repeats_cvd_0.01  


form <- growth_rate ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr1_results_perc_0.01 <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr1_results_perc_0.01



#method 2 - 0.01 sensitivity analysis

tet2_chip_repeats_cvd2_0.01 <- subset(tet2_chip_repeats_single_0.01, select = c(idauniq, AF_2c, AF_4c, AF_6c, AF_8c, AF_9c, indager, indsex, cvd, bb.y, lld, chol, htn, antiplatelets, log_initial_vaf, count))

#log(VAFfollowup - VAFbaseline)/(years between measurements)

#calculate growth rate
tet2_chip_repeats_cvd2_0.01$j_1 <- log(tet2_chip_repeats_cvd2_0.01$AF_9c/tet2_chip_repeats_cvd2_0.01$AF_2c)/15
tet2_chip_repeats_cvd2_0.01$j_2 <- log(tet2_chip_repeats_cvd2_0.01$AF_8c / tet2_chip_repeats_cvd2_0.01$AF_2c)/13
tet2_chip_repeats_cvd2_0.01$j_3 <- log(tet2_chip_repeats_cvd2_0.01$AF_9c / tet2_chip_repeats_cvd2_0.01$AF_4c)/11
tet2_chip_repeats_cvd2_0.01$j_4 <- log(tet2_chip_repeats_cvd2_0.01$AF_8c / tet2_chip_repeats_cvd2_0.01$AF_4c)/9
tet2_chip_repeats_cvd2_0.01$j_5 <- log(tet2_chip_repeats_cvd2_0.01$AF_8c / tet2_chip_repeats_cvd2_0.01$AF_6c)/5
tet2_chip_repeats_cvd2_0.01$j_6 <- log(tet2_chip_repeats_cvd2_0.01$AF_6c / tet2_chip_repeats_cvd2_0.01$AF_2c)/8 


#coalesce delta vaf per year into one column
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- tet2_chip_repeats_cvd2_0.01$j_1
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_cvd2_0.01$growth_rate2), tet2_chip_repeats_cvd2_0.01$j_2, tet2_chip_repeats_cvd2_0.01$growth_rate2)
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_cvd2_0.01$growth_rate2), tet2_chip_repeats_cvd2_0.01$j_3, tet2_chip_repeats_cvd2_0.01$growth_rate2)
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_cvd2_0.01$growth_rate2), tet2_chip_repeats_cvd2_0.01$j_4, tet2_chip_repeats_cvd2_0.01$growth_rate2)
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_cvd2_0.01$growth_rate2), tet2_chip_repeats_cvd2_0.01$j_5, tet2_chip_repeats_cvd2_0.01$growth_rate2)
tet2_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(tet2_chip_repeats_cvd2_0.01$growth_rate2), tet2_chip_repeats_cvd2_0.01$j_6, tet2_chip_repeats_cvd2_0.01$growth_rate2)

#drop if growth rate not calculable
tet2_chip_repeats_cvd2_0.01 <- tet2_chip_repeats_cvd2_0.01 %>% dplyr::filter(!is.na(growth_rate2))

#drop if chol missing
tet2_chip_repeats_cvd2_0.01 <- tet2_chip_repeats_cvd2_0.01 %>% dplyr::filter(!is.na(chol))

#get percentage gr
tet2_chip_repeats_cvd2_0.01$growth_rate2_percentage <- tet2_chip_repeats_cvd2_0.01$growth_rate2*100


#run robust regression for growth rate method 2
dat <- tet2_chip_repeats_cvd2_0.01  


form <- growth_rate2_percentage ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
tet2_gr2_results_perc_0.01 <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


tet2_gr2_results_perc_0.01


#dnmt3a

lod <- subset(dnmt3a_chip_repeats_mutations_coords_depth1, select = c(idauniq, lod))

dnmt3a_chip_repeats_lod2 <- merge(dnmt3a_chip_repeats_comorbs_drugs, lod, by = "idauniq", all.x = TRUE)

dnmt3a_chip_repeats_lod2 <- dnmt3a_chip_repeats_lod2 %>% dplyr::filter(is.na(lod) | lod < 0.02)

#now replace 0.01 VAFs with lod
dnmt3a_chip_repeats_lod2$AF_2 <- ifelse(dnmt3a_chip_repeats_lod2$AF_2 == 0.01 & !is.na(dnmt3a_chip_repeats_lod2$lod) & dnmt3a_chip_repeats_lod2$lod > 0.01, dnmt3a_chip_repeats_lod2$lod, dnmt3a_chip_repeats_lod2$AF_2)
dnmt3a_chip_repeats_lod2$AF_4 <- ifelse(dnmt3a_chip_repeats_lod2$AF_4 == 0.01 & !is.na(dnmt3a_chip_repeats_lod2$lod) & dnmt3a_chip_repeats_lod2$lod > 0.01, dnmt3a_chip_repeats_lod2$lod, dnmt3a_chip_repeats_lod2$AF_4)
dnmt3a_chip_repeats_lod2$AF_6 <- ifelse(dnmt3a_chip_repeats_lod2$AF_6 == 0.01 & !is.na(dnmt3a_chip_repeats_lod2$lod) & dnmt3a_chip_repeats_lod2$lod > 0.01, dnmt3a_chip_repeats_lod2$lod, dnmt3a_chip_repeats_lod2$AF_6)
dnmt3a_chip_repeats_lod2$AF_8 <- ifelse(dnmt3a_chip_repeats_lod2$AF_8 == 0.01 & !is.na(dnmt3a_chip_repeats_lod2$lod) & dnmt3a_chip_repeats_lod2$lod > 0.01, dnmt3a_chip_repeats_lod2$lod, dnmt3a_chip_repeats_lod2$AF_8)
dnmt3a_chip_repeats_lod2$AF_9 <- ifelse(dnmt3a_chip_repeats_lod2$AF_9 == 0.01 & !is.na(dnmt3a_chip_repeats_lod2$lod) & dnmt3a_chip_repeats_lod2$lod > 0.01, dnmt3a_chip_repeats_lod2$lod, dnmt3a_chip_repeats_lod2$AF_9)

#remove duplicates
dnmt3a_chip_repeats_lod2 <- unique(dnmt3a_chip_repeats_lod2)

#wave 6 interview data - HeChMd, HeChMe
setwd("N:/My Documents/ELSA_data/elsa_survey_data/")
w6_data <- read.table("wave_6_clean.txt", header = TRUE, sep = "\t")
#simplify
w6_data_simp <- subset(w6_data, select = c(idauniq, HeChMd, HeChMe))

#merge these data with dnmt3a chip mutation data
dnmt3a_chip_repeats_comorbs_drugs_lod2 <- merge(dnmt3a_chip_repeats_lod2, w6_data_simp, by = "idauniq", all.x = TRUE)

#ensure all cases complete with medication info
dnmt3a_chip_repeats_comorbs_drugs_lod_comp2 <- dnmt3a_chip_repeats_comorbs_drugs_lod2 %>% dplyr::filter(!is.na(lld))

#merge mutation data with myeloid lymphoid ratio data
colnames(fbc_data_simp_w)[1] <- "idauniq"

dnmt3a_chip_repeats_ml2 <- merge(dnmt3a_chip_repeats_comorbs_drugs_lod_comp2, fbc_data_simp_w, by = "idauniq", all.x = TRUE)

#correct allele frequencies for m/l ratio
dnmt3a_chip_repeats_ml2$AF_2c <- dnmt3a_chip_repeats_ml2$AF_2/dnmt3a_chip_repeats_ml2$m_l_ratio.2
dnmt3a_chip_repeats_ml2$AF_4c <- dnmt3a_chip_repeats_ml2$AF_4/dnmt3a_chip_repeats_ml2$m_l_ratio.4
dnmt3a_chip_repeats_ml2$AF_6c <- dnmt3a_chip_repeats_ml2$AF_6/dnmt3a_chip_repeats_ml2$m_l_ratio.6
dnmt3a_chip_repeats_ml2$AF_8c <- dnmt3a_chip_repeats_ml2$AF_8/dnmt3a_chip_repeats_ml2$m_l_ratio.8
dnmt3a_chip_repeats_ml2$AF_9c <- dnmt3a_chip_repeats_ml2$AF_9/dnmt3a_chip_repeats_ml2$m_l_ratio.9

#calculate growth rate- method 1
dnmt3a_chip_repeats_ml2$i_1 <- ((dnmt3a_chip_repeats_ml2$AF_9c / dnmt3a_chip_repeats_ml2$AF_2c)^{1/15} - 1)*100
dnmt3a_chip_repeats_ml2$i_2 <- ((dnmt3a_chip_repeats_ml2$AF_8c / dnmt3a_chip_repeats_ml2$AF_2c)^{1/13} - 1)*100
dnmt3a_chip_repeats_ml2$i_3 <- ((dnmt3a_chip_repeats_ml2$AF_9c / dnmt3a_chip_repeats_ml2$AF_4c)^{1/11} - 1)*100
dnmt3a_chip_repeats_ml2$i_4 <- ((dnmt3a_chip_repeats_ml2$AF_8c / dnmt3a_chip_repeats_ml2$AF_4c)^{1/9} - 1)*100
dnmt3a_chip_repeats_ml2$i_5 <- ((dnmt3a_chip_repeats_ml2$AF_8c / dnmt3a_chip_repeats_ml2$AF_6c)^{1/5} - 1)*100
dnmt3a_chip_repeats_ml2$i_6 <- ((dnmt3a_chip_repeats_ml2$AF_6c / dnmt3a_chip_repeats_ml2$AF_2c)^{1/8} - 1)*100 


#coalesce delta vaf per year into one column
dnmt3a_chip_repeats_ml2$growth_rate <- dnmt3a_chip_repeats_ml2$i_1
dnmt3a_chip_repeats_ml2$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml2$growth_rate), dnmt3a_chip_repeats_ml2$i_2, dnmt3a_chip_repeats_ml2$growth_rate)
dnmt3a_chip_repeats_ml2$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml2$growth_rate), dnmt3a_chip_repeats_ml2$i_3, dnmt3a_chip_repeats_ml2$growth_rate)
dnmt3a_chip_repeats_ml2$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml2$growth_rate), dnmt3a_chip_repeats_ml2$i_4, dnmt3a_chip_repeats_ml2$growth_rate)
dnmt3a_chip_repeats_ml2$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml2$growth_rate), dnmt3a_chip_repeats_ml2$i_5, dnmt3a_chip_repeats_ml2$growth_rate)
dnmt3a_chip_repeats_ml2$growth_rate <- ifelse(is.na(dnmt3a_chip_repeats_ml2$growth_rate), dnmt3a_chip_repeats_ml2$i_6, dnmt3a_chip_repeats_ml2$growth_rate)


#for independent observations select single variant - keep largest growth rate
dnmt3a_chip_repeats_ordered2 <- dnmt3a_chip_repeats_ml2[order(dnmt3a_chip_repeats_ml2$idauniq, abs(dnmt3a_chip_repeats_ml2$growth_rate), decreasing = TRUE),]
dnmt3a_chip_repeats_single_0.01 <- dnmt3a_chip_repeats_ordered2[ !duplicated(dnmt3a_chip_repeats_ordered2$idauniq), ]

#add in wave 6 nurse variables - as mid point, consistent with medication data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
w6_nurse_data <- read.table("wave_6_nurse_clean.txt", sep = "\t", header = TRUE)
w6_nurse_data_simp <- subset(w6_nurse_data, select = c(idauniq, hgb, wbc, mch, hscrp, rtin, BMIVAL, chol, cfib))
dnmt3a_chip_repeats_single_0.01 <- merge(dnmt3a_chip_repeats_single_0.01, w6_nurse_data_simp, by = "idauniq", all.x = TRUE)

#replace censored data for age (>90) with 90
dnmt3a_chip_repeats_single_0.01$indager <- ifelse(dnmt3a_chip_repeats_single_0.01$indager < 0, 90, dnmt3a_chip_repeats_single_0.01$indager)

#count how many missing cholesterol measurement - 13/105 (12%)
nrow(dnmt3a_chip_repeats_single_0.01[dnmt3a_chip_repeats_single_0.01$chol.y <0, ])

#replate missing chol data with na
dnmt3a_chip_repeats_single_0.01$chol <- ifelse(dnmt3a_chip_repeats_single_0.01$chol.y < 0, NA, dnmt3a_chip_repeats_single_0.01$chol.y)


#combine ihd, heart_dise and stroke into cvd
dnmt3a_chip_repeats_single_0.01$cvd <- rep(0, nrow(dnmt3a_chip_repeats_single_0.01))
dnmt3a_chip_repeats_single_0.01 <- dnmt3a_chip_repeats_single_0.01 %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))
dnmt3a_chip_repeats_single_0.01 <- dnmt3a_chip_repeats_single_0.01 %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#make variable for initial VAF
dnmt3a_chip_repeats_single_0.01$initial_vaf <- dnmt3a_chip_repeats_single_0.01$AF_2c
dnmt3a_chip_repeats_single_0.01$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single_0.01$initial_vaf), dnmt3a_chip_repeats_single_0.01$AF_4c, dnmt3a_chip_repeats_single_0.01$initial_vaf)
dnmt3a_chip_repeats_single_0.01$initial_vaf <- ifelse(is.na(dnmt3a_chip_repeats_single_0.01$initial_vaf), dnmt3a_chip_repeats_single_0.01$AF_6c, dnmt3a_chip_repeats_single_0.01$initial_vaf)


#simplify initial vaf as will have outliers - could take log
dnmt3a_chip_repeats_single_0.01$log_initial_vaf <- log(dnmt3a_chip_repeats_single_0.01$initial_vaf)

#simplify to just cvd relevant variables for robust regression
dnmt3a_chip_repeats_cvd_0.01 <- subset(dnmt3a_chip_repeats_single_0.01, select = c("idauniq", "growth_rate", "indager", "indsex", "cvd", "lld", "bb.y", "chol", "htn", "antiplatelets", "RAAS", "current_smoker", "log_initial_vaf", "count"))


#assess missing data
md.pattern(dnmt3a_chip_repeats_cvd_0.01)


#drop those without calculable growth rate 
dnmt3a_chip_repeats_cvd_0.01 <- dnmt3a_chip_repeats_cvd_0.01 %>% dplyr::filter(!is.na(growth_rate))

#drop those without chol 
dnmt3a_chip_repeats_cvd_0.01 <- dnmt3a_chip_repeats_cvd_0.01 %>% dplyr::filter(!is.na(chol))


#run robust regression for growth rate method 1
dat <- dnmt3a_chip_repeats_cvd_0.01  


form <- growth_rate ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr1_results_perc_0.01 <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr1_results_perc_0.01




#method 2

#method 2 - 0.01 sensitivity analysis

dnmt3a_chip_repeats_cvd2_0.01 <- subset(dnmt3a_chip_repeats_single_0.01, select = c(idauniq, AF_2c, AF_4c, AF_6c, AF_8c, AF_9c, indager, indsex, cvd, bb.y, lld, chol, htn, antiplatelets, log_initial_vaf, count))

#log(VAFfollowup - VAFbaseline)/(years between measurements)

#calculate growth rate
dnmt3a_chip_repeats_cvd2_0.01$j_1 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_9c/dnmt3a_chip_repeats_cvd2_0.01$AF_2c)/15
dnmt3a_chip_repeats_cvd2_0.01$j_2 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_8c / dnmt3a_chip_repeats_cvd2_0.01$AF_2c)/13
dnmt3a_chip_repeats_cvd2_0.01$j_3 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_9c / dnmt3a_chip_repeats_cvd2_0.01$AF_4c)/11
dnmt3a_chip_repeats_cvd2_0.01$j_4 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_8c / dnmt3a_chip_repeats_cvd2_0.01$AF_4c)/9
dnmt3a_chip_repeats_cvd2_0.01$j_5 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_8c / dnmt3a_chip_repeats_cvd2_0.01$AF_6c)/5
dnmt3a_chip_repeats_cvd2_0.01$j_6 <- log(dnmt3a_chip_repeats_cvd2_0.01$AF_6c / dnmt3a_chip_repeats_cvd2_0.01$AF_2c)/8 


#coalesce delta vaf per year into one column
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- dnmt3a_chip_repeats_cvd2_0.01$j_1
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_cvd2_0.01$growth_rate2), dnmt3a_chip_repeats_cvd2_0.01$j_2, dnmt3a_chip_repeats_cvd2_0.01$growth_rate2)
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_cvd2_0.01$growth_rate2), dnmt3a_chip_repeats_cvd2_0.01$j_3, dnmt3a_chip_repeats_cvd2_0.01$growth_rate2)
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_cvd2_0.01$growth_rate2), dnmt3a_chip_repeats_cvd2_0.01$j_4, dnmt3a_chip_repeats_cvd2_0.01$growth_rate2)
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_cvd2_0.01$growth_rate2), dnmt3a_chip_repeats_cvd2_0.01$j_5, dnmt3a_chip_repeats_cvd2_0.01$growth_rate2)
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2 <- ifelse(is.na(dnmt3a_chip_repeats_cvd2_0.01$growth_rate2), dnmt3a_chip_repeats_cvd2_0.01$j_6, dnmt3a_chip_repeats_cvd2_0.01$growth_rate2)

#drop if growth rate not calculable
dnmt3a_chip_repeats_cvd2_0.01 <- dnmt3a_chip_repeats_cvd2_0.01 %>% dplyr::filter(!is.na(growth_rate2))

#drop if chol missing
dnmt3a_chip_repeats_cvd2_0.01 <- dnmt3a_chip_repeats_cvd2_0.01 %>% dplyr::filter(!is.na(chol))

#get percentage gr
dnmt3a_chip_repeats_cvd2_0.01$growth_rate2_percentage <- dnmt3a_chip_repeats_cvd2_0.01$growth_rate2*100


#run robust regression for growth rate method 2
dat <- dnmt3a_chip_repeats_cvd2_0.01  


form <- growth_rate2_percentage ~ indager + indsex + cvd + lld + chol + htn + log_initial_vaf + count

#Fit once to get baseline estimates 
fit0 <- rlm(form, data = dat, maxit = 200)
coef_names <- names(coef(fit0))       
p <- length(coef_names)

#Bootstrap function
boot_fun <- function(d, i) {
  dd  <- d[i, , drop = FALSE]
  fit <- try(suppressWarnings(rlm(form, data = dd, maxit = 200)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, p))
  cc <- coef(fit)
  out <- rep(NA_real_, p); names(out) <- coef_names
  out[names(cc)] <- cc
  out
}



set.seed(161)
B <- 1000  # number of bootstrap resamples
bt <- boot(dat, statistic = boot_fun, R = B)

#Percentile-tail p-value 
p_emp_two_sided <- function(v, theta0 = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  2 * min(mean(v <= theta0), mean(v >= theta0))
}


#Build results
dnmt3a_gr2_results_perc_0.01 <- do.call(rbind, lapply(seq_along(coef_names), function(j) {
  bj <- bt$t[, j]
  # Percentile CI via boot.ci (fallback to simple quantiles if needed)
  ci_obj  <- try(boot.ci(bt, index = j, type = "perc"), silent = TRUE)
  if (!inherits(ci_obj, "try-error") && !is.null(ci_obj$percent)) {
    ci_low  <- ci_obj$percent[4]
    ci_high <- ci_obj$percent[5]
  } else {
    qs <- quantile(bj[is.finite(bj)], probs = c(0.025, 0.975), names = FALSE)
    ci_low <- qs[1]; ci_high <- qs[2]
  }
  
  data.frame(
    term      = coef_names[j],
    estimate  = unname(coef(fit0)[j]),
    se_boot   = sd(bj, na.rm = TRUE),
    ci_low    = ci_low,
    ci_high   = ci_high,
    p_emp     = p_emp_two_sided(bj, theta0 = 0),
    n_eff     = sum(is.finite(bj)),
    row.names = NULL
  )
}))


dnmt3a_gr2_results_perc_0.01


