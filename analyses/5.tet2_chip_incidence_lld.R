#script to look at whether statins have impact on incidence of tet2 chip
library(dplyr)
library(readxl)
library(tidyr)
library(survival)
library(survminer)

#load growth rate data 
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions/growth_rates")
load("growth_rate_calcs_data.RData")

#take those with repeated sampling for tet2 mutations where lod < 0.02

incident_tet2 <- tet2_chip_repeats_comorbs_drugs_lod %>% dplyr::filter(lod < 0.02)

#simplify

incident_tet2_simp <- subset(incident_tet2, select = c(idauniq, AF_2, AF_4, AF_6, AF_8, AF_9))

#one individual loses detectable tet2 clone - remove. 
incident_tet2_simp <- incident_tet2_simp %>% dplyr::filter(idauniq != "[censored]")

#add time variable for wave of incident chip

incident_tet2_simp$time <- NA
incident_tet2_simp$time <- ifelse(incident_tet2_simp$AF_4 >= 0.02 & !is.na(incident_tet2_simp$AF_4), 4, incident_tet2_simp$time)
incident_tet2_simp$time <- ifelse(incident_tet2_simp$AF_6 >= 0.02 & !is.na(incident_tet2_simp$AF_6), 6, incident_tet2_simp$time)
#one individual had incident chip at wave 6
incident_tet2_simp$time <- ifelse(incident_tet2_simp$AF_8 >= 0.02 & !is.na(incident_tet2_simp$AF_8) & is.na(incident_tet2_simp$time), 8, incident_tet2_simp$time)
incident_tet2_simp$time <- ifelse(incident_tet2_simp$AF_9 >= 0.02 & !is.na(incident_tet2_simp$AF_9) & is.na(incident_tet2_simp$time), 9, incident_tet2_simp$time)

#add outcome
incident_tet2_simp$outcome <- "tet2_chip"

#simplify just to idauniq, time, outcome
incident_tet2_simp2 <- subset(incident_tet2_simp, select = c(idauniq, time, outcome))

#remove duplicates 
incident_tet2_simp2 <- unique(incident_tet2_simp2)

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

#combine tet2 chip cases and controls
repeats_for_coxph <- rbind(incident_tet2_simp2, control_repeats_simp)

#load cumulative incidence data of main comorbidities and statin treatment
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions/cvd_cumulative_incidence")
load("revisions_cvd_ci.RData")

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

#2 individuals have NA at wave 6 for statins but have never been prescribed statin: 106668, 161328
statin_id <- as.data.frame(statin_id)
statin_id[statin_id$ids == "106668",][2] <- 0
statin_id[statin_id$ids == "161328",][2] <- 0
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

#two individuals again missing wave 6 - 106668 and 161328
#106668 - 64 and sex 2
#161328 - 64 and sex 2
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "106668",]$indager <- 64
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "106668",]$indsex <- 2

repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "161328",]$indager <- 64
repeats_for_coxph_fw_statins_age[repeats_for_coxph_fw_statins_age$idauniq == "161328",]$indsex <- 2

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
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 6, 2012, NA)
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 8, 2017, repeats_for_coxph_fw_statins_age$time_year)
repeats_for_coxph_fw_statins_age$time_year <- ifelse(repeats_for_coxph_fw_statins_age$time == 9, 2019, repeats_for_coxph_fw_statins_age$time_year)

#time in years between start and censoring
repeats_for_coxph_fw_statins_age$time_in_years <- repeats_for_coxph_fw_statins_age$time_year - repeats_for_coxph_fw_statins_age$year_1

#check for baseline differences - age, sex, cholesterol, smoking, diabetes, htn
#age - no diff
tet2_incidence_data <- repeats_for_coxph_fw_statins_age

tet2_incidence_cases <- tet2_incidence_data %>% dplyr::filter(outcome == "tet2_chip")
tet2_incidence_controls <- tet2_incidence_data %>% dplyr::filter(outcome == "control")

t.test(tet2_incidence_cases$indager, tet2_incidence_controls$indager)
sd(tet2_incidence_cases$indager)
sd(tet2_incidence_controls$indager)

#sex 
sex.m <- as.matrix(table(tet2_incidence_data$outcome, tet2_incidence_data$indsex))
chisq.test(sex.m)

#cholesterol 
t.test(tet2_incidence_cases$chol, tet2_incidence_controls$chol)
sd(tet2_incidence_cases$chol)
sd(tet2_incidence_controls$chol)

#smoking 
smoke.m <- as.matrix(table(tet2_incidence_data$outcome, tet2_incidence_data$smoke))
chisq.test(tet2_incidence_data$outcome, tet2_incidence_data$smoke)

#diabetes 
dm.m <- as.matrix(table(tet2_incidence_data$outcome, tet2_incidence_data$diabetes))
chisq.test(tet2_incidence_data$outcome, tet2_incidence_data$diabetes)
fisher.test(tet2_incidence_data$outcome, tet2_incidence_data$diabetes)

#htn 
htn.m <- as.matrix(table(tet2_incidence_data$outcome, tet2_incidence_data$htn))
chisq.test(tet2_incidence_data$outcome, tet2_incidence_data$htn)

#cvd
tet2_incidence_data$cvd <- ifelse(tet2_incidence_data$ihd == 1 | tet2_incidence_data$stroke == 1, 1, 0)
cvd.m <- as.matrix(table(tet2_incidence_data$outcome, tet2_incidence_data$cvd))
chisq.test(tet2_incidence_data$outcome, tet2_incidence_data$cvd)
fisher.test(tet2_incidence_data$outcome, tet2_incidence_data$cvd)

#follow up
tet2_incidence_cases$follow_up <- tet2_incidence_cases$year_f - tet2_incidence_cases$year_1
tet2_incidence_controls$follow_up <- tet2_incidence_controls$year_f - tet2_incidence_controls$year_1


t.test(tet2_incidence_cases$follow_up, tet2_incidence_controls$follow_up)
sd(tet2_incidence_cases$follow_up)
sd(tet2_incidence_controls$follow_up)

#do coxph 
tet2_incidence_data$status <- ifelse(tet2_incidence_data$outcome == "tet2_chip", TRUE, FALSE)
mod_tet2_incidence <- Surv(tet2_incidence_data$time_in_years, tet2_incidence_data$status)

tet2_coxph_model <- coxph(Surv(time_in_years, status) ~ cvd + statin, data = tet2_incidence_data)
cox.zph(tet2_coxph_model)
summary(tet2_coxph_model)

tet2_coxph_model2 <- coxph(Surv(time_in_years, status) ~ cvd + strata(statin), data = tet2_incidence_data)
summary(tet2_coxph_model2)

fit <- survfit(formula = mod_tet2_incidence ~ cvd + strata(statin), data = tet2_incidence_data)
fit2 <- survfit(tet2_coxph_model2)

ggsurvplot(fit2, data = tet2_incidence_data, conf.int = TRUE, risk.table = TRUE, risk.table.col = "strata", fun = "event",
           break.time.by = 1, 
           ylab = "TET2 CHIP cumulative events",
           xlab = "Time (years)",
           legend.labs = c("No statin", "Statin"))


