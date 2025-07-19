#script to repeat statin primary prevention cumulative incidence 

library(readxl)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(mice)
library(ggmice)
library(MatchIt)
library(adjustedCurves)

setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions")

chip_vars <- read.table("anon_chip_variant_data_clean_asxl1_fixed.txt", header = TRUE, sep = "\t")

#need to remove TET2 Q1274L - check whether any of these were single mutations and therefore whether they should be removed from list of cases
Q1274L <- chip_vars %>% dplyr::filter(NonsynOI == "Q1274L")
no_Q1274L <- chip_vars %>% dplyr::filter(NonsynOI != "Q1274L")

#simplify to get sample number
Q1274L_simp <- subset(Q1274L, select = c(Sample, NonsynOI))
Q1274L_simp$other_vars <- ifelse(Q1274L_simp$Sample %in% no_Q1274L$Sample, TRUE, FALSE)


setwd("N:/My Documents/ELSA_data/")
cases <- readxl::read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Chip cases.xlsx")

#if no other vars detected for that sample then remove from cases (27 out of 35 cases)
Q1274L_only <- Q1274L_simp %>% dplyr::filter(Q1274L_simp$other_vars == FALSE)

cases_noQ1274L <- cases %>% dplyr::filter(!cases$Sample %in% Q1274L_only$Sample)

#import elsa survey data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")


#survey data
wave_1 <- read.csv("wave_1_core_data_v3.tab", header = TRUE, sep = "\t")
wave_2 <- read.csv("wave_2_clean.txt", header = TRUE, sep = "\t")
wave_4 <- read.csv("wave_4_clean.txt", header = TRUE, sep = "\t")
wave_6 <- read.csv("wave_6_clean.txt", header = TRUE, sep = "\t")
wave_8 <- read.csv("wave_8_clean_new.txt", header = TRUE, sep = "\t")
wave_9 <- read.csv("wave_9_clean.txt", header = TRUE, sep = "\t")

#extra waves
wave_3 <- read.csv("wave_3_elsa_data_v4.tab", header = TRUE, sep = "\t")
wave_5 <- read.csv("wave_5_elsa_data_v4.tab", header = TRUE, sep = "\t")
wave_7 <- read.csv("wave_7_elsa_data.tab", header = TRUE, sep = "\t")

#wave 6 and 8/9 medications
setwd("N:/My Documents/ELSA_data/drug_regression_analysis")
wave_6_meds <- read.csv("wave_6_drug_categories.txt", header = TRUE, sep = "\t")

#nurse wave data
setwd("N:/My Documents/ELSA_data/elsa_survey_data")
wave_2_nurse <- read.csv("wave_2_nurse_clean.txt", header = TRUE, sep = "\t")
wave_4_nurse <- read.csv("wave_4_nurse_clean.txt", header = TRUE, sep = "\t")
wave_6_nurse <- read.csv("wave_6_nurse_clean.txt", header = TRUE, sep = "\t")
wave_8_nurse <- read.csv("wave_8_nurse_clean.txt", header = TRUE, sep = "\t")
wave_8_9_nurse <- read.csv("elsa_nurse_w8w9_data_sl_protect.tab", header = TRUE, sep = "\t")

#keep only data corresponding to cases and controls

#import controls
setwd("N:/My Documents/ELSA_data")
controls <- read_xlsx("Waves 2,4,6,8,9 Lab IDs and Idauniq -Control cases.xlsx")

#now drop unnecessary variables
controls_simp <- subset(controls, select=-c(bldrec, consn, CONBST, ConStorB))

#####CASES
#now drop unnecessary variables

cases_simp <- subset(cases_noQ1274L, select=-c(bldrec, consn, CONBST))

#drop observations where idauniq is NA

cases_simp <- drop_na(cases_simp, idauniq)

#create label
case_label <- c(rep("case", length(cases_simp$idauniq)))
cases_simp$label <- case_label

control_label <- c(rep("ctl", length(controls_simp$idauniq)))
controls_simp$label <- control_label

#create df of elsa survey data from just cases
wave_2_cases <- wave_2[(wave_2$idauniq %in% cases_simp$idauniq),]
wave_4_cases <- wave_4[(wave_4$idauniq %in% cases_simp$idauniq),]
wave_6_cases <- wave_6[(wave_6$idauniq %in% cases_simp$idauniq),]
wave_8_cases <- wave_8[(wave_8$idauniq %in% cases_simp$idauniq),]
wave_9_cases <- wave_9[(wave_9$idauniq %in% cases_simp$idauniq),]


#dfs of elsa survey data from just controls
wave_2_controls <- wave_2[(wave_2$idauniq %in% controls_simp$idauniq),]
wave_4_controls <- wave_4[(wave_4$idauniq %in% controls_simp$idauniq),]
wave_6_controls <- wave_6[(wave_6$idauniq %in% controls_simp$idauniq),]
wave_8_controls <- wave_8[(wave_8$idauniq %in% controls_simp$idauniq),]
wave_9_controls <- wave_9[(wave_9$idauniq %in% controls_simp$idauniq),]


#create df to allow ascertainment of incident cv disease across the cohort by wave
#ihd - just MI
#wave 1 - ever been told had a heart attack
wave_1_summary <- wave_1
wave_1_summary$ihd1 <- 0
wave_1_summary$ihd1 <- ifelse(wave_1_summary$hedia01 == 3 |  wave_1_summary$hedia02 == 3 |  wave_1_summary$hedia03 == 3 |  wave_1_summary$hedia04 == 3 | wave_1_summary$hedia05 == 3 | wave_1_summary$hedia06 == 3 | wave_1_summary$hedia07 == 3 | wave_1_summary$hedia08 == 3 | wave_1_summary$hedia09 == 3 | wave_1_summary$hedia10 == 3, 1, 0)
wave_1_simp <- subset(wave_1_summary, select = c(idauniq, ihd1))



#wave 2 - newly diagnosed
wave_2_summary <- wave_2
wave_2_summary$ihd2 <- 0
wave_2_summary$ihd2 <- ifelse( wave_2_summary$hedia01 == 3 |  wave_2_summary$hedia02 == 3 |  wave_2_summary$hedia03 == 3 | wave_2_summary$hedia04 == 3 | wave_2_summary$hedia05 == 3 | wave_2_summary$hedia06 == 3 | wave_2_summary$hedia07 == 3 | wave_2_summary$hedia08 == 3 | wave_2_summary$hedia09 == 3, 1, 0)
wave_2_simp <- subset(wave_2_summary, select = c(idauniq, ihd2))

#wave 3
wave_3_summary <- wave_3
wave_3_summary$ihd3 <- 0
wave_3_summary$ihd3 <- ifelse(wave_3_summary$hediami == 1, 1, 0)
wave_3_simp <- subset(wave_3_summary, select = c(idauniq, ihd3))


#wave 4

wave_4_summary <- wave_4
wave_4_summary$ihd4 <- 0
wave_4_summary$ihd4 <- ifelse(wave_4_summary$hediami == 1, 1, 0)
wave_4_simp <- subset(wave_4_summary, select = c(idauniq, ihd4))

#wave 5
wave_5_summary <- wave_5
wave_5_summary$ihd5 <- 0
wave_5_summary$ihd5 <- ifelse(wave_5_summary$hediami == 1, 1, 0)
wave_5_simp <- subset(wave_5_summary, select = c(idauniq, ihd5))

#wave 6
wave_6_summary <- wave_6
wave_6_summary$ihd6 <- 0
wave_6_summary$ihd6 <- ifelse(wave_6_summary$hediami == 1, 1, 0)
wave_6_simp <- subset(wave_6_summary, select = c(idauniq, ihd6))

#wave 7
wave_7_summary <- wave_7
wave_7_summary$ihd7 <- 0
wave_7_summary$ihd7 <- ifelse(wave_7_summary$hediami == 1, 1, 0)
wave_7_simp <- subset(wave_7_summary, select = c(idauniq, ihd7))

#wave 8
wave_8_summary <- wave_8
wave_8_summary$ihd8 <- 0
wave_8_summary$ihd8 <- ifelse(wave_8_summary$hediami == 1 , 1, 0)
wave_8_simp <- subset(wave_8_summary, select = c(idauniq, ihd8))

#wave 9
wave_9_summary <- wave_9
wave_9_summary$ihd9 <- 0
wave_9_summary$ihd9 <- ifelse(wave_9_summary$hediami == 1 , 1, 0)
wave_9_simp <- subset(wave_9_summary, select = c(idauniq, ihd9))

#combine ihd
ihd_longitudinal <- merge(wave_1_simp, wave_2_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_3_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_4_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_5_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_6_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_7_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_8_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
ihd_longitudinal <- merge(ihd_longitudinal, wave_9_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)



#stroke
#wave 1
wave_1_summary$stroke1 <- 0
wave_1_summary$stroke1 <- ifelse(wave_1_summary$hedia01 == 8 | wave_1_summary$hedia02 == 8 | wave_1_summary$hedia03 == 8 | wave_1_summary$hedia04 == 8 | wave_1_summary$hedia05 == 8 | wave_1_summary$hedia06 == 8 | wave_1_summary$hedia07 == 8 | wave_1_summary$hedia08 == 8 | wave_1_summary$hedia09 == 8 | wave_1_summary$hedia10 == 8, 1, 0)
wave_1_stroke_simp <- subset(wave_1_summary, select = c(idauniq, stroke1))


#wave 2
wave_2_summary$stroke2 <- 0
wave_2_summary$stroke2 <- ifelse(wave_2_summary$hedia01 == 8 | wave_2_summary$hedia02 == 8 | wave_2_summary$hedia03 == 8 | wave_2_summary$hedia04 == 8 | wave_2_summary$hedia05 == 8 | wave_2_summary$hedia06 == 8 | wave_2_summary$hedia07 == 8 | wave_2_summary$hedia08 == 8 | wave_2_summary$hedia09 == 8, 1, 0)
wave_2_stroke_simp <- subset(wave_2_summary, select = c(idauniq, stroke2))

#wave_3
wave_3_summary$stroke3 <- 0
wave_3_summary$stroke3 <- ifelse(wave_3_summary$hediast == 1, 1, 0)
wave_3_stroke_simp <- subset(wave_3_summary, select = c(idauniq, stroke3))


#wave_4
wave_4_summary$stroke4 <- 0
wave_4_summary$stroke4 <- ifelse(wave_4_summary$hediast == 1, 1, 0)
wave_4_stroke_simp <- subset(wave_4_summary, select = c(idauniq, stroke4))

#wave_5
wave_5_summary$stroke5 <- 0
wave_5_summary$stroke5 <- ifelse(wave_5_summary$hediast == 1, 1, 0)
wave_5_stroke_simp <- subset(wave_5_summary, select = c(idauniq, stroke5))


#wave_6
wave_6_summary$stroke6 <- 0
wave_6_summary$stroke6 <- ifelse(wave_6_summary$hediast == 1, 1, 0)
wave_6_stroke_simp <- subset(wave_6_summary, select = c(idauniq, stroke6))

#wave_7
wave_7_summary$stroke7 <- 0
wave_7_summary$stroke7 <- ifelse(wave_7_summary$hediast == 1, 1, 0)
wave_7_stroke_simp <- subset(wave_7_summary, select = c(idauniq, stroke7))


#wave_8
wave_8_summary$stroke8 <- 0
wave_8_summary$stroke8 <- ifelse(wave_8_summary$hedimst == 1, 1, 0)
wave_8_stroke_simp <- subset(wave_8_summary, select = c(idauniq, stroke8))

#wave_9
wave_9_summary$stroke9 <- 0
wave_9_summary$stroke9 <- ifelse(wave_9_summary$hediast == 1, 1, 0)
wave_9_stroke_simp <- subset(wave_9_summary, select = c(idauniq, stroke9))

#combine stroke
stroke_longitudinal <- merge(wave_1_stroke_simp, wave_2_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_3_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_4_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_5_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_6_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_7_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_8_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
stroke_longitudinal <- merge(stroke_longitudinal, wave_9_stroke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)



#htn
#wave 1
wave_1_summary$htn1 <- 0
wave_1_summary$htn1 <- ifelse( wave_1_summary$hedia01 == 1 | wave_1_summary$hedia02 == 1 | wave_1_summary$hedia03 == 1 | wave_1_summary$hedia04 == 1 | wave_1_summary$hedia05 == 1 | wave_1_summary$hedia06 == 1 | wave_1_summary$hedia07 == 1 | wave_1_summary$hedia08 == 1 | wave_1_summary$hedia09 == 1 | wave_1_summary$hedia10 == 1, 1, 0)
wave_1_htn_simp <- subset(wave_1_summary, select = c(idauniq, htn1))


#wave 2
wave_2_summary$htn2 <- 0
wave_2_summary$htn2 <- ifelse(wave_2_summary$hedia01 == 1 | wave_2_summary$hedia02 == 1 | wave_2_summary$hedia03 == 1 | wave_2_summary$hedia04 == 1 | wave_2_summary$hedia05 == 1 | wave_2_summary$hedia06 == 1 | wave_2_summary$hedia07 == 1 | wave_2_summary$hedia08 == 1 | wave_2_summary$hedia09 == 1, 1, 0)
wave_2_htn_simp <- subset(wave_2_summary, select = c(idauniq, htn2))

#wave_3
wave_3_summary$htn3 <- 0
wave_3_summary$htn3 <- ifelse(wave_3_summary$hediabp == 1, 1, 0)
wave_3_htn_simp <- subset(wave_3_summary, select = c(idauniq, htn3))


#wave_4
wave_4_summary$htn4 <- 0
wave_4_summary$htn4 <- ifelse(wave_4_summary$hediabp == 1, 1, 0)
wave_4_htn_simp <- subset(wave_4_summary, select = c(idauniq, htn4))

#wave_5
wave_5_summary$htn5 <- 0
wave_5_summary$htn5 <- ifelse(wave_5_summary$hediabp == 1, 1, 0)
wave_5_htn_simp <- subset(wave_5_summary, select = c(idauniq, htn5))


#wave_6
wave_6_summary$htn6 <- 0
wave_6_summary$htn6 <- ifelse(wave_6_summary$hediabp == 1, 1, 0)
wave_6_htn_simp <- subset(wave_6_summary, select = c(idauniq, htn6))

#wave_7
wave_7_summary$htn7 <- 0
wave_7_summary$htn7 <- ifelse(wave_7_summary$hediabp == 1, 1, 0)
wave_7_htn_simp <- subset(wave_7_summary, select = c(idauniq, htn7))


#wave_8
wave_8_summary$htn8 <- 0
wave_8_summary$htn8 <- ifelse(wave_8_summary$hediabp == 1, 1, 0)
wave_8_htn_simp <- subset(wave_8_summary, select = c(idauniq, htn8))

#wave_9
wave_9_summary$htn9 <- 0
wave_9_summary$htn9 <- ifelse(wave_9_summary$hediabp == 1, 1, 0)
wave_9_htn_simp <- subset(wave_9_summary, select = c(idauniq, htn9))

#combine htn
htn_longitudinal <- merge(wave_1_htn_simp, wave_2_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_3_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_4_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_5_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_6_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_7_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_8_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
htn_longitudinal <- merge(htn_longitudinal, wave_9_htn_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)


#diabetes
#wave 1
wave_1_summary$diabetes1 <- 0
wave_1_summary$diabetes1 <- ifelse(wave_1_summary$hedia01 == 7 | wave_1_summary$hedia02 == 7 | wave_1_summary$hedia03 == 7 | wave_1_summary$hedia04 == 7 | wave_1_summary$hedia05 == 7 | wave_1_summary$hedia06 == 7 | wave_1_summary$hedia07 == 7 | wave_1_summary$hedia08 == 7 | wave_1_summary$hedia09 == 7 | wave_1_summary$hedia10 == 7, 1, 0)
wave_1_diabetes_simp <- subset(wave_1_summary, select = c(idauniq, diabetes1))


#wave 2
wave_2_summary$diabetes2 <- 0
wave_2_summary$diabetes2 <- ifelse(wave_2_summary$HeDiaC7 == 1 | wave_2_summary$hedia01 == 7 | wave_2_summary$hedia02 == 7 | wave_2_summary$hedia03 == 7 | wave_2_summary$hedia04 == 7 | wave_2_summary$hedia05 == 7 | wave_2_summary$hedia06 == 7 | wave_2_summary$hedia07 == 7 | wave_2_summary$hedia08 == 7 | wave_2_summary$hedia09 == 7, 1, 0)
wave_2_diabetes_simp <- subset(wave_2_summary, select = c(idauniq, diabetes2))

#wave 3
wave_3_summary$diabetes3 <- 0
wave_3_summary$diabetes3 <- ifelse(wave_3_summary$hediadi == 1, 1, 0)
wave_3_diabetes_simp <- subset(wave_3_summary, select = c(idauniq, diabetes3))


#wave_4
wave_4_summary$diabetes4 <- 0
wave_4_summary$diabetes4 <- ifelse(wave_4_summary$hediadi == 1, 1, 0)
wave_4_diabetes_simp <- subset(wave_4_summary, select = c(idauniq, diabetes4))

#wave_5
wave_5_summary$diabetes5 <- 0
wave_5_summary$diabetes5 <- ifelse(wave_5_summary$hediadi == 1, 1, 0)
wave_5_diabetes_simp <- subset(wave_5_summary, select = c(idauniq, diabetes5))

#wave_6
wave_6_summary$diabetes6 <- 0
wave_6_summary$diabetes6 <- ifelse(wave_6_summary$hediadi == 1, 1, 0)
wave_6_diabetes_simp <- subset(wave_6_summary, select = c(idauniq, diabetes6))

#WAVE7
wave_7_summary$diabetes7 <- 0
wave_7_summary$diabetes7 <- ifelse(wave_7_summary$hediadi == 1, 1, 0)
wave_7_diabetes_simp <- subset(wave_7_summary, select = c(idauniq, diabetes7))


#wave_8
wave_8_summary$diabetes8 <- 0
wave_8_summary$diabetes8 <- ifelse(wave_8_summary$hediadi == 1, 1, 0)
wave_8_diabetes_simp <- subset(wave_8_summary, select = c(idauniq, diabetes8))

#wave_9
wave_9_summary$diabetes9 <- 0
wave_9_summary$diabetes9 <- ifelse(wave_9_summary$hediadi == 1, 1, 0)
wave_9_diabetes_simp <- subset(wave_9_summary, select = c(idauniq, diabetes9))

#combine diabetes
diabetes_longitudinal <- merge(wave_1_diabetes_simp, wave_2_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_3_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_4_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_5_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_6_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_7_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_8_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
diabetes_longitudinal <- merge(diabetes_longitudinal, wave_9_diabetes_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)



#statin
#wave 1 - no data
wave_1_summary$statin1 <- NA
wave_1_statin_simp <- subset(wave_1_summary, select = c(idauniq, statin1))


#wave 2 - no data
wave_2_summary$statin2 <- NA
wave_2_statin_simp <- subset(wave_2_summary, select = c(idauniq, statin2))

#wave 3
wave_3_summary$statin3 <- 0
wave_3_summary$statin3 <- ifelse(wave_3_summary$hechmd == 1, 1, NA)
wave_3_summary$statin3 <- ifelse(wave_3_summary$hechmd == 2 | wave_3_summary$hechmd == -1, 0, wave_3_summary$statin3)
wave_3_statin_simp <- subset(wave_3_summary, select = c(idauniq, statin3))


#wave_4
wave_4_summary$statin4 <- 0
wave_4_summary$statin4 <- ifelse(wave_4_summary$hechmd == 1, 1, NA)
wave_4_summary$statin4 <- ifelse(wave_4_summary$hechmd == 2 | wave_4_summary$hechmd == -1, 0, wave_4_summary$statin4)
wave_4_statin_simp <- subset(wave_4_summary, select = c(idauniq, statin4))

#wave_5
wave_5_summary$statin5 <- 0
wave_5_summary$statin5 <- ifelse(wave_5_summary$hechmd == 1, 1, NA)
wave_5_summary$statin5 <- ifelse(wave_5_summary$hechmd == 2 | wave_5_summary$hechmd == -1, 0, wave_5_summary$statin5)
wave_5_statin_simp <- subset(wave_5_summary, select = c(idauniq, statin5))

#wave_6
wave_6_summary$statin6 <- 0
wave_6_summary$statin6 <- ifelse(wave_6_summary$HeChMd == 1, 1, NA)
wave_6_summary$statin6 <- ifelse(wave_6_summary$HeChMd == 2 | wave_6_summary$HeChMd == -1, 0, wave_6_summary$statin6)
wave_6_statin_simp <- subset(wave_6_summary, select = c(idauniq, statin6))

wave_6_statin_nurse <- subset(wave_6_meds, select = c(idauniq, lld))
wave_6_statin_comparison <- merge(wave_6_statin_nurse, wave_6_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)


#wave_7
wave_7_summary$statin7 <- 0
wave_7_summary$statin7 <- ifelse(wave_7_summary$HeChMd == 1, 1, NA)
wave_7_summary$statin7 <- ifelse(wave_7_summary$HeChMd == 2 | wave_7_summary$HeChMd == -1, 0, wave_7_summary$statin7)
wave_7_statin_simp <- subset(wave_7_summary, select = c(idauniq, statin7))

#wave_8 - get from medication prescription data
wave_8_summary$statin8 <- 0
wave_8_summary$statin8 <- ifelse(wave_8_summary$hechmd == 1, 1, NA)
wave_8_summary$statin8 <- ifelse(wave_8_summary$hechmd == 2 | wave_8_summary$hechmd == -1, 0, wave_8_summary$statin8)
wave_8_statin_simp <- subset(wave_8_summary, select = c(idauniq, statin8))

#wave_9 - get from medication prescription data
wave_9_summary$statin9 <- 0
wave_9_summary$statin9 <- ifelse(wave_9_summary$hechmd == 1, 1, NA)
wave_9_summary$statin9 <- ifelse(wave_9_summary$hechmd == 2 | wave_9_summary$hechmd == -1, 0, wave_9_summary$statin9)
wave_9_statin_simp <- subset(wave_9_summary, select = c(idauniq, statin9))

#combine statin
statin_longitudinal <- merge(wave_1_statin_simp, wave_2_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_3_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_4_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_5_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_6_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_7_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_8_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
statin_longitudinal <- merge(statin_longitudinal, wave_9_statin_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)



#cholesterol
#wave 1
wave_1_chol <- as.data.frame(wave_1$idauniq)
colnames(wave_1_chol) <- "idauniq"
wave_1_chol$chol1 <- NA


#wave_2
wave_2_chol <- subset(wave_2_nurse, select = c(idauniq, chol))
wave_2_chol$chol <- ifelse(wave_2_chol$chol < 0, NA, wave_2_chol$chol)
colnames(wave_2_chol)[2] <- "chol2"

#wave 3
wave_3_chol <- as.data.frame(wave_3$idauniq)
colnames(wave_3_chol) <- "idauniq"
wave_3_chol$chol3 <- NA

#wave 4
wave_4_chol <- subset(wave_4_nurse, select = c(idauniq, chol))
wave_4_chol$chol <- ifelse(wave_4_chol$chol < 0, NA, wave_4_chol$chol)
colnames(wave_4_chol)[2] <- "chol4"

#wave 5
wave_5_chol <- as.data.frame(wave_5$idauniq)
colnames(wave_5_chol) <- "idauniq"
wave_5_chol$chol5 <- NA

#wave_6
wave_6_chol <- subset(wave_6_nurse, select = c(idauniq, chol))
wave_6_chol$chol <- ifelse(wave_6_chol$chol < 0, NA, wave_6_chol$chol)
colnames(wave_6_chol)[2] <- "chol6"

#wave 7
wave_7_chol <- as.data.frame(wave_7$idauniq)
colnames(wave_7_chol) <- "idauniq"
wave_7_chol$chol7 <- NA

#wave_8
wave_8_chol <- subset(wave_8_nurse, select = c(idauniq, chol))
wave_8_chol$chol <- ifelse(wave_8_chol$chol < 0, NA, wave_8_chol$chol)
colnames(wave_8_chol)[2] <- "chol8"

#wave_9
columns <- c("idauniq", "chol")
wave_9_chol <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(wave_9_chol) <- columns

for (i in (1:nrow(wave_8_9_nurse))){
  if(!(wave_8_9_nurse$idauniq[i] %in% wave_8_nurse$idauniq)){
    a <- wave_8_9_nurse[i,]
    b <- subset(a, select = c(idauniq, chol))
    wave_9_chol <- rbind(wave_9_chol, b)
    rm(a)
    rm(b)
  }}

wave_9_chol$chol <- ifelse(wave_9_chol$chol < 0, NA, wave_9_chol$chol)
colnames(wave_9_chol)[2] <- "chol9"


#combine cholesterol
chol_longitudinal <- merge(wave_1_chol, wave_2_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_3_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_4_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_5_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_6_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_7_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_8_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)
chol_longitudinal <- merge(chol_longitudinal, wave_9_chol, by = "idauniq", all.x = TRUE, all.y = TRUE)

#drop individuals with no cholesterol measurements across all waves
chol_longitudinal$missing <- rowSums(is.na(chol_longitudinal))

chol_longitudinal_nomissing <- chol_longitudinal %>% dplyr::filter(missing < 8)

#haematological malignancy 
#wave 1 - no data
wave_1_summary$hm1 <- NA
wave_1_hm_simp <- subset(wave_1_summary, select = c(idauniq, hm1))


#wave 2
wave_2_summary$hm2 <- 0
wave_2_summary$hm2 <- ifelse(wave_2_summary$HeCana == 4 | wave_2_summary$HeCana == 5, 1, 0)
wave_2_hm_simp <- subset(wave_2_summary, select = c(idauniq, hm2))

#wave 3
wave_3_summary$hm3 <- 0
wave_3_summary$hm3 <- ifelse(wave_3_summary$heleuk == 1 | wave_3_summary$heleuk == 2 | wave_3_summary$hecanaa == 4 | wave_3_summary$hecanaa == 5, 1, 0)
wave_3_hm_simp <- subset(wave_3_summary, select = c(idauniq, hm3))


#wave_4
wave_4_summary$hm4 <- 0
wave_4_summary$hm4 <- ifelse(wave_4_summary$heleuk == 1 | wave_4_summary$heleuk == 2 | wave_4_summary$hecanaa == 4 | wave_4_summary$hecanaa == 5, 1, 0)
wave_4_hm_simp <- subset(wave_4_summary, select = c(idauniq, hm4))

#wave_5
wave_5_summary$hm5 <- 0
wave_5_summary$hm5 <- ifelse(wave_5_summary$heleuk == 1 | wave_5_summary$heleuk == 2 | wave_5_summary$hecanaa == 4 | wave_5_summary$hecanaa == 5, 1, 0)
wave_5_hm_simp <- subset(wave_5_summary, select = c(idauniq, hm5))

#wave_6
wave_6_summary$hm6 <- 0
wave_6_summary$hm6 <- ifelse(wave_6_summary$HeLeuk == 1 | wave_6_summary$HeLeuk == 2 | wave_6_summary$HeCanaa == 4 | wave_6_summary$HeCanaa == 5, 1, 0)
wave_6_hm_simp <- subset(wave_6_summary, select = c(idauniq, hm6))

#wave_7
wave_7_summary$hm7 <- 0
wave_7_summary$hm7 <- ifelse(wave_7_summary$HeLeuk == 1 | wave_7_summary$HeLeuk == 2 | wave_7_summary$HeCanaa == 4 | wave_7_summary$HeCanaa == 5, 1, 0)
wave_7_hm_simp <- subset(wave_7_summary, select = c(idauniq, hm7))

#wave_8
wave_8_summary$hm8 <- 0
wave_8_summary$hm8 <- ifelse(wave_8_summary$hecanbb == 1, 1, 0)
wave_8_hm_simp <- subset(wave_8_summary, select = c(idauniq, hm8))

#wave_9
wave_9_summary$hm9 <- 0
wave_9_summary$hm9 <- ifelse(wave_9_summary$hecanbb == 1, 1, 0)
wave_9_hm_simp <- subset(wave_9_summary, select = c(idauniq, hm9))

#combine hm
hm_longitudinal <- merge(wave_1_hm_simp, wave_2_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_3_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_4_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_5_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_6_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_7_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_8_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
hm_longitudinal <- merge(hm_longitudinal, wave_9_hm_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)


#age at each wave
wave_1_age <- subset(wave_1_summary, select = c(idauniq, indager))
colnames(wave_1_age)[2] <- "indager1"

wave_2_age <- subset(wave_2_summary, select = c(idauniq, indager))
colnames(wave_2_age)[2] <- "indager2"

wave_3_age <- subset(wave_3_summary, select = c(idauniq, indager))
colnames(wave_3_age)[2] <- "indager3"
wave_4_age <- subset(wave_4_summary, select = c(idauniq, indager))
colnames(wave_4_age)[2] <- "indager4"
wave_5_age <- subset(wave_5_summary, select = c(idauniq, indager))
colnames(wave_5_age)[2] <- "indager5"
wave_6_age <- subset(wave_6_summary, select = c(idauniq, indager))
colnames(wave_6_age)[2] <- "indager6"
wave_7_age <- subset(wave_7_summary, select = c(idauniq, indager))
colnames(wave_7_age)[2] <- "indager7"
wave_8_age <- subset(wave_8_summary, select = c(idauniq, indager))
colnames(wave_8_age)[2] <- "indager8"
wave_9_age <- subset(wave_9_summary, select = c(idauniq, indager))
colnames(wave_9_age)[2] <- "indager9"

#combine age
age_longitudinal <- merge(wave_1_age, wave_2_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_3_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_4_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_5_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_6_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_7_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_8_age, by = "idauniq", all.x = TRUE, all.y = TRUE)
age_longitudinal <- merge(age_longitudinal, wave_9_age, by = "idauniq", all.x = TRUE, all.y = TRUE)


#smoking status - current smoker
#wave 1
wave_1_summary$smoke1 <- 0
wave_1_summary$smoke1 <- ifelse(wave_1_summary$heska == 1, 1, 0)
wave_1_smoke_simp <- subset(wave_1_summary, select = c(idauniq, smoke1))



#wave 2
wave_2_summary$smoke2 <- 0
wave_2_summary$smoke2 <- ifelse(wave_2_summary$HESka == 1, 1, 0)
wave_2_smoke_simp <- subset(wave_2_summary, select = c(idauniq, smoke2))

#wave 3
wave_3_summary$smoke3 <- 0
wave_3_summary$smoke3 <- ifelse(wave_3_summary$heska == 1, 1, 0)
wave_3_smoke_simp <- subset(wave_3_summary, select = c(idauniq, smoke3))

#wave 4
wave_4_summary$smoke4 <- 0
wave_4_summary$smoke4 <- ifelse(wave_4_summary$heska == 1, 1, 0)
wave_4_smoke_simp <- subset(wave_4_summary, select = c(idauniq, smoke4))

#wave 5
wave_5_summary$smoke5 <- 0
wave_5_summary$smoke5 <- ifelse(wave_5_summary$heska == 1, 1, 0)
wave_5_smoke_simp <- subset(wave_5_summary, select = c(idauniq, smoke5))

#wave 6
wave_6_summary$smoke6 <- 0
wave_6_summary$smoke6 <- ifelse(wave_6_summary$HESka == 1, 1, 0)
wave_6_smoke_simp <- subset(wave_6_summary, select = c(idauniq, smoke6))

#wave 7
wave_7_summary$smoke7 <- 0
wave_7_summary$smoke7 <- ifelse(wave_7_summary$HESka == 1, 1, 0)
wave_7_smoke_simp <- subset(wave_7_summary, select = c(idauniq, smoke7))

#wave 8
wave_8_summary$smoke8 <- 0
wave_8_summary$smoke8 <- ifelse(wave_8_summary$heska == 1, 1, 0)
wave_8_smoke_simp <- subset(wave_8_summary, select = c(idauniq, smoke8))

#wave 9
wave_9_summary$smoke9 <- 0
wave_9_summary$smoke9 <- ifelse(wave_9_summary$heska == 1, 1, 0)
wave_9_smoke_simp <- subset(wave_9_summary, select = c(idauniq, smoke9))

#combine smoke
smoke_longitudinal <- merge(wave_1_smoke_simp, wave_2_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_3_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_4_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_5_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_6_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_7_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_8_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)
smoke_longitudinal <- merge(smoke_longitudinal, wave_9_smoke_simp, by = "idauniq", all.x = TRUE, all.y = TRUE)


#chip - just try binary present or absent at wave 2 
#wave 2 chip cases
wave_2_chip_cases <- cases_noQ1274L %>% dplyr::filter(elsa_wave == 2)
wave_2_chip_cases <- wave_2_chip_cases %>% dplyr::filter(!(is.na(idauniq)))

wave_2_controls <- controls %>% dplyr::filter(elsa_wave == 2)
wave_2_controls <- wave_2_controls %>% dplyr::filter(!(is.na(idauniq)))

#check for overlap
table(wave_2_chip_cases$idauniq %in% wave_2_controls$idauniq)
#FALSE 
#1064 (10 removed with excluding Q1274L)


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

#restrict to just wave 2 cases and controls
required_ids <- c(wave_2_chip_cases$idauniq, wave_2_controls$idauniq)

longitudinal_data <- longitudinal_data %>% dplyr::filter(idauniq %in% required_ids)

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


#cvd_wave 
ordered_long_data$cvd_wave <- NA
for (i in (1:nrow(ordered_long_data))) {
  ordered_long_data$cvd_wave[i] <- ifelse(ordered_long_data$ihd[i] == 1 | ordered_long_data$stroke[i] == 1, ordered_long_data$wave[i], NA)
}




#multiple imputation of cholesterol based on idauniq, sex, age and statin

chol <- subset(ordered_long_data, select = c(idauniq, wave, indager, chol, statin))
pred <- make.predictorMatrix(chol)


chol_imps <- mice(chol, predictorMatrix = pred)

chol_imputed <- complete(chol_imps)

#drop statin and rename vars
colnames(chol_imputed)[4] <- "chol_imp"
chol_imputed <- subset(chol_imputed, select = c(idauniq, wave, chol_imp))

#merge
ordered_long_data <- merge(ordered_long_data, chol_imputed, by = c("idauniq", "wave"))


#follow up wave 2 individuals - add variable for whether chip or control at wave 2
ordered_long_data$chip <- ifelse(ordered_long_data$idauniq %in% wave_2_chip_cases$idauniq, 1, 0)
ordered_long_data$chip <- as.factor(ordered_long_data$chip)
ordered_long_data$indsex <- as.factor(ordered_long_data$indsex)
ordered_long_data$statin <- as.factor(ordered_long_data$statin)






#prior statin - statin before cvd event
ordered_long_data$prior_statin <- 0
for (i in (1:nrow(ordered_long_data))) {
  ordered_long_data$prior_statin[i] <- ifelse(ordered_long_data$cvd[i] == 0 & ordered_long_data$statin[i] == 1, 1, 0)
}

#reorganise data with time to event
#just consider new diagnoses from wave 2 - not those pre existing first chip measurement 
wave1_cvd <- ordered_long_data %>% dplyr::filter(cvd_wave == 1)
ordered_long_data_w2on <- ordered_long_data %>% dplyr::filter(!(idauniq %in% wave1_cvd$idauniq))

#for each individual, specify idauniq, indsex, cvd, cvd_wave, wave 2 age, wave 2 cholesterol, htn_cont, diabetes_cont, smoke_cont, chip, prior statin cont
ids <- unique(ordered_long_data_w2on$idauniq)


#need to specify cvd wave and also last wave recorded
last_wave <- subset(ordered_long_data_w2on, select = c(idauniq, wave, indager))

last_wave_var <- data.frame(matrix(ncol = 2, nrow = 0))


for (i in ids) {
  a <- ordered_long_data_w2on %>% dplyr::filter(idauniq == i)
  row <- which.max(a$indager)
  last <- a$wave[row]
  vec <- c(i, last)
  last_wave_var <- rbind(last_wave_var, vec)
  rm(a)
  rm(row)
  rm(last)
  rm(vec)
}

colnames(last_wave_var) <- c("idauniq", "last_wave_recorded")

ordered_long_data_w2on <- merge(ordered_long_data_w2on, last_wave_var, by = "idauniq", all.x = TRUE)


ordered_long_data_w2on$time <- ifelse(is.na(ordered_long_data_w2on$cvd_wave), ordered_long_data_w2on$last_wave_recorded, ordered_long_data_w2on$cvd_wave)


df <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(df) <- c("idauniq", "indsex", "cvd", "time", "indager", "chol", "htn_cont", "diabetes_cont", "smoke_cont", "chip", "prior_statin")


for (i in ids) {
  a <- ordered_long_data_w2on %>% dplyr::filter(idauniq == i)
  vars <- c()
  vars[1] <- i
  vars[2] <- a$indsex[1]
  vars[3] <- max(a$cvd, na.rm = TRUE) 
  vars[4] <- min(a$time, na.rm = TRUE)
  vars[5] <- a$indager[1]
  vars[6] <- a$chol[2]
  vars[7] <- a$htn_cont[1]
  vars[8] <- a$diabetes_cont[1]
  vars[9] <- a$smoke_cont[1]
  vars[10] <- a$chip[1]
  vars[11] <- ifelse(mean(a$prior_statin, na.rm = TRUE) > 0, 1, 0)
  vars.df <- as.data.frame(t(vars))
  colnames(vars.df) <- c("idauniq", "indsex", "cvd", "cvd_wave", "indager", "chol", "htn_cont", "diabetes_cont", "smoke_cont", "chip", "prior_statin")
  df <- rbind(df, vars)
  rm(a)
  rm(vars)
  rm(vars.df)
}

colnames(df) <- c("idauniq", "indsex", "cvd", "time", "indager", "chol", "htn_cont", "diabetes_cont", "smoke_cont", "chip", "prior_statin")


#remove nas in prior statin
df_comp <- df %>% dplyr::filter(!is.na(prior_statin))

df_comp$chip_fac <- as.factor(ifelse(df_comp$chip == 0, 1, 2))
df_comp$prior_statin_fac <- as.factor(ifelse(df_comp$prior_statin == 0, 1, 2))
df_comp$cvd2 <- as.numeric(as.character(df_comp$cvd))


#separate out those who received prior statins as primary prevention, with and without chip
prior_stat <- df_comp %>% dplyr::filter(prior_statin == 1)

prior_stat_fac <- prior_stat
prior_stat_fac$indsex <- as.factor(prior_stat_fac$indsex)
prior_stat_fac$htn_cont <- as.factor(prior_stat_fac$htn_cont)
prior_stat_fac$diabetes_cont <- as.factor(prior_stat_fac$diabetes_cont)
prior_stat_fac$smoke_cont <- as.factor(prior_stat_fac$smoke_cont)
prior_stat_fac$chip_fac <- as.factor(prior_stat_fac$chip)
prior_stat_fac$prior_statin_fac <- as.factor(prior_stat_fac$prior_statin)

#remove if chol NA
prior_stat_fac <- prior_stat_fac %>% dplyr::filter(!is.na(chol))


#whole case cohort

#check for baseline differences - age, sex, cholesterol, smoking, diabetes, htn
prior_stat_fac_cases <- prior_stat_fac %>% dplyr::filter(chip == 2)
prior_stat_fac_ctls <- prior_stat_fac %>% dplyr::filter(chip == 1)

#age 
t.test(prior_stat_fac_cases$indager, prior_stat_fac_ctls$indager)
sd(prior_stat_fac_ctls$indager)
sd(prior_stat_fac_cases$indager)


#sex 
sex.m <- as.matrix(table(prior_stat_fac$indsex, prior_stat_fac$chip))
chisq.test(sex.m)


#cholesterol 
t.test(prior_stat_fac_cases$chol, prior_stat_fac_ctls$chol)
sd(prior_stat_fac_cases$chol)
sd(prior_stat_fac_ctls$chol)


#smoking 
smoke.m <- as.matrix(table(prior_stat_fac$smoke_cont, prior_stat_fac$chip))
chisq.test(smoke.m)

#htn 
htn.m <- as.matrix(table(prior_stat_fac$htn_cont, prior_stat_fac$chip))
chisq.test(htn.m)


#diabetes 
dm.m <- as.matrix(table(prior_stat_fac$diabetes_cont, prior_stat_fac$chip))
chisq.test(dm.m)

cox_model <-  coxph(Surv(time, cvd) ~ chip_status + indager, data = prior_stat_fac, x = TRUE)

ph_check_whole <- cox.zph(cox_model)

treatment_mod <- glm(chip_fac ~ indager + indsex + smoke_cont, data = prior_stat_fac, family = binomial())
set.seed(145)
adjsurv <- adjustedsurv(data = prior_stat_fac, variable = "chip_fac",
                        ev_time = "time",
                        event = "cvd",
                        method = "iptw_km",
                        outcome_model = cox_model,
                        treatment_model = treatment_mod,
                        conf_int = TRUE,
                        times = c(3,4,5,6,7,8,9),
                        bootstrap = TRUE,
                        n_boot = 1000,
                        model = cox_model)

adjusted_curve_test(adjsurv, from = 3, to = 9)

plot(adjsurv, conf_int = TRUE, 
     ylab = "Adjusted ASCVD-free Probability",
     x_breaks = c(2,3,4,5,6,7,8,9),
     legend.title = "CHIP status",
     risk_table = TRUE,
     risk_table_stratify = TRUE,
     xlab = "ELSA wave")



