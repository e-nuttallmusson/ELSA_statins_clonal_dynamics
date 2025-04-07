#script to look at high vaf vs low vaf logistic regression for TET2 and DNMT3A


library(readxl)
library(dplyr)
library(glmtoolbox)

#import wave_8_9 chip mutation data, wave_8_9_comorbidities, drugs and blood parameters
setwd("N:/My Documents/ELSA_data/upgrade/")
chip_cases_w8_9_metadata <- read_xlsx("wave_8_and_9_chip_cases_comorbs_drugs_bloods.xlsx")

#add recent_vaf variable
chip_cases_w8_9_metadata$recent_vaf <- chip_cases_w8_9_metadata$AF_9
chip_cases_w8_9_metadata$recent_vaf <- ifelse(is.na(chip_cases_w8_9_metadata$recent_vaf), chip_cases_w8_9_metadata$AF_8, chip_cases_w8_9_metadata$recent_vaf)

#add mutation count variable
chip_cases_w8_9_metadata <- chip_cases_w8_9_metadata %>% group_by(idauniq) %>% mutate(count = n())



************************************TET2

#just consider TET2 chip - consider largest vaf if more than one tet2 mutation present
chip_tet2_mut <- chip_cases_w8_9_metadata %>% dplyr::filter(grepl("TET2", Gene.refGene))

chip_tet2_mut_dups <- chip_tet2_mut %>% group_by(idauniq) %>% filter(n() > 1) %>% arrange(idauniq, desc(recent_vaf))

chip_tet2_mut_dups2 <- chip_tet2_mut_dups[!duplicated(chip_tet2_mut_dups$idauniq),] 

#get non duplicated data
chip_tet2_mut_no_dups <- chip_tet2_mut %>% group_by(idauniq) %>% filter(n() == 1) %>% ungroup()

#recombine

chip_tet2_mut_no_dups <- rbind(chip_tet2_mut_no_dups, chip_tet2_mut_dups2)

#add variable for variant type - nonsense, fs etc
chip_tet2_mut_no_dups$variant_type <- NA
chip_tet2_mut_no_dups$variant_type[grep("fs*", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[grep("X", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[grep("nan", chip_tet2_mut_no_dups$NonsynOI)] <- "splice_variant"
chip_tet2_mut_no_dups$variant_type[grep("del", chip_tet2_mut_no_dups$NonsynOI)] <- "indel"
chip_tet2_mut_no_dups$variant_type[grep("\\*$", chip_tet2_mut_no_dups$NonsynOI)] <- "nonsense"
chip_tet2_mut_no_dups$variant_type[is.na(chip_tet2_mut_no_dups$variant_type)] <- "snv"

#replace censored ages
chip_tet2_mut_no_dups$indager <- replace(chip_tet2_mut_no_dups$indager, which(chip_tet2_mut_no_dups$indager == -7), 90)
chip_tet2_mut_no_dups$mapval <- replace(chip_tet2_mut_no_dups$mapval, which(chip_tet2_mut_no_dups$mapval < 0), NA)

#convert boolean to 1/0
cols <- sapply(chip_tet2_mut_no_dups, is.logical)
chip_tet2_mut_no_dups[,cols]<- lapply(chip_tet2_mut_no_dups[,cols], as.numeric)

#drop wave AF vals as working from recent_vaf
chip_tet2_mut_no_dups <- subset(chip_tet2_mut_no_dups, select = -c(AF_2, AF_4, AF_6, AF_8, AF_9))


#5 individuals without cholesterol value - replace with median
imps <- lapply(chip_tet2_mut_no_dups, median, na.rm = TRUE)

for (i in colnames(chip_tet2_mut_no_dups)) {
  chip_tet2_mut_no_dups[,i][is.na(chip_tet2_mut_no_dups[,i])] <- as.numeric(imps[i])
}

#drop those with blood disorder - 7 individuals
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% dplyr::filter(blood_dis == 0)

#remove 3 samples who are early wave cases with no detectable variants in w8/9
summary(chip_tet2_mut_no_dups$recent_vaf)
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% dplyr::filter(recent_vaf != 0.01)

#make cvd variable combining ihd, heart_dis and stroke.

chip_tet2_mut_no_dups$cvd <- rep(0, nrow(chip_tet2_mut_no_dups))
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))

chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#combine acei and arb to RAAS
chip_tet2_mut_no_dups$RAAS <- rep(0, nrow(chip_tet2_mut_no_dups))
chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% mutate(RAAS = case_when(acei == 1 | arb == 1 ~ 1, is.na(RAAS) ~ 0))

chip_tet2_mut_no_dups <- chip_tet2_mut_no_dups %>% mutate(RAAS = ifelse(is.na(RAAS), 0, RAAS))

#demographics of cohort
summary(chip_tet2_mut_no_dups$indager)
sd(chip_tet2_mut_no_dups$indager)

table(chip_tet2_mut_no_dups$indsex)

table(chip_tet2_mut_no_dups$current_smoker)

table(chip_tet2_mut_no_dups$vaf_cat)

table(chip_tet2_mut_no_dups$ihd)

table(chip_tet2_mut_no_dups$heart_dis)

table(chip_tet2_mut_no_dups$stroke)

table(chip_tet2_mut_no_dups$htn)

table(chip_tet2_mut_no_dups$diabetes)

table(chip_tet2_mut_no_dups$cancer)

table(chip_tet2_mut_no_dups$statins)

table(chip_tet2_mut_no_dups$antiplatelets)

table(chip_tet2_mut_no_dups$bb)

table(chip_tet2_mut_no_dups$RAAS)

#import fbc results
setwd("N:/My Documents/ELSA_data/FBC_results")
fbc <- read.csv("fbc_clean_20250123.csv", header = TRUE)

#keep waves 8 and 9

fbc_8 <- fbc %>% dplyr::filter(wave == 8)
fbc_9 <- fbc %>% dplyr::filter(wave == 9)

table(fbc_8$ID %in% fbc_9$ID)

fbc_8_tet2 <- fbc_8 %>% dplyr::filter(ID %in% chip_tet2_mut_no_dups$idauniq)
fbc_9_tet2 <- fbc_9 %>% dplyr::filter(ID %in% chip_tet2_mut_no_dups$idauniq)

late_wave_fbc <- rbind(fbc_8_tet2, fbc_9_tet2)

colnames(late_wave_fbc)[1] <- "idauniq"

chip_tet2_mut_no_dups_fbc <- merge(chip_tet2_mut_no_dups, late_wave_fbc, by = "idauniq", all.x = TRUE)

chip_tet2_mut_no_dups_fbc$HGB_g.L <- as.numeric(chip_tet2_mut_no_dups_fbc$HGB_g.L)
summary(chip_tet2_mut_no_dups_fbc$HGB_g.L)
sd(chip_tet2_mut_no_dups_fbc$HGB_g.L, na.rm = TRUE)

chip_tet2_mut_no_dups_fbc$NEUA_10.9.L <- as.numeric(chip_tet2_mut_no_dups_fbc$NEUA_10.9.L)
summary(chip_tet2_mut_no_dups_fbc$NEUA_10.9.L)
sd(chip_tet2_mut_no_dups_fbc$NEUA_10.9.L, na.rm = TRUE)

chip_tet2_mut_no_dups_fbc$PLT_10.9.L <- as.numeric(chip_tet2_mut_no_dups_fbc$PLT_10.9.L)
summary(chip_tet2_mut_no_dups_fbc$PLT_10.9.L)
sd(chip_tet2_mut_no_dups_fbc$PLT_10.9.L, na.rm = TRUE)

summary(chip_tet2_mut_no_dups_fbc$chol)
sd(chip_tet2_mut_no_dups_fbc$chol)

##################logistic regression
#binary outcome - low vaf high vaf
chip_tet2_mut_no_dups$vaf_cat <- ifelse(chip_tet2_mut_no_dups$recent_vaf <0.1, 0, 1)
#make outcome variable a 2 level factor
chip_tet2_mut_no_dups$vaf_cat <- as.factor(chip_tet2_mut_no_dups$vaf_cat)


#############statin focussed logistic regression 
#include: statins, age, sex, current smoker, htn, antiplatelet, cvd, bb, raas, 

tet2_statin_model <- glm(vaf_cat ~ indager + indsex + cvd + statins + htn + chol + antiplatelets + bb + RAAS, family = binomial(link = "logit"), data = chip_tet2_mut_no_dups)
summary(tet2_statin_model)

exp(cbind(OR=coef(tet2_statin_model), confint(tet2_statin_model)))



#check assumptions
#check for multicolinearity

gvif(tet2_statin_model)

#test for linearity between continuous vars and logits 

logit_check <- subset(chip_tet2_mut_no_dups, select = c(indager, chol))
predictors <- colnames(logit_check)

#save predicted probabilities
logit_check$probabilities <- statin_model$fitted.values

library(tidyr)
#calculate logit values and tidy data for plot
logit_check <- logit_check %>% mutate(logit = log(probabilities/(1-probabilities))) %>%
  select(-probabilities)%>%
  gather(key = "predictors", value = "predictor.value", -logit)

#plot
ggplot(logit_check, aes(y = logit, x = predictor.value)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") +
  theme_bw() +
  facet_wrap(~predictors, scales = "free_x")


#age adjusted univariate analyses
#statins
statin_uni_model <- glm(vaf_cat ~ indager + statins, family = binomial(link = "logit"), data = chip_tet2_mut_no_dups)
summary(statin_uni_model)
exp(cbind(OR=coef(statin_uni_model), confint(statin_uni_model)))

#antiplatelets
antiplt_uni_model <- glm(vaf_cat ~ indager + antiplatelets, family = binomial(link = "logit"), data = chip_tet2_mut_no_dups)
summary(antiplt_uni_model)
exp(cbind(OR=coef(antiplt_uni_model), confint(antiplt_uni_model)))



**************************************DNMT3A

#just consider dnmt3a chip - consider largest vaf if more than one dnmt3a mutation present
chip_dnmt3a_mut <- chip_cases_w8_9_metadata %>% dplyr::filter(grepl("DNMT3A", Gene.refGene))

chip_dnmt3a_mut_dups <- chip_dnmt3a_mut %>% group_by(idauniq) %>% filter(n() > 1) %>% arrange(idauniq, desc(recent_vaf))

chip_dnmt3a_mut_dups2 <- chip_dnmt3a_mut_dups[!duplicated(chip_dnmt3a_mut_dups$idauniq),] 

#get non duplicated data
chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut %>% group_by(idauniq) %>% filter(n() == 1) %>% ungroup()

#recombine

chip_dnmt3a_mut_no_dups <- rbind(chip_dnmt3a_mut_no_dups, chip_dnmt3a_mut_dups2)

#add variable for variant type - nonsense, fs etc
chip_dnmt3a_mut_no_dups$variant_type <- NA
chip_dnmt3a_mut_no_dups$variant_type[grep("fs*", chip_dnmt3a_mut_no_dups$NonsynOI)] <- "nonsense"
chip_dnmt3a_mut_no_dups$variant_type[grep("X", chip_dnmt3a_mut_no_dups$NonsynOI)] <- "nonsense"
chip_dnmt3a_mut_no_dups$variant_type[grep("nan", chip_dnmt3a_mut_no_dups$NonsynOI)] <- "splice_variant"
chip_dnmt3a_mut_no_dups$variant_type[grep("del", chip_dnmt3a_mut_no_dups$NonsynOI)] <- "indel"
chip_dnmt3a_mut_no_dups$variant_type[grep("\\*$", chip_dnmt3a_mut_no_dups$NonsynOI)] <- "nonsense"
chip_dnmt3a_mut_no_dups$variant_type[is.na(chip_dnmt3a_mut_no_dups$variant_type)] <- "snv"

#replace censored ages
chip_dnmt3a_mut_no_dups$indager <- replace(chip_dnmt3a_mut_no_dups$indager, which(chip_dnmt3a_mut_no_dups$indager == -7), 90)
chip_dnmt3a_mut_no_dups$mapval <- replace(chip_dnmt3a_mut_no_dups$mapval, which(chip_dnmt3a_mut_no_dups$mapval < 0), NA)

#convert boolean to 1/0
cols <- sapply(chip_dnmt3a_mut_no_dups, is.logical)
chip_dnmt3a_mut_no_dups[,cols]<- lapply(chip_dnmt3a_mut_no_dups[,cols], as.numeric)

#drop wave AF vals as working from recent_vaf
chip_dnmt3a_mut_no_dups <- subset(chip_dnmt3a_mut_no_dups, select = -c(AF_2, AF_4, AF_6, AF_8, AF_9))

#replace nas with median
#4 without cholesterol measurements
imps <- lapply(chip_dnmt3a_mut_no_dups, median, na.rm = TRUE)

for (i in colnames(chip_dnmt3a_mut_no_dups)) {
  chip_dnmt3a_mut_no_dups[,i][is.na(chip_dnmt3a_mut_no_dups[,i])] <- as.numeric(imps[i])
}

#drop those with blood disorder - 3 individuals
chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% dplyr::filter(blood_dis == 0)

#remove 6 samples who are early wave cases with no detectable variants in w8/9
summary(chip_dnmt3a_mut_no_dups$recent_vaf)
chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% dplyr::filter(recent_vaf != 0.01)

#make cvd variable combining ihd, heart_dis and stroke.

chip_dnmt3a_mut_no_dups$cvd <- rep(0, nrow(chip_dnmt3a_mut_no_dups))
chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% mutate(cvd = case_when(ihd == 1 | stroke == 1 | heart_dis == 1 ~ 1, is.na(cvd) ~ 0))

chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% mutate(cvd = ifelse(is.na(cvd), 0, cvd))

#combine acei and arb to RAAS
chip_dnmt3a_mut_no_dups$RAAS <- rep(0, nrow(chip_dnmt3a_mut_no_dups))
chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% mutate(RAAS = case_when(acei == 1 | arb == 1 ~ 1, is.na(RAAS) ~ 0))

chip_dnmt3a_mut_no_dups <- chip_dnmt3a_mut_no_dups %>% mutate(RAAS = ifelse(is.na(RAAS), 0, RAAS))

#demographics
summary(chip_dnmt3a_mut_no_dups$indager)
sd(chip_dnmt3a_mut_no_dups$indager)

table(chip_dnmt3a_mut_no_dups$indsex)

table(chip_dnmt3a_mut_no_dups$current_smoker)

table(chip_dnmt3a_mut_no_dups$vaf_cat)

table(chip_dnmt3a_mut_no_dups$ihd)

table(chip_dnmt3a_mut_no_dups$stroke)

table(chip_dnmt3a_mut_no_dups$htn)

table(chip_dnmt3a_mut_no_dups$diabetes)

table(chip_dnmt3a_mut_no_dups$cancer)

table(chip_dnmt3a_mut_no_dups$statins)

table(chip_dnmt3a_mut_no_dups$antiplatelets)

table(chip_dnmt3a_mut_no_dups$bb)

table(chip_dnmt3a_mut_no_dups$RAAS)

#add in fbc data
fbc_8_dnmt3a <- fbc_8 %>% dplyr::filter(ID %in% chip_dnmt3a_mut_no_dups$idauniq)
fbc_9_dnmt3a <- fbc_9 %>% dplyr::filter(ID %in% chip_dnmt3a_mut_no_dups$idauniq)

late_wave_fbc <- rbind(fbc_8_dnmt3a, fbc_9_dnmt3a)

colnames(late_wave_fbc)[1] <- "idauniq"

chip_dnmt3a_mut_no_dups_fbc <- merge(chip_dnmt3a_mut_no_dups, late_wave_fbc, by = "idauniq", all.x = TRUE)

chip_dnmt3a_mut_no_dups_fbc$HGB_g.L <- as.numeric(chip_dnmt3a_mut_no_dups_fbc$HGB_g.L)
summary(chip_dnmt3a_mut_no_dups_fbc$HGB_g.L)
sd(chip_dnmt3a_mut_no_dups_fbc$HGB_g.L, na.rm = TRUE)

chip_dnmt3a_mut_no_dups_fbc$NEUA_10.9.L <- as.numeric(chip_dnmt3a_mut_no_dups_fbc$NEUA_10.9.L)
summary(chip_dnmt3a_mut_no_dups_fbc$NEUA_10.9.L)
sd(chip_dnmt3a_mut_no_dups_fbc$NEUA_10.9.L, na.rm = TRUE)

chip_dnmt3a_mut_no_dups_fbc$PLT_10.9.L <- as.numeric(chip_dnmt3a_mut_no_dups_fbc$PLT_10.9.L)
summary(chip_dnmt3a_mut_no_dups_fbc$PLT_10.9.L)
sd(chip_dnmt3a_mut_no_dups_fbc$PLT_10.9.L, na.rm = TRUE)

summary(chip_dnmt3a_mut_no_dups_fbc$chol)
sd(chip_dnmt3a_mut_no_dups_fbc$chol)

##################logistic regression
#binary outcome - low vaf high vaf
chip_dnmt3a_mut_no_dups$vaf_cat <- ifelse(chip_dnmt3a_mut_no_dups$recent_vaf <0.1, 0, 1)

table(chip_dnmt3a_mut_no_dups$cvd, chip_dnmt3a_mut_no_dups$vaf_cat)

#make outcome variable a 2 level factor
chip_dnmt3a_mut_no_dups$vaf_cat <- as.factor(chip_dnmt3a_mut_no_dups$vaf_cat)
chip_dnmt3a_mut_no_dups$indsex_cat <- as.factor(chip_dnmt3a_mut_no_dups$indsex)

#############statin focussed logistic regression with kmer cross validation
#include: statins, age, sex, current smoker, htn, antiplatelet, cvd, bb, raas, 

#run logistic regression

dnmt3a_statin_model <- glm(vaf_cat ~ indager + indsex_cat + cvd + statins + htn + chol + antiplatelets + bb + RAAS, family = binomial(link = "logit"), data = chip_dnmt3a_mut_no_dups)
summary(dnmt3a_statin_model)

exp(cbind(OR=coef(dnmt3a_statin_model), confint(dnmt3a_statin_model)))


#check assumptions
#check for multicolinearity
gvif(dnmt3a_statin_model)


#test for linearity between continuous vars and logits 

logit_check <- subset(chip_dnmt3a_mut_no_dups, select = c(indager, chol))
predictors <- colnames(logit_check)

#save predicted probabilities
logit_check$probabilities <- statin_model$fitted.values

library(tidyr)
#calculate logit values and tidy data for plot
logit_check <- logit_check %>% mutate(logit = log(probabilities/(1-probabilities))) %>%
  select(-probabilities)%>%
  gather(key = "predictors", value = "predictor.value", -logit)

#plot
ggplot(logit_check, aes(y = logit, x = predictor.value)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") +
  theme_bw() +
  facet_wrap(~predictors, scales = "free_x")


