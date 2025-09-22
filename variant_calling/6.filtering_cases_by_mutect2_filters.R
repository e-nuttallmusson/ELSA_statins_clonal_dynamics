#script to filter variants by mutect2 filters, integrating information from longitudinal sampling

library(dplyr)
library(tidyr)
library(readxl)

#read in mutation data with variants carrying flags

setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions")
vars_inc_flags <- read.table("anon_chip_variant_data_clean_asxl1_fixed.txt", header = TRUE, sep = "\t")

#remove those with TET2 Q1274L
vars_inc_flags <- vars_inc_flags %>% dplyr::filter(NonsynOI != "Q1274L")

#separate out those with flags and those that pass (those that pass including asxl1 multiallelic)

#vars passing
vars_pass <- vars_inc_flags %>% dplyr::filter(Otherinfo10 == "PASS")

#asxl1 multiallelic - consider also pass as already filtered based on vaf
asxl1_multi_pass <- vars_inc_flags %>% dplyr::filter(Gene.refGene == "ASXL1")
asxl1_multi_pass <- asxl1_multi_pass %>% dplyr::filter(grepl("multiallelic", Otherinfo10))
#keep only those with vaf > 10%
asxl1_multi_pass <- asxl1_multi_pass %>% dplyr::filter(AF >= 0.1)

#vars flagged
vars_flagged <- vars_inc_flags %>% dplyr::filter(Otherinfo10 != "PASS")

#combine asxl1 with vars passing
vars_pass <- rbind(vars_pass, asxl1_multi_pass)

#vars flagged and not included
vars_not_passing <- anti_join(vars_flagged, vars_pass)

#read in variants that have been integrated overtime to show persistent drivers, independent of mutect2 flags
setwd("N:/My Documents/ELSA_data/leukaemia_letter_statins/revisions_2/exploring_flagged_vars")
persistent_drivers <- read.table("all_repeat_variant_flags_no_controls.txt", header = TRUE, sep = " ")

#identify those samples with flag that should be included for persistent presence - by sample, gene, NonsynOI and Otherinfo10
persistent_drivers_mixed <- persistent_drivers %>% dplyr::filter(status == "mixed")

#add which wave is nonpassing
persistent_drivers_mixed <- persistent_drivers_mixed %>%
  rowwise() %>%
  mutate(non_pass_timepoint = {
    vals <- c_across(all_of(oi_cols))
    vals[vals %in% c("NA","na","N/A","")] <- NA_character_
    tps <- timepts[!is.na(vals) & vals != "PASS"]
   if (length(tps) == 0) NA_integer_ else min(tps)  # earliest wave as integer
    }) %>%
   ungroup()

timepts <- c(2, 4, 6, 8, 9)
# 1) Add a row id and derive the wave per row 
with_wave <- persistent_drivers_mixed %>%
  mutate(.rid = row_number(),
         wave_chosen = map_int(
           str_extract_all(as.character(non_pass_timepoint), "\\d+"),
           ~ if (length(.x)) min(as.integer(.x)) else NA_integer_
         ))

# 2) Pivot AF_/Sample_/Otherinfo10_ to long (one row per row_id per timepoint)
long <- persistent_drivers_mixed %>%
  mutate(.rid = row_number()) %>%
  pivot_longer(
    cols = matches(paste0("^(AF|Sample|Otherinfo10)_(", paste(timepts, collapse="|"), ")$")),
    names_to = c(".value","wave"),
    names_pattern = "^(AF|Sample|Otherinfo10)_(\\d+)$"
  ) %>%
  mutate(wave = as.integer(wave))

# 3) Keep only the selected timepoint per original row
mixed_vars_to_keep<- long %>%
  inner_join(with_wave %>% select(.rid, wave_chosen),
             by = ".rid") %>%
  filter(wave == wave_chosen) %>%
  select(idauniq, Gene.refGene, NonsynOI, wave, AF, Sample, Otherinfo10)

keys <- c("Sample", "Gene.refGene", "NonsynOI", "AF")

mixed_to_keep_annovar <- vars_flagged %>%
  semi_join(
    mixed_vars_to_keep %>% select(all_of(keys)) %>% distinct() ,
    by = keys
  )


#persistent drivers will flags for all samples
persistent_drivers_allnonpass <- persistent_drivers %>% dplyr::filter(status == "all_nonPASS")

persistent_drivers_allnonpass_long <- persistent_drivers_allnonpass %>%
  pivot_longer(
    cols = matches("^(AF|Sample|Otherinfo10)_\\d+$"),
    names_to = c(".value","wave"),
    names_pattern = "^(AF|Sample|Otherinfo10)_(\\d+)$",
    values_drop_na = FALSE
  ) %>%
  mutate(Otherinfo10 = na_if(Otherinfo10, "NA")) %>%
  filter(!is.na(Sample), Sample != "") %>%
  select(Sample, Gene.refGene, NonsynOI, AF, Otherinfo10)  


keys <- c("Sample", "Gene.refGene", "NonsynOI", "AF")

nonpass_to_keep_annovar <- vars_flagged %>%
  semi_join(
    persistent_drivers_allnonpass_long %>% select(all_of(keys)) %>% distinct() ,
    by = keys
  )


 #read in manually reviewed
dnmt3a_nonpassing_to_keep <- read_xlsx("dnmt3a_variants_w_ctl_sample_for_igv_review.xlsx", sheet = "dnmt3a_vars_for_inclusion")
tet2_nonpassing_to_keep <- read_xlsx("tet2_variants_w_ctl_sample_for_igv_review.xlsx", sheet = "tet2_vars_for_inclusion")

#combine
manually_reviewed_to_keep <- rbind(dnmt3a_nonpassing_to_keep, tet2_nonpassing_to_keep)

#get full annovar listings
keys <- c("Sample", "Gene.refGene", "Chr", "Start", "End", "Ref", "Alt")

manual_nonpass_to_keep_annovar <- vars_flagged %>%
  semi_join(
    manually_reviewed_to_keep %>% select(all_of(keys)) %>% distinct() ,
    by = keys
  )


#now combine all vars to keep

vars_include_flags_reviewed <- rbind(vars_pass, mixed_to_keep_annovar, nonpass_to_keep_annovar, manual_nonpass_to_keep_annovar)


#output these variants
write.table(vars_include_flags_reviewed, "included_variants_after_flags_reviewed.txt", quote = FALSE, row.names = FALSE, sep = "\t")


