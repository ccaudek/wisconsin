# Script name: 001_gen_input_hddm_stim.R
# Project: WCST
# Script purpose: create single file with all raw PRL data.
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed Jun  2 10:07:50 2021
# Last Modified Date: Sat Jun 18 10:04:13 2022
#
# Convergence problems:
# ch_pi_2001_10_08_418_f (AN)

# Prelims
suppressPackageStartupMessages({
  library("here")
  library("tidyverse")
  library("stringi")
  library("readxl")
})

# Increase max print
options(max.print = .Machine$integer.max)

# source(here("code", "functions", "funs_gen_data_for_hddm.R"))
source(here::here("src", "R", "functions", "funs_prl.R"))


# 1. Merge psytoolkit PRL files  ------------------------------------------

# Read individual psytoolkit PRL files and generate 4 RDS files.
# E.g. here("data", "processed", "prl", "data_social", 
# "controls_with_psychtoolkit_code.rds")

## load_psychtoolkit_files()


# 2.  Create data list ----------------------------------------------------

# Read the RDS files files with raw PRL data for each group and condition and 
# create a data_list. We also add the participant's identifier (es., 
# ca_po_2002_05_25_700_f).

# The data_list includes one data.frame for each condition:
#
# data_list[[1]] : df_food_patients
# data_list[[2]] : df_food_controls
# data_list[[3]] : df_social_patients
# data_list[[4]] : df_social_controls

## data_list <- write_prl_raw_data_list()


# 3. Corrected subj_name --------------------------------------------------

# For each element of the list, we correct the subject's identifier.

## d_list <- correct_subj_names(data_list)

# unique(d_list[[1]]$subj_name)
# unique(d_list[[2]]$subj_name)
# unique(d_list[[3]]$subj_name)
# unique(d_list[[4]]$subj_name)
# The returned list has the same structure as data_list. 


# 4. Clean and save single file --------------------------------------------

# Remove PRL sessions with too many NAs, bind the data frames in list and save 
# cleaned file in 
# here("data", "processed", "prl", "complete_cleaned_raw_data", 
# "raw_data_prl_both_groups.rds")

## binding_cleaned_data_frames(d_list)


# 5.  Add diagnostic category ---------------------------------------------

# here::here("data", "processed", "prl", "raw_prl_data", "prl_tot_raw_data.rds")
## add_diagnostic_category()

# All subjects (patients and controls)
prl_dat <- readRDS(
  here::here(
    "data", "processed", "prl", "raw_prl_data", "prl_tot_raw_data.rds")
  )

length(unique(prl_dat$subj_name))
# [1] 338


# Data cleaning ----------

# Recode feedback
prl_dat$feedback <- ifelse(prl_dat$feedback == 3, NA, prl_dat$feedback)
prl_dat$feedback <- ifelse(prl_dat$feedback == 2, 0, prl_dat$feedback)

prl_dat$feedback <- ifelse(
  prl_dat$rt < 200 | prl_dat$rt > 2499, NA, prl_dat$feedback
)

# Reaction times
prl_dat$rt <- ifelse(
  prl_dat$rt < 200 | prl_dat$rt > 2499, NA, prl_dat$rt
)
prl_dat$rt <- ifelse(prl_dat$trial == 1, NA, prl_dat$rt)

# For diag_cat, replace NA with HC
my_vector <- prl_dat$diag_cat
my_vector[is.na(my_vector)] <- "HC"
prl_dat$diag_cat <- my_vector
prl_dat$diag_cat <- factor(prl_dat$diag_cat)

# Multiple imputation
data_for_imp <- prl_dat %>% 
  dplyr::select(
    subj_name, stimulus_type, epoch, target_position, 
    which_image_is_rewarded_more,
    keypress, rt, feedback, is_target_img_rewarded_in_first_epoch, 
    is_target_rewared_in_present_epoch, position_target_img, resp, 
    is_target_img_chosen, trial, group, diag_cat
  )

set.seed(123)
miceObj <- miceRanger::miceRanger(
  data_for_imp
  , m = 1
  , returnModels = TRUE
  , verbose = TRUE
)

dataList <- miceRanger::completeData(miceObj)

d <- dataList[[1]]

# These are the cleaned data
prl_dat$feedback <- d$feedback
prl_dat$rt <- d$rt

prl_dat |> 
  group_by(diag_cat) |> 
  summarize(
    n = n_distinct(subj_name)
  )

# ----------

unique(prl_dat[prl_dat$group == "patients", ]$subj_name) |>
  sort()
# [1] "al_ca_1996_03_27_621_f" "al_ro_1989_04_25_160_f" "al_zu_1997_04_02_880_f"
# [4] "am_gu_1999_02_11_937_f" "an_am_1996_05_12_176_f" "an_de_1998_11_10_289_f"
# [7] "ar_ce_2005_04_20_937_f" "ar_co_1996_12_27_348_f" "as_ga_2005_06_15_329_f"
# [10] "au_ru_1998_09_21_806_f" "be_ma_1999_06_15_475_f" "bi_an_2001_09_16_735_f"
# [13] "bi_di_2006_04_20_725_f" "ca_fa_1996_03_26_092_f" "ca_po_2002_05_25_700_f"
# [16] "ca_so_2001_01_09_118_f" "ch_be_1990_12_20_153_f" "ch_br_1993_10_04_623_f"
# [19] "ch_ca_2000_09_26_406_f" "ch_ma_2001_10_27_332_f" "ch_na_2007_06_23_908_f"
# [22] "ch_pi_2001_10_08_418_f" "ch_pi_2004_02_25_126_f" "ch_ri_1993_05_05_564_f"
# [25] "cl_pu_2007_05_24_423_f" "cr_gi_1994_10_14_378_f" "cr_pa_1969_04_12_179_f"
# [28] "da_de_1998_08_15_141_m" "de_sc_1992_07_02_116_f" "el_ma_1986_06_14_839_f"
# [31] "em_al_1989_07_27_200_f" "em_bi_2007_12_28_766_f" "em_gr_2002_08_25_628_f"
# [34] "em_or_2003_01_02_101_f" "es_bo_2004_07_23_474_f" "fe_al_1988_05_06_180_f"
# [37] "fe_ma_1998_06_29_257_f" "fe_sa_2002_05_09_008_f" "fr_au_1987_12_16_221_f"
# [40] "fr_bo_1993_09_09_170_f" "fr_la_2004_05_17_363_f" "fr_ro_1982_08_15_048_f"
# [43] "ga_gi_2003_02_09_229_f" "gi_ba_2008_01_31_376_f" "gi_ca_2006_10_14_101_f"
# [46] "gi_ma_1999_09_26_585_f" "gi_to_1996_02_02_043_f" "gi_va_1992_04_14_174_f"
# [49] "gi_za_1992_09_07_575_f" "gr_bo_1996_07_31_547_f" "gr_de_2002_09_21_426_f"
# [52] "he_ha_2006_04_21_874_f" "il_fu_2002_12_30_306_f" "ir_bo_1981_03_29_325_f"
# [55] "ir_pi_2002_01_22_765_f" "ir_to_2007_08_01_838_f" "ir_ve_2004_02_09_500_f"
# [58] "lu_mu_1997_03_18_059_f" "lu_te_1990_10_28_496_f" "ma_ba_1995_05_25_321_f"
# [61] "ma_be_1997_09_01_726_f" "ma_va_1998_07_04_538_f" "ma_za_2002_02_28_051_f"
# [64] "ra_al_2002_10_05_370_f" "sa_ca_2004_06_11_885_f" "sa_ta_2003_11_14_150_f"
# [67] "so_be_2008_12_15_399_f"

# patients_prl_names <- unique(prl_dat[prl_dat$group == "patients", ]$subj_name)

# controls_prl_names <- 
#   unique(prl_dat[prl_dat$group == "controls", ]$subj_name) |> 
#   sort()


# Read the file with lists of participants (patients and controls) that 
# completed both the WCST and the PRL

source(here::here(
  "src", "R", "scripts", "scripts_prl_wcst", 
  "01_preprocessing", "002_select_subjects.R")
)

# Select PRL controls both experiments (prl and wcst)
controls_both_exps_df <- 
  prl_dat[prl_dat$subj_name %in% controls_both_exps_names, ]

controls_both_exps_df$group <- factor(controls_both_exps_df$group)


# Select control subjects who performed above chance level in the PRL task

bysubj_feedback <- controls_both_exps_df |> 
  group_by(subj_name) |> 
  summarize(
    p = mean(feedback)
  )

bysubj_feedback_good <- bysubj_feedback[bysubj_feedback$p > 0.5, ]

subj_above_threshold_prl_names <- bysubj_feedback_good$subj_name

controls_both_exps_df[
  controls_both_exps_df$subj_name %in% subj_above_threshold_prl_names, ] |> 
  group_by(diag_cat) |> 
  summarize(
    n = n_distinct(subj_name)
  )

# Subjects at risk
prl_ri <- controls_both_exps_df[
  (controls_both_exps_df$subj_name %in% subj_above_threshold_prl_names) & 
    controls_both_exps_df$diag_cat == "RI", 
  ]

at_risk_names <- unique(prl_ri$subj_name) |> sort()
at_risk_names
# [1] "al_be_1997_03_10_966_f" "an_re_2001_08_28_633_f" "ca_mi_2001_06_16_988_f"
# [4] "ch_lo_2000_09_25_565_f" "ch_ma_1995_08_28_639_f" "ed_sc_2001_09_26_034_m"
# [7] "el_li_1999_09_08_687_f" "gi_ga_2001_07_20_277_f" "gi_gi_1990_03_28_384_f"
# [10] "gi_me_2001_03_31_627_f" "gi_po_1998_11_07_576_f" "gi_se_1999_07_09_402_f"
# [13] "gi_sp_1995_10_16_533_f" "ha_ri_2001_07_07_704_f" "li_li_2001_12_04_406_f"
# [16] "ma_la_2001_09_12_609_f" "ma_ma_2001_07_10_611_f" "ma_pa_2001_06_11_636_f"
# [19] "ma_pi_2001_05_11_566_f" "ma_ta_2001_05_23_401_f" "ot_na_1999_03_18_271_f"
# [22] "sa_la_1994_11_13_963_f" "se_pi_2001_01_22_920_f" "vi_mi_2000_08_21_472_f"


# Select HC subjects that completed both tasks and with a PRL performance > 0.5
prl_hc <- controls_both_exps_df[
  (controls_both_exps_df$subj_name %in% subj_above_threshold_prl_names) & 
    controls_both_exps_df$diag_cat == "HC", 
]

prl_hc_names <- unique(prl_hc$subj_name)

# birth_date <- substr(prl_hc_names, 7, 10) |> as.numeric()
# 
# age <- 2022 - birth_date
# 
# age_df <- data.frame(
#   subj_name = prl_hc_names, age
# )
# 
# prl_hc <- left_join(prl_hc, age_df, by = "subj_name")

# add gender
sex <- stringr::str_sub(prl_hc_names, - 1, - 1) 

sex_df <- data.frame(
  subj_name = prl_hc_names, sex
)

prl_hc <- left_join(prl_hc, sex_df, by = "subj_name")

# Remove males
prl_hc <- prl_hc |> 
  dplyr::filter(sex != "m")

prl_hc_names <- unique(prl_hc$subj_name) 

set.seed(12345)
controls_55_names <- sample(prl_hc_names, 55)

# Select names of 52 patients, 55 controls who performed both tasks, and 24 RI
balanced_subjs_prl_task_names <- c(
  controls_55_names, patients_both_experiments, at_risk_names
)

# Get the raw data for the 52 + 55 + 24 participants
patients_controls_ri_prl_df <- 
  prl_dat[prl_dat$subj_name %in% balanced_subjs_prl_task_names, ]

length(unique(patients_controls_ri_prl_df$subj_name))

patients_controls_ri_prl_df$group <- 
  as.character(patients_controls_ri_prl_df$group)

patients_controls_ri_prl_df$group <- ifelse(
  patients_controls_ri_prl_df$group == "controls" & 
    patients_controls_ri_prl_df$diag_cat == "RI", "atrisk", 
  patients_controls_ri_prl_df$group
)


# Remove subject ch_pi_2001_10_08_418_f because Gelman-Rubin Rhat is 1.2.
clean_data <- patients_controls_ri_prl_df[
  !(patients_controls_ri_prl_df$subj_name %in% "ch_pi_2001_10_08_418_f"), 
]


# pazienti da buttare: 
bad_patients <- c(
  "ch_be_1990_12_20_153_f",
  "ca_fa_1996_03_26_092_f",
  "ma_ba_1995_05_25_321_f",
  "de_sc_1992_07_02_116_f",
  "ma_be_1997_09_01_726_f",
  "fr_mo_1982_08_15_048_f"
)

# controlli da buttare: 
bad_controls <- c(
  "ir_po_1993_02_15_077_f",
  "an_am_1993_05_20_789_f",
  "ar_mo_1999_08_01_051_f",
  "ca_co_1995_11_12_414_f",
  "ca_va_2001_08_28_797_f"
)

subjects_to_remove <- c(bad_patients, bad_controls)
clean2_data <- clean_data[!(clean_data$subj_name %in% subjects_to_remove), ]

clean2_data |> 
  group_by(diag_cat) |> 
  summarize(
    n = n_distinct(subj_name)
  )

clean2_data$group <- factor(clean2_data$group)

# At risk df ------
atrisk_df <- clean2_data |> 
  dplyr::filter(group == "atrisk")

length(unique(atrisk_df$subj_name))

# Patients df ------
patients_df <- clean2_data |> 
  dplyr::filter(group == "patients")

length(unique(patients_df$subj_name))

# HC df -------

hc_df <- prl_dat |> 
  dplyr::filter(group == "controls" & diag_cat != "RI")

length(unique(hc_df$subj_name))

hc_performance <- hc_df |> 
  group_by(subj_name) |> 
  summarise(
    acc = mean(feedback),
    n = n()
  ) 

hc_performance_good <- hc_performance |> 
  dplyr::filter(acc > 0.52 & n == 320)
  
names_hc_performance_good <- hc_performance_good$subj_name

filtered_females_names <- 
  names_hc_performance_good[grepl("f$", names_hc_performance_good)]


filtered1_hc <- hc_df[hc_df$subj_name %in%filtered_females_names, ]

hc1_performance <- filtered1_hc |> 
  group_by(subj_name) |> 
  summarise(
    acc = mean(feedback),
    n = n()
  ) 

hc_names <- unique(filtered1_hc$subj_name)

set.seed(1)
hc_to_keep <- sample(hc_names, 50)

hc_df <- prl_dat[prl_dat$subj_name %in% hc_to_keep, ]

mean(hc_df$feedback)


# ### I need to remove the HC that only completed one condition:
# hc_with_only_one_condition <- c(
#   'gi_ma_2000_06_20_421_f', 'gi_pr_2000_03_11_167_f', 'ir_bo_1995_04_03_290_f', 
#   'ma_ma_2000_12_06_306_f', 'sa_st_1998_04_21_184_f', 'so_ia_2000_03_21_569_f'
# )
# 
# # I will replace them with HC that completed both conditions and have a level
# # of performance similar to the mean of the whole group.
# 
# temp_hc <- prl_dat |> 
#   dplyr::filter(diag_cat == "HC") # 
# 
# hc_grp_avg <- mean(temp_hc$feedback, trim = 0.1)
# hc_grp_avg
# # [1] 0.5490772
# 
# temp_bysubj_fdk <- temp_hc |> 
#   group_by(subj_name) |> 
#   summarize(
#     avg_fdk = mean(feedback)
#   )
# 
# temp_bysubj_avg_fdk <- temp_bysubj_fdk[
#     temp_bysubj_fdk$avg_fdk > 0.52, 
#   ]
# 
# candidate_hc <- temp_bysubj_avg_fdk$subj_name
# candidate_hc_good <- c(
#   "al_la_2001_07_14_104_f", 
#   "ca_sa_2001_06_01_608_f",  
#   "fr_pl_2002_02_14_755_f",
#   "gi_ma_2000_06_20_421_f", "la_sa_2001_11_15_307_f",
#   "le_la_2001_11_11_647_f", "ma_ce_2002_04_17_755_f",
#   "ma_na_2001_12_03_678_f",
#   "mi_pr_1997_07_03_575_f", "mo_ru_2001_04_11_516_f", 
#   "re_ve_2001_03_28_201_f", "sa_mu_1999_09_30_649_f",
#   "sa_pa_2001_05_14_311_f", "sa_po_2001_07_11_268_f",
#   "so_sa_2001_08_07_938_f", "so_tr_2001_01_04_222_f" 
# )
# 
# # These are the 50 HC with both conditions that had been originally selected.
# original_hct50 <- c(
#   'ma_pe_2000_09_06_022_f', 'ma_sp_2000_08_01_464_f',
#   'gr_po_2002_02_12_110_f', 'co_pr_2002_03_03_902_f',
#   'el_ma_2001_07_17_978_f', 'ma_ro_1996_10_30_794_f',
#   'ma_na_2001_12_03_678_f', 'vi_te_2001_06_08_644_f',
#   'la_sc_2001_11_22_662_f', 'ma_ca_1995_10_01_691_f',
#   'el_ri_2000_10_18_444_f', 'co_mi_2001_06_21_110_f',
#   'al_zu_2001_03_12_239_f', 'fe_de_1996_03_19_136_f',
#   'va_st_2001_05_08_382_f', 'gi_sa_2001_05_25_391_f',
#   'ma_ma_2000_12_06_306_f', 'al_an_1996_06_03_205_f',
#   'fe_te_2000_05_17_086_f', 'ch_ca_1998_01_31_179_f',
#   'so_mo_2001_10_27_943_f', 'ca_ma_2000_10_22_906_f',
#   'fr_ta_1996_03_28_515_f', 'ir_ma_2001_08_24_646_f',
#   'gi_pe_2000_01_28_173_f', 'ma_ma_2001_12_21_380_f',
#   'oa_sp_1989_01_05_303_f', 'vi_ze_2001_01_02_521_f',
#   'il_ma_2001_07_30_601_f', 'ca_bu_1980_05_03_267_f',
#   'fe_ro_1999_10_25_558_f', 'el_ma_2001_12_21_179_f',
#   'el_ot_2001_09_05_187_f', 'gi_ma_2000_06_20_421_f',
#   'so_ta_1996_08_06_280_f', 'ma_pe_2000_11_24_865_f',
#   'la_sa_2001_11_15_307_f', 'al_lo_2001_02_10_286_f',
#   'ch_ma_2001_05_30_271_f', 'so_lu_2001_04_19_331_f',
#   'al_ma_2001_03_01_678_f', 'zo_me_2000_05_18_333_f',
#   'sa_pa_2001_05_14_311_f', 'so_ia_2000_03_21_569_f',
#   'ag_no_2000_02_12_330_f', 'fr_ba_1997_10_29_663_f',
#   'gi_pr_2000_03_11_167_f', 'sa_st_1998_04_21_184_f',
#   'ir_bo_1995_04_03_290_f', 'sa_li_2001_12_08_953_f'
# )
# 
# # Candidates HC that are not present in the original set
# setdiff(candidate_hc_good, original_hct50)
# 
# # set.seed(345)
# # names_hc_to_be_added <- sample(setdiff(candidate_hc_good, original_hct50), 6)
# 
# names_hc_to_be_added <- c(
#   "mo_ru_2001_04_11_516_f", "fr_pl_2002_02_14_755_f", "ma_ce_2002_04_17_755_f",
#   "ca_sa_2001_06_01_608_f", "le_la_2001_11_11_647_f", "re_ve_2001_03_28_201_f"
# )
# 
# # Selected randomly v2
# names_hc_to_be_added <- c(
#   "ag_no_2000_02_12_330_f", "al_an_1996_06_03_205_f", "al_la_2001_07_14_104_f",
#   "al_lo_2001_02_10_286_f", "al_ma_2001_03_01_678_f", "al_me_2001_06_24_456_f"
# )
# 
# # Check whether they completed both conditions (360 trials total)
# candidate_hc_df <- prl_dat[prl_dat$subj_name %in% names_hc_to_be_added, ]
# 
# candidate_hc_df |> 
#   group_by(subj_name) |> 
#   summarize(
#     n = n()
#   )
# 
# # To be added to the 44 HC
# # hc_new <- c(
# #   "mo_ru_2001_04_11_516_f", "fr_pl_2002_02_14_755_f", 
# #   "ma_ce_2002_04_17_755_f", "ca_sa_2001_06_01_608_f", 
# #   "le_la_2001_11_11_647_f", "re_ve_2001_03_28_201_f"
# # )
# hc_new <- names_hc_to_be_added

# Having found the 6 HC who completed both conditions, I remove the 6 HC that
# completed only one condition and add the 6 HC who completed both conditions.

# All data without the 6 HC who completed only one condition
# clean3_data <- clean2_data[
#   !(clean2_data$subj_name %in% hc_with_only_one_condition), ]
# 
# # Add the 6 HC that had been selected above (hc_new)
# to_be_added_df <- prl_dat[prl_dat$subj_name %in% hc_new, ]
# 
# clean4_data <- rbind(clean3_data, to_be_added_df)

clean4_data <- rbind(hc_df, atrisk_df, patients_df)

clean4_data |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_name),
    p = mean(feedback)
  )

clean4_data |> 
  group_by(group, subj_name) |> 
  summarize(
    p = mean(feedback)
  ) |> 
  as.data.frame()


saveRDS(
  clean4_data,
  here::here(
    "data", "processed", "prl", "raw_prl_data", 
    "balanced_patients_controls_risk_data.RDS")
)


# 6. Create input for HDDMrl ----------------------------------------------

source(here::here("src", "R", "functions", "funs_prl.R"))
write_input_for_hddmrl_wcst()


# ----- eof -------





