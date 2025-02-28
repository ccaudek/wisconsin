#' Script for generating the RDS file with the behavioral statistics 
#' (perseverative errors, ...) for the three groups of participants in the
#' project.


suppressPackageStartupMessages({
  library("tidyverse")
  library("here")
  library("brms")
  library("cmdstanr")
})

source(here::here("data", "list_participants.R"))


temp1 <- readRDS(
  here::here("data", "external", "correspondence_table_psytoolkit_code.RDS")
) 

temp1$subj_name <- temp1$subj_name |> 
  dplyr::recode(
    "em_or_2003_01_01_101_f" = "em_or_2003_01_02_101_f",
    "chiara_benazzi_1990_12_20_153_f" = "ch_be_1990_12_20_153_f",
    "asia_milani_1999_03_27_854_f" = "as_mi_1999_03_27_854_f",
    "giulia_toma_1997_07_30_762_f" = "gi_to_1997_07_30_762_f",
    "chiara_gravina_1992_09_24_323_f" = "ch_gr_1992_09_24_323_f",
    "ma a_te_2001_05_31_333_m" = "ma_te_2001_05_31_333_m",
    "giulia_polichetti_2001_05_17_867_f" = "gi_po_2001_05_17_867_f",
    "cristina_cinque_1999_08_21_931_f" = "cr_ci_1999_08_21_931_f"
  )

temp2 <- readRDS(
  here::here("data", "external", "behavioral_measures_wcst_248.RDS")
) |> 
  dplyr::rename(
    "code_psytoolkit"  = "subj_name"
  )

d <- left_join(temp2, temp1, by = "code_psytoolkit")
d$subj_name <- factor(d$subj_name)

d_an <- d[d$subj_name %in% participants_list[[1]], ]
d_an$group <- "an"
d_hc <- d[d$subj_name %in% participants_list[[2]], ]
d_hc$group <- "hc"
d_ri <- d[d$subj_name %in% participants_list[[3]], ]
d_ri$group <- "ri"

d <- rbind(d_an, d_hc, d_ri)

saveRDS(
  d,
  here::here(
    "data", "processed", "wcst", "behavioral_stats.RDS"
  )
)


# Read raw WCST data

raw_data_df <- rio::import(
  here::here(
    "data", "external", "complete_raw_data.csv"
  )
)
length(unique(raw_data_df$subj_name))


# Recover raw data of the the AN group -----------------------------------------
corresp_table_an <- temp1[temp1$subj_name %in% participants_list[[1]], ]

raw_an <- 
  raw_data_df[raw_data_df$subj_name %in% corresp_table_an$code_psytoolkit, ]
length(unique(raw_an$subj_name))
raw_an$group <- "an"

# Recover raw data of the the HC group -----------------------------------------
corresp_table_hc <- temp1[temp1$subj_name %in% participants_list[[2]], ]

raw_hc <- 
  raw_data_df[raw_data_df$subj_name %in% corresp_table_hc$code_psytoolkit, ]
length(unique(raw_hc$subj_name))
raw_hc$group <- "hc"

# Recover raw data of the the RI group -----------------------------------------
corresp_table_ri <- temp1[temp1$subj_name %in% participants_list[[3]], ]

raw_ri <- 
  raw_data_df[raw_data_df$subj_name %in% corresp_table_ri$code_psytoolkit, ]
length(unique(raw_ri$subj_name))
raw_ri$group <- "ri"

# Combine the three data sets

raw_project <- rbind(raw_an, raw_hc, raw_ri)

raw_project |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_name)
  )

rio::export(
  raw_project,
  here::here(
    "data", "processed", "wcst", "raw_data_project.csv"
  )
)


# eof ----

