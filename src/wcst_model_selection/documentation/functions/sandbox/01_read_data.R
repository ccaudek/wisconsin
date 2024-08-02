# Script name: 30_recoding_for_stan.R
# Project: Eating disorders Montecatini
# Script purpose: Generate list for PRL Stan model.
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed May 25 05:57:16 2022
# Last Modified Date: Thu Jan 18 15:09:18 2024
#
# 👉 This script generates the input data for Steinke's algorithm of the WCST.

# rule_choice : which category is rewarded.
#   color  : 1
#   shape  : 2
#   number : 3

# resp_choice rew resp_color resp_shape resp_number

# Prelims
library("here")
library("tidyverse")
library("stringi")

# Source functions -------------------------------------------------------------
source(here::here("src", "functions", "funs_wcst.R"))

# - get_one_subj_data_for_stan()
source(here::here("src", "functions", "funs_input_for_stan_wcst.R"))


# Generate RDS raw data for patients or controls -------------------------------

# Select group.
GROUP <- "controls"  # "controls" "patients"

dir <- here("data", "raw", GROUP)

if (GROUP == "patients") {
  file_names <- as.character(list.files(path=dir, pattern="wcst_pazienti"))
} else {
  file_names <- as.character(list.files(path=dir, pattern="wcst_eds1"))
}

n_files <- length(file_names)
n_files

d_list <- list()

for (i in 1:n_files) {
  
  d  <- read.table(
    here("data", "raw", GROUP, file_names[i]), 
    header = FALSE
  )
  
  d$subj_name <- file_names[i]
  d$block <- rep(1:6, each = 10)
  
  d$card_shown <- d$V1
  d$correct_card <- d$V2
  d$card_chosen_if_perseveration <- d$V3
  d$trial_in_a_sequence <- d$V4 
  d$name_of_task <- d$V5
  d$card_shape <- d$V6
  d$card_number_of_symbols <- d$V7
  d$card_color <- d$V8
  d$rt <- d$V9
  d$is_correct <- d$V10 # 1=correct, 2=wrong card, 3=too slow
  d$chosen_card <- d$V11 # a number between 1 and 4, or 0 if none clicked
  d$is_error <- d$V12 
  # If 1, this trial was an error (otherwise 0); 0 or 1
  d$is_perseverative_error <- d$V13 
  # If 1, this trial was a perseveration error (otherwise 0); 0, 1
  d$is_non_perseverative_error <- d$V14
  # If 1, this trial was not a perseveration error (otherwise 0); 0, 1
  
  d_list[[i]] <- d
}

# convert list into data.frame
df0 <- do.call(rbind.data.frame, d_list)

length(unique(df0$subj_name))

# Examine accuracy.
bysubj_acc <- df0 %>% 
  group_by(subj_name) %>% 
  summarise(
    error_rate = mean(is_error)
  ) %>% as.data.frame()

bad_subj_df <- bysubj_acc %>%
  dplyr::filter(error_rate > 0.5)
 
bad_subjects <- bad_subj_df$subj_name

# # Remove bad subjects.
df1 <- df0[!(df0$subj_name %in% bad_subjects), ]
length(unique(df1$subj_name))

df <- df1[!is.na(df1$chosen_card), ]

# Recoding as required by Steinke's algorithm. 
df <- df %>% 
  dplyr::rename(
    card_number = card_number_of_symbols
  )

length(unique(df$subj_name))

if (GROUP == "patients") {
  rio::export(
    df,
    here::here("src", "wcst_model_selection", "data", "raw_data_patients.csv")
  )
} else {
  rio::export(
    df,
    here::here("src", "wcst_model_selection", "data", "raw_data_controls.csv")
  )
}


# eof ----

