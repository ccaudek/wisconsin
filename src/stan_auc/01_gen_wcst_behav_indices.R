# Overview ----------------------------------------------------------------
# Associated project:  
# Script purpose: 
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2024-07-25
# Last update: 
# Status: In progress
# Notes: 


# Load necessary libraries ------------------------------------------------

suppressPackageStartupMessages({
  library("here")
  library("tidyverse")
})



# Perseverative errors ----------------------------------------------------
pers_err_df <- readRDS(
  here("data", "processed", "wcst", "behavioral_stats.rds")
) |> 
  dplyr::select(
    code_psytoolkit, subj_name, group,
    prop_pers_err, prop_non_pers_err
  )


# Perseverative responses -------------------------------------------------
pers_resp_df <- rio::import(
  here::here(
    "data", "processed", "wcst", "raw_data_project.csv"
  )
) 

length(unique(raw_df$subj_name))

raw_df$group <- factor(raw_df$group)
raw_df$group <- relevel(raw_df$group, ref = "hc")

raw_df$is_correct <- ifelse(
  raw_df$is_correct == 2, 0, 
  ifelse(raw_df$is_correct == 1, 1, NA)
)

perseverative_df <- raw_df %>%
  group_by(subj_name) %>%
  mutate(
    prev_name_of_task = lag(name_of_task),
    rule_changed = if_else(name_of_task != prev_name_of_task & !is.na(prev_name_of_task), TRUE, FALSE),
    post_rule_change_period = cumsum(rule_changed),
    trial_num = row_number()  # Add a column for trial number within each subject
  ) %>%
  mutate(
    prev_chosen_card = lag(chosen_card),
    is_pers_resp = if_else(
      chosen_card == prev_chosen_card & post_rule_change_period == lag(post_rule_change_period) & 
        !is.na(prev_chosen_card), 1, 0, missing = 0),
    is_pers_err = if_else(
      is_pers_resp == 1 & is_correct == 0 & trial_num > 10, 1, 0, missing = 0)
  ) %>%
  ungroup()

# Calculate proportion of perseverative responses for each participant
proportion_perseverative_df <- perseverative_df %>%
  group_by(group, subj_name) %>%
  summarize(
    total_responses = n(),
    pers_resp = sum(is_pers_resp, na.rm = TRUE),
    pers_err  = sum(is_pers_err, na.rm = TRUE),
    prop_pers_resp = pers_resp / total_responses,
    prop_pers_err = pers_err / 50,
  ) %>%
  ungroup() |> 
  dplyr::select(group, subj_name, prop_pers_resp) |> 
  dplyr::rename(code_psytoolkit = subj_name) |> 
  dplyr::select(-group)


# Join --------------------------------------------------------------------
behav_indices_df <- left_join(
  pers_err_df, proportion_perseverative_df, by = "code_psytoolkit"
) |> 
  dplyr::select(-code_psytoolkit) |> 
  dplyr::filter(group != "ri")


# Save CSV file -----------------------------------------------------------
rio::export(
  behav_indices_df,
  here::here(
    "data", "processed", "wcst", "classification", "wcst_behav_indices.csv"
  )
)



