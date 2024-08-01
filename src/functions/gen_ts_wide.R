# Overview ----------------------------------------------------------------
# Associated project: WCST, Task Switching, and PRL
# Script purpose: Convert in wide format the parms file
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2024-08-01
# Last update: 
# Status: In progress
# Notes: 


# Load necessary libraries ------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(rio)
})

d_long <- rio::import(
  here::here("src", "stan_auc", "models_params", "task_switching_params_long.csv")
) |> 
  dplyr::select(-c(params, subj_id))

# Convert the switch column
d_long <- d_long %>%
  mutate(switch = ifelse(switch == 1, "switch", "repeat"))

d_long$is_patient <- ifelse(d_long$group == "AN", 1, 0)
d_long$group <- NULL

# Convert data frame from long to wide format
d_wide <- pivot_wider(d_long, names_from = c(par, switch), values_from = val)

# Save CSV file
rio::export(
  d_wide, 
  here::here(
    here::here("src", "stan_auc", "models_params", "task_switching_params.csv")
  )
)


