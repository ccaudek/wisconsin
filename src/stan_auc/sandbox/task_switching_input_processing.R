library(tidyverse)

d <- rio::import(
  here::here(
    "src", "stan_auc", "sandbox", "subj_code_hddm_params.csv"
    )
  ) |> 
  dplyr::select(-c(params, subj_id))

d$is_patient <- ifelse(d$group == "AN", 1, 0)
d$group <- NULL

d_wide <- d |> 
  pivot_wider(names_from = c(par, switch), values_from = val)

# Rename columns
d_wide <- d_wide |> 
  rename_with(~ gsub("_0", "_repeat", .), contains("_0")) |> 
  rename_with(~ gsub("_1", "_switch", .), contains("_1")) |> 
  dplyr::select(-subj_idx) |> 
  select(-is_patient, everything(), is_patient)


rio::export(
  d_wide, 
  here::here(
    "src", "stan_auc", "models_params", "task_switching_params.csv"
  )
)
