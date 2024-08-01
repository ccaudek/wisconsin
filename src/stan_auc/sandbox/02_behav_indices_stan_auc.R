# Overview ----------------------------------------------------------------
# Associated project: WCST, PRL, and Task Switching for discriminating
#  AN patients and healthy controls.
# Script purpose: Compute two posterior AUC differences:
#  - AUC1 = AUC(prl) - AUC(steinke)
#  - AUC2 = AUC(steinke) - AUC(wcst_behav_indices)
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: Thu Jul 25 15:24:54 2024
# Last update:
# Pre-requisites: In "src/stan_auc/models_params" there must be 3 files:
#  - wcst_steinke_params.csv
#  - prl_params.csv
#  - wcst_behav_indices.csv
# Status: In progress
# Notes:
#  TODO: Add the Task Switching DDM parameters and change the Stan
#  script to include this condition.


# Set up ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(cmdstanr)
  library(pROC)
  library(missRanger)
})

# Find the root of the project
proj_root <- find_root(has_file("wisconsin.Rproj"))

# Load and prepare data ---------------------------------------------------
# 1: WCST
data1 <- rio::import(
  file.path(
    proj_root, "src", "stan_auc", "models_params", "wcst_steinke_params.csv"
  )
) |>
  dplyr::select(-index)

# 2: PRL
data2 <- rio::import(
  file.path(
    proj_root, "src", "stan_auc", "models_params", "prl_params.csv"
  )
) |>
  dplyr::select(-index)

# 3: behavioral
data3 <- rio::import(
  file.path(
    proj_root, "src", "stan_auc", "models_params", "wcst_behav_indices.csv"
  )
) |>
  mutate(is_patient = ifelse(group == "an", 1, 0)) |>
  dplyr::select(-c(subj_name, group))

# Impute missing values
data3 <- missRanger(data3, pmm.k = 5, num.trees = 100)

# Standardize all columns except the last one (is_patient)
data1_st <- data1 |>
  mutate(across(-is_patient, ~ scale(.) |> as.numeric()))

data2_st <- data2 |>
  mutate(across(-is_patient, ~ scale(.) |> as.numeric()))

data3_st <- data3 |>
  mutate(across(-is_patient, ~ scale(.) |> as.numeric()))

# Prepare data for Stan ---------------------------------------------------
stan_data <- list(
  N1 = nrow(data1_st),
  N2 = nrow(data2_st),
  N3 = nrow(data3_st),
  K1 = 8, # number of predictors in model 1
  K2 = 5, # number of predictors in model 2
  K3 = 3, # number of predictors in model 3
  X1 = as.matrix(
    data1_st[, c(
      "MB_Arew", "MB_Apun", "MB_gamma", "MF_Arew", "MF_Apun", 
      "MF_gamma", "temp", "w"
    )]
  ),
  X2 = as.matrix(data2_st[, c("a", "alpha", "pos_alpha", "t", "v")]),
  X3 = as.matrix(
    data3_st[, c("prop_pers_err", "prop_non_pers_err", "prop_pers_resp")]
  ),
  y1 = data1_st$is_patient,
  y2 = data2_st$is_patient,
  y3 = data3_st$is_patient
)

# Compile the Stan model --------------------------------------------------
model_path <- here::here("src", "stan_auc", "auc_3model_comparison.stan")
model <- cmdstan_model(model_path)

# Fit the model -----------------------------------------------------------
fit <- model$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 20000,
  iter_sampling = 20000
)

# Extract and analyze results ---------------------------------------------
posterior <- fit$draws(
  c("auc1", "auc2", "auc3", "auc_diff1", "auc_diff2"),
  format = "df"
)

# Summary statistics
summary <- posterior |>
  pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
  group_by(parameter) |>
  summarize(
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975)
  )

print(summary)

# Plot posterior distributions
ggplot(posterior, aes(x = auc_diff1)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior Distribution of AUC Difference",
    x = "AUC(PRL) - AUC(WCST)",
    y = "Density"
  )

ggplot(posterior, aes(x = auc_diff2)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior Distribution of AUC Difference",
    x = "AUC(PRL) - AUC(WCST)",
    y = "Density"
  )

prob_auc_diff1_positive <- mean(posterior$auc_diff1 > 0)
print(prob_auc_diff1_positive)

prob_auc_diff2_positive <- mean(posterior$auc_diff2 > 0)
print(prob_auc_diff2_positive)


## eof
