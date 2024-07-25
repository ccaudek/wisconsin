# Overview ----------------------------------------------------------------
# Associated project: WCST, PRL, and task switching
# Script purpose: Compute difference in AUC for classification with PRL 
#  and WCST tasks.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2024-07-24
# Last update: 
# Status: In progress
# Notes: 


# Load necessary libraries ------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cmdstanr)
  library(brms)
  library(workflows)
  library(tune)
  library(pROC)
  library(tidymodels)  
  library(stringr)
  library(readr)       
  library(tidyr)
  library(vip)
})


# Load and prepare data ---------------------------------------------------
data_files <- list(
  wcst = "src/stan_auc/models_params/wcst_params_classification_wide.csv",
  prl = "src/stan_auc/models_params/hddmrl_params_wide.csv"
)

data_list <- lapply(data_files, function(file) {
  tryCatch(
    rio::import(here::here(file)),
    error = function(e) {
      stop(paste("Error loading file:", file, "\n", e))
    }
  )
})

data1 <- data_list$wcst
data2 <- data_list$prl

data1$index <- NULL

# Standardize all columns except the last one (is_patient)
data1_st <- data1 %>%
  mutate(across(-is_patient, ~ scale(.) %>% as.numeric()))

data2_st <- data2 |> 
  dplyr::select(-index) |> 
  mutate(across(-is_patient, ~ scale(.) %>% as.numeric()))


# Prepare data for Stan ---------------------------------------------------
stan_data <- list(
  N1 = nrow(data1_st),
  N2 = nrow(data2_st),
  K1 = 7,  # number of predictors in data1
  K2 = 5,  # number of predictors in data2
  X1 = as.matrix(data1_st[, c("MB_Arew", "MB_Apun", "MB_gamma", "MF_Arew", "MF_Apun", "MF_gamma", "temp")]),
  X2 = as.matrix(data2_st[, c("a", "alpha", "pos_alpha", "t", "v")]),
  y1 = data1_st$is_patient,
  y2 = data2_st$is_patient
)

# Fit the model or load existing fit --------------------------------------
fit_file_path <- here::here(
  "src/stan_auc/saved_fits/fit_classification_wcst_prl.RDS"
)

if (!file.exists(fit_file_path)) {
  model_path <- here::here("src/stan_auc/auc_comparison_model_2.stan")
  model <- tryCatch(
    cmdstan_model(model_path),
    error = function(e) {
      stop(paste("Error compiling Stan model:", model_path, "\n", e))
    }
  )
  
  fit <- tryCatch(
    model$sample(
      data = stan_data,
      seed = 123,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 2000,
      iter_sampling = 2000
    ),
    error = function(e) {
      stop("Error during model fitting: ", e)
    }
  )
  
  fit$save_object(file = fit_file_path)
  rm(fit)  # Free memory
}

fit <- readRDS(fit_file_path)

# Extract and analyze results ---------------------------------------------
posterior <- fit$draws(c("auc1", "auc2", "auc_diff"), format = "df")

# Summary statistics
summary <- posterior %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarize(
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    q2.5 = quantile(value, 0.025),
    q97.5 = quantile(value, 0.975)
  )

print(summary)
# parameter      mean   median       sd    q2.5    q97.5
# <chr>         <dbl>    <dbl>    <dbl>   <dbl>    <dbl>
# 1 .chain        2.5      2.5      1.12    1        4    
# 2 .draw      4000.    4000.    2310.    201.    7800.   
# 3 .iteration 1000.    1000.     577.     51.0   1950.   
# 4 auc1          0.552    0.560    0.199   0.166    0.902
# 5 auc2          0.801    0.838    0.149   0.438    0.985
# 6 auc_diff      0.249    0.250    0.248  -0.252    0.700

# Plot posterior distributions
ggplot(posterior, aes(x = auc_diff)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Posterior Distribution of AUC Difference",
    x = "AUC(PRL) - AUC(WCST)",
    y = "Density"
  )

prob_auc_diff_positive <- mean(posterior$auc_diff > 0)
print(prob_auc_diff_positive)
