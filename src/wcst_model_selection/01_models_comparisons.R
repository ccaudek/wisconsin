# Overview ----------------------------------------------------------------
# Associated project:  WCST, PRL, and Task Switching
# Script purpose: preprocessing raw data and generating input for Stan models
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2024-08-02
# Last update: 
# Status: In progress
# Notes: 
#
# rule_choice : which category is rewarded.
#   color  : 1
#   shape  : 2
#   number : 3
#
# resp_choice rew resp_color resp_shape resp_number


# Load necessary libraries.

suppressPackageStartupMessages({
  library("tidyverse")
  library("cmdstanr")
  library("posterior")
  library("bayesplot")
  library("insight")
  library("here")
  library("pROC")
  library("stringr")
  library("pROC")
  library("loo")
})


# ---- Source functions ---- #

source(here::here(
  "src", "wcst_model_selection", "documentation", "functions", "funs_wcst.R")
)

source(here::here(
  "src", "wcst_model_selection", "documentation", "functions", 
  "funs_model_selection_wcst.R")
)

# process_and_prepare_stan_data()
source(
  here::here("src", "wcst_model_selection", "documentation", "functions", 
             "funs_input_for_stan_wcst.R")
  )


# ---- Data wrangling ---- #

generate_csv("controls")
generate_csv("patients")

# Generate the input for the Stan models
stan_data <- process_and_prepare_stan_data()


# ---- Chose model ---- #

# Define the flag for model selection
model_flag <- 1 # Change this value to choose a different model

# Initialize variables for file path and parameters
file <- ""
params_mod <- c()

# Select model and parameters based on the flag
switch(as.character(model_flag),
       '1' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "03_MBRL.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_temp",
                     "MB_Arew", "MB_Apun", "MB_gamma", "temp", "log_lik", "y_pred")
       },
       '2' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "04_PRL.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew",
                     "mu_MF_Apun", "mu_MF_gamma", "mu_temp", "MB_Arew", "MB_Apun",
                     "MB_gamma", "MF_Arew", "MF_Apun", "MF_gamma", "temp", "log_lik",
                     "y_pred")
       },
       '3' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "05_MBRL_without_inertia.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_temp",
                     "MB_Arew", "MB_Apun", "temp", "log_lik", "y_pred")
       },
       '4' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "06_PRL_without_inertia.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun",
                     "mu_temp", "MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun", "temp",
                     "log_lik", "y_pred")
       },
       '5' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "07_PRL_weighting.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew",
                     "mu_MF_Apun", "mu_MF_gamma", "mu_temp", "mu_w", "MB_Arew",
                     "MB_Apun", "MB_gamma", "MF_Arew", "MF_Apun", "MF_gamma", "temp",
                     "w", "log_lik", "y_pred")
       },
       '6' = {
         file <- file.path("src", "wcst_model_selection", "documentation", "stan", "08_PRL_weighting_without_inertia.stan")
         params_mod <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun",
                     "mu_temp", "mu_w", "MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun",
                     "temp", "log_lik", "y_pred")
       },
       {
         stop("Invalid model flag. Please choose a number between 1 and 6.")
       }
)

# Print chosen model and parameters
print(paste("Chosen model file:", file))
print("Parameters:")
print(params_mod)

# Fix old syntax.
# The original stan files are written with an old syntax.
# To fix them, follow the instructions here:
# https://mc-stan.org/cmdstanr/reference/model-method-format.html

# set compile=FALSE then call format to fix old syntax
mod <- cmdstan_model(file, compile = FALSE)
# mod$format(canonicalize = list("deprecations"))
# overwrite the original file instead of just printing it
mod$format(canonicalize = list("deprecations"), overwrite_file = TRUE)

# Compile model.
mod <- cmdstan_model(file)

# Sampling.
if (0) {
  fit_mcmc <- mod$sample(
    data = stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 50 # print update every 500 iters
  )
}

# Variational inference
fit <- mod$pathfinder(
  data = stan_data,
  seed = 1234
)


# ---- Save fit ---- #

# Change the name of the output file!!
fit$save_object(
  file = here::here("src", "wcst_model_selection", "data", "fits", "fit1.RDS")
)

fit <- readRDS(
  here::here("src", "wcst_model_selection", "data", "fits", "fit7_OLD.RDS")
)

# ---- Get individual parameters ---- #

draws <- fit$draws(variables = params_mod, format = "data.frame")

# Trasforma i dati in un formato lungo per facilitare il calcolo delle medie per soggetto
draws_long <- draws %>%
  pivot_longer(cols = starts_with("mu_"), names_to = "parameter", values_to = "value")

# Estrai l'indice del soggetto da ogni parametro e calcola la media
draws_long <- draws_long %>%
  mutate(subject = str_extract(parameter, "\\[\\d+\\]") %>% str_remove_all("\\[|\\]"),
         parameter = str_remove(parameter, "\\[\\d+\\]")) %>%
  group_by(subject, parameter) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            sd_value = sd(value, na.rm = TRUE),
            lower_95 = quantile(value, 0.025, na.rm = TRUE),
            upper_95 = quantile(value, 0.975, na.rm = TRUE),
            .groups = 'drop')

# View the result
print(draws_long)

# Get names parameters
params <- names(mod$variables()$parameters)

fit$summary(
  variables = params,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)

# Obtain a posterior mode (penalized maximum likelihood) estimate.
# fit_mle <- mod$optimize(data = data_list, seed = 123)

res_mcmc <- fit$summary(params) 
traces_df <- res_mcmc %>% 
  as.data.frame()

traces_df <- traces_df[!grepl("y_pred", traces_df$variable), ]
traces_df <- traces_df[!grepl("log_lik", traces_df$variable), ]

traces_clean_df <- traces_df |> 
  dplyr::select(variable, mean)

# Use the separate function to split the variable column
traces_clean_df <- traces_clean_df %>%
  separate(variable, into = c("params", "subj_idx"), sep = "\\[", convert = TRUE) %>%
  mutate(subj_idx = as.numeric(gsub("\\]", "", subj_idx)))

# Assuming traces_clean_df is already separated into params and subj_idx
traces_clean_df <- traces_clean_df %>%
  mutate(subj_idx = as.integer(subj_idx)) # Ensure subj_idx is integer

# Create a complete dataset with all combinations of subj_idx and params
traces_complete_df <- traces_clean_df %>%
  complete(subj_idx, params, fill = list(mean = 0)) # Fill missing mean values with 0 or another appropriate value

# Now pivot to wide format
traces_wide_df <- traces_complete_df %>%
  pivot_wider(names_from = params, values_from = mean)

# traces_wide_df <- traces_wide_df[-c(92:115), ]

# The first 46 subjects are patients
traces_wide_df$is_patient <- ifelse(traces_wide_df$subj_idx < 46, 1, 0)

mydat <- traces_wide_df |> 
  dplyr::select(-subj_idx)

# Convert the vector to a single string with parameters separated by '+'
params_string <- paste(params, collapse = " + ")

# Create the formula string for the glm function
formula_string <- paste("is_patient ~", params_string)

# Convert the string to a formula
formula <- as.formula(formula_string)

# Build the glm model using the formula
fm <- glm(formula, family = binomial(), data = mydat)

test_prob = predict(fm, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)

# Compute WAIC -----------------------------------------------------------------
log_lik <- fit$draws("log_lik", format = "matrix")
waic_result <- waic(log_lik)
print(waic_result)


# MODEL 1:
#           Estimate   SE
# elpd_waic  -2492.3 32.1
# p_waic         0.7  0.3
# waic        4984.7 64.3

# glm AUC = 0.629

# MODEL 7:
#           Estimate   SE
# elpd_waic  -2488.8 32.4
# p_waic         6.0  1.0
# waic        4977.5 64.7

# glm AUC = 0.725




## eof

