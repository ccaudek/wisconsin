# Overview ----------------------------------------------------------------
# Associated project: WCST, Task Switching, and PRL  
# Script purpose: Comparisons of the Steinke's mdels
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 08:43:25 2024
# Last update: Sat Jun 15 11:09:20 2024
# Status: In progress
# Notes: 


# Load necessary libraries ------------------------------------------------

suppressPackageStartupMessages({
  library("tidyverse")
  library("cmdstanr")
  library("posterior")
  library("bayesplot")
  library("recipes")
  library("MCMCpack")
  library("bayestestR")
  library("insight")
  library("here")
  library("tidyverse")
  library("cmdstanr")
  library("brms")
  library("workflows")
  library("tune")
  library("pROC")
  library("tidymodels")  
  library("stringr")
  # Helper packages
  library("readr")       # for importing data
  library("tidyr")
  library("vip") 
  library("pROC")
})



# Get the data of the two groups ------------------------------------------

stan_data <- readRDS(
  here::here(
    "data", "processed", "wcst", "stan_data_4classification.RDS"
  )
)



# Chose model -------------------------------------------------------------

# file <- file.path("src", "wcst", "models", "stan", "01_AU_rpdi_bis.stan") 
# file <- file.path("src", "wcst", "models", "stan", "02_AU_1pd1.stan") 
# file <- file.path("src", "wcst", "models", "stan", "03_MBRL.stan") 
# file <- file.path("src", "wcst", "models", "stan", "04_PRL.stan") 
# file <- file.path("src", "wcst", "models", "stan", "05_MBRL_without_inertia.stan") #
# file <- file.path("src", "wcst", "models", "stan", "06_PRL_without_inertia.stan")
file <- file.path("src", "wcst", "models", "stan", "07_PRL_weighting.stan") 
# file <- file.path("src", "wcst", "models", "stan", "07bis_PRL_weighting.stan") 
# file <- file.path("src", "wcst", "models", "stan", "08_PRL_weighting_without_inertia.stan") 

# file <- file.path("src", "wcst", "models", "stan", "revised_07_c.stan") 

# Parmeters of interest
# params_mod1 <- c("mu_r", "mu_p", "mu_d", "mu_i", "r", "p", "d", "i")
# 
# params_mod2 <- c("mu_p", "mu_d", "p", "d")
# 
# params_mod3 <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_temp",
#                  "MB_Arew", "MB_Apun", "MB_gamma","temp","log_lik","y_pred")
# 
# params_mod5 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_temp",
#                   "MB_Arew", "MB_Apun", "temp","log_lik","y_pred")  #### PRL
# 
# params_mod4 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew",
#                   "mu_MF_Apun", "mu_MF_gamma","mu_temp", "MB_Arew", "MB_Apun",
#                   "MB_gamma","MF_Arew", "MF_Apun", "MF_gamma","temp","log_lik",
#                   "y_pred")
#   
# params_mod6 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun",
#                   "mu_temp","MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun", "temp",
#                   "log_lik","y_pred")

params_mod7 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew",
                  "mu_MF_Apun", "mu_MF_gamma","mu_temp","mu_w", "MB_Arew",
                  "MB_Apun", "MB_gamma","MF_Arew", "MF_Apun", "MF_gamma","temp",
                  "w","log_lik","y_pred")

# params_mod7bis <- c(
#   "mu_MB_Arew[1]", "mu_MB_Apun[1]", "mu_MB_gamma[1]",
#   "mu_MF_Arew[1]", "mu_MF_Apun[1]", "mu_MF_gamma[1]",
#   "mu_temp[1]", "mu_w[1]",
#   "mu_MB_Arew[2]", "mu_MB_Apun[2]", "mu_MB_gamma[2]",
#   "mu_MF_Arew[2]", "mu_MF_Apun[2]", "mu_MF_gamma[2]",
#   "mu_temp[2]", "mu_w[2]",
#   "mu_MB_Arew[3]", "mu_MB_Apun[3]", "mu_MB_gamma[3]",
#   "mu_MF_Arew[3]", "mu_MF_Apun[3]", "mu_MF_gamma[3]",
#   "mu_temp[3]", "mu_w[3]",
#   "MB_Arew", "MB_Apun", "MB_gamma",
#   "MF_Arew", "MF_Apun", "MF_gamma",
#   "temp", "w",
#   "log_lik", "y_pred"
# )
# 
# params_mod8 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun",
#                   "mu_temp","mu_w", "MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun",
#                   "temp","log_lik","y_pred")


# Fix old syntax ---------------------------------------------------------------
# The original stan files are written with an old syntax.
# To fix them, follow the instructions here:
# https://mc-stan.org/cmdstanr/reference/model-method-format.html

# set compile=FALSE then call format to fix old syntax
mod <- cmdstan_model(file, compile = FALSE)
# mod$format(canonicalize = list("deprecations"))
# overwrite the original file instead of just printing it
mod$format(canonicalize = list("deprecations"), overwrite_file = TRUE)



# compile model
mod <- cmdstan_model(file)

mod$print()
# mod$exe_file()

if (0) {
  fit_mcmc <- mod$sample(
    data = stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 50 # print update every 500 iters
  )
  
  # Save the object to a file.
  qs::qsave(x = fit_mcmc, file = "fit_steinke_mod7.qs")
}


# Read the object.
fit_vi <- qs::qread(
  here::here(
    "src", "wcst", "models", "traces", "fit_steinke_mod7.qs"
  )
)

if (0) {
  # Variational inference ---
  fit_vi <- mod$variational(
    data = stan_data,
    seed = 1234
  )
}

draws <- fit_vi$draws(variables = params_mod7, format = "data.frame")
# draws <- fit_mcmc$draws(variables = params_mod7, format = "data.frame")

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
draws_long %>% tibble::as_tibble() %>% print(n=500)




# Trasforma i dati in un formato lungo per facilitare il calcolo delle medie per soggetto
draws_long <- draws %>%
  pivot_longer(cols = !starts_with("mu_"), names_to = "parameter", values_to = "value")


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
draws_long %>% tbl_df %>% print(n=200)





# fit_mcmc <- readRDS(
#   here::here("src", "traces", "fit_wcst_mcmc.rds")
# )

# Get names parameters
params <- names(mod$variables()$parameters)

fit_vi$summary(
  variables = params,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)


# Model 8, patients
# variable   mean median    sd    mad      q5    q95   q2.75  q97.5
# 1 mu_p[1]   2.48   2.43  0.492 0.463   1.77    3.37   1.67    3.61 
# 2 mu_p[2]  -3.11  -3.05  0.417 0.394  -3.86   -2.51  -4.03   -2.43 
# 3 mu_p[3]  -3.85  -3.81  0.347 0.316  -4.50   -3.37  -4.65   -3.30 
# 4 mu_p[4]  -3.13  -3.06  0.372 0.351  -3.85   -2.64  -4.00   -2.57 
# 5 mu_p[5]  -2.01  -1.99  0.108 0.0978 -2.21   -1.86  -2.27   -1.85 
# 6 mu_p[6]  -1.54  -1.52  0.140 0.133  -1.80   -1.35  -1.86   -1.32 


# Obtain a posterior mode (penalized maximum likelihood) estimate.
# fit_mle <- mod$optimize(data = data_list, seed = 123)



params_mod <- params_mod7 # <------------ Change this for each of the 8 models.

res_mcmc <- fit_mcmc$summary(params_mod) 
traces_df <- res_mcmc %>% as.data.frame()

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

traces_wide_df$is_patient <- ifelse(traces_wide_df$subj_idx < 46, 1, 0)

mydat <- traces_wide_df |> 
  dplyr::select(-subj_idx)




fit <- glm(
  is_patient ~ MB_Apun + MB_Arew + MB_gamma + MF_Apun + MF_gamma + w + temp, 
  family = binomial(),
  data = mydat
)

test_prob = predict(fit, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)




log_lik <- fit_mcmc$draws("log_lik", format = "matrix")

library(loo)

waic_result <- waic(log_lik)
print(waic_result)




# Model 3 -----------------------------------------------------------------

fit3 <- glm(
  is_patient ~ MB_Arew + MB_Apun + MB_gamma + temp, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit3, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.632


# Model 4 -----------------------------------------------------------------

fit4 <- glm(
  is_patient ~ MB_Arew + MB_Apun + MB_gamma + MF_Arew + MF_Apun + MF_gamma + temp, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit4, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.640


# Model 5  ----------------------------------------------------------------

# Rescorla-Wagner

fit5 <- glm(
  is_patient ~ MB_Apun + MB_Arew + temp, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit5, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.511


# Model 6 -----------------------------------------------------------------

fit6 <- glm(
  is_patient ~ MB_Arew + MB_Apun + MF_Arew + MF_Apun + temp, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit6, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.627


# Model 7 -----------------------------------------------------------------

params_mod
fit7 <- glm(
  is_patient ~ MB_Arew + MB_Apun + MB_gamma + MF_Arew + 
    MF_Apun + MF_gamma + temp + w, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit7, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.672


# Model 7bis --------------------------------------------------------------

params_mod
fit7bis <- glm(
  is_patient ~ MB_Arew + MB_Apun + MB_gamma + MF_Arew + MF_Apun + MF_gamma +
    temp + w, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit7bis, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: NON CONVERGE


# Model 8 -----------------------------------------------------------------

params_mod
fit8 <- glm(
  is_patient ~ MB_Arew + MB_Apun + MF_Arew + MF_Apun + temp, 
  family = binomial(),
  data = mydat
)
test_prob = predict(fit8, newdata = mydat, type = "response")
test_roc = roc(mydat$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.662







traces_clean_df <- traces_clean_df[-c(1:8), ]

# Splitting the variable column into two
traces_clean_df <- traces_clean_df %>%
  separate(variable, into = c("variable_name", "index"), sep = "\\[") %>%
  mutate(index = gsub("\\]", "", index), # Remove the closing bracket
         index = as.integer(index)) # Convert index to integer


traces_wide_df <- traces_clean_df %>%
  pivot_wider(names_from = variable_name, values_from = mean, id_cols = index)

rio::export(traces_wide_df, "wcst_params_classification.csv")



# eof ----




max(res_mcmc$rhat)
# [1] 1.011177

mean(res_mcmc$rhat)
# [1] 1.000409


# Group 1: patients
# Group 2: hc
# group 3: ri

# Extract posterior samples for mu_MB_Arew -------------------------------------
mu_MB_Arew_1_samples <- fit_mcmc$draws("mu_MB_Arew[1]")
mu_MB_Arew_2_samples <- fit_mcmc$draws("mu_MB_Arew[2]")
mu_MB_Arew_3_samples <- fit_mcmc$draws("mu_MB_Arew[3]")

mean(mu_MB_Arew_1_samples > mu_MB_Arew_2_samples)
# [1] 0.4685
mean(mu_MB_Arew_1_samples > mu_MB_Arew_3_samples)
# [1] 0.414

# Extract posterior samples for mu_MB_Apun -------------------------------------
mu_MB_Apun_1_samples <- fit_mcmc$draws("mu_MB_Apun[1]")
mu_MB_Apun_2_samples <- fit_mcmc$draws("mu_MB_Apun[2]")
mu_MB_Apun_3_samples <- fit_mcmc$draws("mu_MB_Apun[3]")

mean(mu_MB_Apun_1_samples > mu_MB_Apun_2_samples)
# [1] 0.4635
mean(mu_MB_Apun_1_samples > mu_MB_Apun_3_samples)
# [1] 0.3405

# Extract posterior samples for mu_MB_gamma ------------------------------------
mu_MB_gamma_1_samples <- fit_mcmc$draws("mu_MB_gamma[1]")
mu_MB_gamma_2_samples <- fit_mcmc$draws("mu_MB_gamma[2]")
mu_MB_gamma_3_samples <- fit_mcmc$draws("mu_MB_gamma[3]")

mean(mu_MB_gamma_1_samples > mu_MB_gamma_2_samples)
# [1] 0.67925
mean(mu_MB_gamma_1_samples > mu_MB_gamma_3_samples)
# [1] 0.36075


# Extract posterior samples for MF_Arew ----------------------------------------
mu_MF_Arew_1_samples <- fit_mcmc$draws("mu_MF_Arew[1]")
mu_MF_Arew_2_samples <- fit_mcmc$draws("mu_MF_Arew[2]")
mu_MF_Arew_3_samples <- fit_mcmc$draws("mu_MF_Arew[3]")

mean(mu_MF_Arew_1_samples > mu_MF_Arew_2_samples)
# [1] 0.484
mean(mu_MF_Arew_1_samples > mu_MF_Arew_3_samples)
# [1] 0.51225

# Extract posterior samples for MF_Arew ----------------------------------------
mu_MF_Apun_1_samples <- fit_mcmc$draws("mu_MF_Apun[1]")
mu_MF_Apun_2_samples <- fit_mcmc$draws("mu_MF_Apun[2]")
mu_MF_Apun_3_samples <- fit_mcmc$draws("mu_MF_Apun[3]")

mean(mu_MF_Apun_1_samples > mu_MF_Apun_2_samples)
# [1] 0.454
mean(mu_MF_Apun_1_samples > mu_MF_Apun_3_samples)
# [1] 0.321

# Extract posterior samples for mu_MF_gamma ------------------------------------
mu_MF_gamma_1_samples <- fit_mcmc$draws("mu_MF_gamma[1]")
mu_MF_gamma_2_samples <- fit_mcmc$draws("mu_MF_gamma[2]")
mu_MF_gamma_3_samples <- fit_mcmc$draws("mu_MF_gamma[3]")

mean(mu_MF_gamma_1_samples > mu_MF_gamma_2_samples)
# [1] 0.5575
mean(mu_MF_gamma_1_samples > mu_MF_gamma_3_samples)
# [1] 0.312

# Extract posterior samples for mu_MF_temp ------------------------------------
mu_temp_1_samples <- fit_mcmc$draws("mu_temp[1]")
mu_temp_2_samples <- fit_mcmc$draws("mu_temp[2]")
mu_temp_3_samples <- fit_mcmc$draws("mu_temp[3]")

mean(mu_temp_1_samples > mu_temp_2_samples)
# [1] 0.487
mean(mu_temp_1_samples > mu_temp_3_samples)
# [1] 0.31925

# Extract posterior samples for mu_MF_temp ------------------------------------
mu_w_1_samples <- fit_mcmc$draws("mu_w[1]")
mu_w_2_samples <- fit_mcmc$draws("mu_w[2]")
mu_w_3_samples <- fit_mcmc$draws("mu_w[3]")

mean(mu_w_1_samples > mu_w_2_samples)
# [1] 0.546
mean(mu_w_1_samples > mu_w_3_samples)
# [1] 0.50025


# eof ----


