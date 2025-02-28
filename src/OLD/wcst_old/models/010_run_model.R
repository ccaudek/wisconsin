# Script name: 31_run_model.R
# Project: Eating disorders Montecatini
# Script purpose: run the Stan models for the WCST data
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: date_created
# Last Modified Date: Thu Jan 18 15:10:54 2024
#
# 👉 

library("tidyverse")
library("cmdstanr")
library("posterior")
library("bayesplot")
library("recipes")
library("MCMCpack")
library("bayestestR")
library("insight")


# Get the data of the three groups.
stan_data <- readRDS(
  here::here(
    "data", "processed", "wcst", "stan_data_for_wcst_model.RDS"
  )
)


# chose model
file <- file.path("src", "wcst", "models", "stan", "01_AU_rpdi_bis.stan") 
file <- file.path("src", "wcst", "models", "stan", "02_AU_1pd1.stan") 
file <- file.path("src", "wcst", "models", "stan", "03_MBRL.stan") 
file <- file.path("src", "wcst", "models", "stan", "04_PRL.stan") 
file <- file.path("src", "wcst", "models", "stan", "05_MBRL_without_inertia.stan") #
file <- file.path("src", "wcst", "models", "stan", "06_PRL_without_inertia.stan")
file <- file.path("src", "wcst", "models", "stan", "07_PRL_weighting.stan") 
file <- file.path("src", "wcst", "models", "stan", "07bis_PRL_weighting.stan") 
file <- file.path("src", "wcst", "models", "stan", "08_PRL_weighting_without_inertia.stan") 


# Parmeters of interest
params_mod1 <- c("mu_r", "mu_p", "mu_d", "mu_i", "r", "p", "d", "i")

params_mod2 <- c("mu_p", "mu_d", "p", "d")

params_mod3 <- c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_temp",
                 "MB_Arew", "MB_Apun", "MB_gamma","temp","log_lik","y_pred")

params_mod5 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_temp",
                  "MB_Arew", "MB_Apun", "temp","log_lik","y_pred")  #### PRL

params_mod4 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew", 
                  "mu_MF_Apun", "mu_MF_gamma","mu_temp", "MB_Arew", "MB_Apun", 
                  "MB_gamma","MF_Arew", "MF_Apun", "MF_gamma","temp","log_lik",
                  "y_pred")
  
params_mod6 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun", 
                  "mu_temp","MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun", "temp",
                  "log_lik","y_pred")

params_mod7 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MB_gamma", "mu_MF_Arew", 
                  "mu_MF_Apun", "mu_MF_gamma","mu_temp","mu_w", "MB_Arew", 
                  "MB_Apun", "MB_gamma","MF_Arew", "MF_Apun", "MF_gamma","temp",
                  "w","log_lik","y_pred")

params_mod7bis <- c(
  "mu_MB_Arew[1]", "mu_MB_Apun[1]", "mu_MB_gamma[1]", 
  "mu_MF_Arew[1]", "mu_MF_Apun[1]", "mu_MF_gamma[1]", 
  "mu_temp[1]", "mu_w[1]",
  "mu_MB_Arew[2]", "mu_MB_Apun[2]", "mu_MB_gamma[2]", 
  "mu_MF_Arew[2]", "mu_MF_Apun[2]", "mu_MF_gamma[2]", 
  "mu_temp[2]", "mu_w[2]",
  "mu_MB_Arew[3]", "mu_MB_Apun[3]", "mu_MB_gamma[3]", 
  "mu_MF_Arew[3]", "mu_MF_Apun[3]", "mu_MF_gamma[3]", 
  "mu_temp[3]", "mu_w[3]",
  "MB_Arew", "MB_Apun", "MB_gamma", 
  "MF_Arew", "MF_Apun", "MF_gamma", 
  "temp", "w", 
  "log_lik", "y_pred"
)

params_mod8 <-  c("mu_MB_Arew", "mu_MB_Apun", "mu_MF_Arew", "mu_MF_Apun",
                  "mu_temp","mu_w", "MB_Arew", "MB_Apun", "MF_Arew", "MF_Apun", 
                  "temp","log_lik","y_pred")


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

fit_mcmc <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 50 # print update every 500 iters
)

draws <- fit_mcmc$draws(variables = params_mod7, format = "data.frame")


# Trasforma i dati in un formato lungo per facilitare il calcolo delle medie per soggetto
draws_long <- draws %>%
  pivot_longer(cols = starts_with("mu_"), names_to = "parameter", values_to = "value")


# Estrai l'indice del soggetto da ogni parametro e calcola la media
draws_long <- draws_long %>%
  mutate(subject = stringr::str_extract(parameter, "\\[\\d+\\]") %>% stringr::str_remove_all("\\[|\\]"),
         parameter = stringr::str_remove(parameter, "\\[\\d+\\]")) %>%
  group_by(subject, parameter) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')

print(draws_long)

dim(draws_long)


# fit_mcmc <- readRDS(
#   here::here("src", "traces", "fit_wcst_mcmc.rds")
# )

# Get names parameters
params <- names(mod$variables()$parameters)

fit_mcmc$summary(
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



params_mod <- params_mod7bis # <------------ Change this for each of the 8 models.

res_mcmc <- fit_mcmc$summary(params_mod) 
traces_df <- res_mcmc %>% as.data.frame()

traces_df <- traces_df[!grepl("y_pred", traces_df$variable), ]
traces_df <- traces_df[!grepl("log_lik", traces_df$variable), ]

# variable         mean       median           sd          mad           q5
# 1    mu_MB_Arew[1] 0.5105776908 0.4641085000 0.1998476209 0.1464727257 2.605691e-01
# 2    mu_MB_Apun[1] 0.0003968613 0.0002108090 0.0005468703 0.0002527072 7.058052e-06
# 3   mu_MB_gamma[1] 0.0075089395 0.0047184000 0.0082852409 0.0053919197 1.760638e-04
# 4    mu_MF_Arew[1] 0.6432045364 0.7198755000 0.2878466412 0.3212015835 1.394712e-01
# 5    mu_MF_Apun[1] 0.0125371270 0.0065866050 0.0167774477 0.0079630520 2.399984e-04
# 6   mu_MF_gamma[1] 0.0137338375 0.0080181000 0.0241200730 0.0091742991 3.718174e-04
# 7       mu_temp[1] 0.1389853274 0.1406805000 0.0301778367 0.0302539356 8.654822e-02
# 8          mu_w[1] 0.7850709997 0.8608895000 0.1709890304 0.0722952825 4.000340e-01
# 9    mu_MB_Arew[2] 0.5242867507 0.4840245000 0.1953239783 0.1543557099 2.657629e-01
# 10   mu_MB_Apun[2] 0.0004822560 0.0002587975 0.0006345636 0.0003156119 8.138855e-06
# 11  mu_MB_gamma[2] 0.0032982523 0.0020184800 0.0037649853 0.0023592844 9.141218e-05
# 12   mu_MF_Arew[2] 0.6652723994 0.7389945000 0.2755105119 0.2980493019 1.738048e-01
# 13   mu_MF_Apun[2] 0.0155013118 0.0084592100 0.0193684039 0.0100929329 2.855893e-04
# 14  mu_MF_gamma[2] 0.0100660345 0.0061198650 0.0115031015 0.0070859014 2.810295e-04
# 15      mu_temp[2] 0.1405409471 0.1429055000 0.0305855397 0.0291953592 8.663962e-02
# 16         mu_w[2] 0.7798119845 0.8468265000 0.1578052140 0.0731159016 4.196666e-01
# 17   mu_MB_Arew[3] 0.5630100865 0.5292580000 0.1963316975 0.1853235174 2.853224e-01
# 18   mu_MB_Apun[3] 0.0009257392 0.0004968545 0.0012349995 0.0006062625 1.325006e-05
# 19  mu_MB_gamma[3] 0.0135295286 0.0090988300 0.0140363626 0.0100931776 4.364699e-04
# 20   mu_MF_Arew[3] 0.6354909986 0.6839335000 0.2681009838 0.3262068411 1.674330e-01
# 21   mu_MF_Apun[3] 0.0501803092 0.0174951000 0.1125359649 0.0215191680 6.879350e-04
# 22  mu_MF_gamma[3] 0.0655894705 0.0203012000 0.1432084921 0.0238367313 9.126665e-04
# 23      mu_temp[3] 0.1608304775 0.1612600000 0.0366822443 0.0369211878 9.679360e-02
# 24         mu_w[3] 0.7985933822 0.8452800000 0.1470469785 0.0923652387 4.728481e-01

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


