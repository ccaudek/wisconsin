# Overview ----------------------------------------------------------------
# Associated project: WCST
# Script purpose: Fit to the WCST data a Rescorla-Wagner model.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-02-28
# Last update:
# Status: In progress
# Notes:
# rule_choice : which category is rewarded.
#   color  : 1
#   shape  : 2
#   number : 3
#
# resp_choice rew resp_color resp_shape resp_number


# Load necessary libraries ------------------------------------------------

if (!requireNamespace("pacman")) install.packages("pacman")

pacman::p_load(
  here, tidyverse, cmdstanr, posterior, bayesplot, insight,
  pROC, stringr, loo
)


# Source functions --------------------------------------------------------

source(here::here(
  "src", "wcst", "documentation", "functions", "funs_wcst.R"
))

source(here::here(
  "src", "wcst", "documentation", "functions", "funs_model_selection_wcst.R"
))

# process_and_prepare_stan_data()
source(
  here::here(
    "src", "wcst", "documentation", "functions",
    "funs_input_for_stan_wcst.R"
  )
)

# Data preparation ----------------------------------------------------------

generate_csv("controls")
generate_csv("patients")

# Generate the input for the Stan models
stan_data <- process_and_prepare_stan_data()
str(stan_data)


# Compiling ---------------------------------------------------------------

file <- here::here(
  "src", "wcst", "documentation", "stan", "01_rescorla_wagner.stan"
)

# Compile model
mod <- cmdstan_model(file)

fit1 <- mod$sample(
  data = stan_data,
  seed = 1234
)

# Change the name of the output file!!
fit1$save_object(
  file = here::here("src", "wcst", "fits", "fit1.RDS")
)


# Compute data frame for AUC ----------------------------------------------

# Summaries for alpha[i]
alpha_summ <- fit1$summary("alpha")
# This should have rows named alpha[1], alpha[2], ..., alpha[N]
# and columns including "mean", "median", "sd", etc.

# Summaries for beta[i]
beta_summ <- fit1$summary("beta")

# Number of subjects
N <- stan_data$N

# Create a data frame with subject ID
mydat <- data.frame(
  subject_id = seq_len(N),

  # Convert group to a 0/1 indicator:
  # group=1 => is_patient=0, group=2 => is_patient=1
  is_patient = as.integer(stan_data$group == 2),

  # alpha: posterior mean (or median) for each subject
  alpha = alpha_summ$mean,

  # beta: posterior mean (or median)
  beta = beta_summ$mean
)

head(mydat)


# Compute AUC -------------------------------------------------------------

# Specify which parameters to include
params <- c("alpha", "beta")

# Build a logistic regression formula on the fly

# Join the params into "alpha + beta + ..."
params_string <- paste(params, collapse = " + ")
# Then paste into a formula "is_patient ~ alpha + beta"
formula_string <- paste("is_patient ~", params_string)
formula_logit <- as.formula(formula_string)
formula_logit

# Fit the logistic regression (GLM)
fm <- glm(
  formula = formula_logit,
  family  = binomial(link = "logit"),
  data    = mydat
)

# Predict probabilities and compute AUC

# Predict the probability that is_patient = 1 for each row of mydat
test_prob <- predict(fm, newdata = mydat, type = "response")

# Compute ROC curve and AUC using pROC
roc_obj <- roc(mydat$is_patient, test_prob, plot = TRUE, print.auc = TRUE)

# If you want to store the numeric AUC value separately:
auc_value <- roc_obj$auc
cat("AUC =", auc_value, "\n")



## eof
