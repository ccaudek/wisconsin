# Overview ----------------------------------------------------------------
# Associated project: WCST
# Script purpose: Fit a Rescorla-Wagner model to the WCST data.
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-02-28
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

num_cores <- parallel::detectCores() # Automatically detect the number of cores
# Ensure we don't exceed the number of available cores
parallel_chains <- min(4, num_cores) # Limit to 4 chains or available cores


# Rescorla-Wagner model ---------------------------------------------------

file <- here::here(
  "src", "wcst", "documentation", "stan", "01_rescorla_wagner.stan"
)

# Compile model
mod <- cmdstan_model(file)

fit1 <- mod$sample(
  data = stan_data,
  seed = 1234,
  chains = 4, # Use 4 chains
  parallel_chains = parallel_chains, # Run chains in parallel
  threads_per_chain = 1 # Optional: Use 1 thread per chain (default)
)

fit1 <- mod$variational(
  data = stan_data,
  seed = 1234
)

# Save fit.
fit1$save_object(
  file = here::here("src", "wcst", "fits", "fit1.RDS")
)

alpha_summ <- fit1$summary("alpha")
beta_summ <- fit1$summary("beta")

mydat <- data.frame(
  subject_id = seq_len(N),

  # Convert group to a 0/1 indicator:
  # group=1 => is_patient=0, group=2 => is_patient=1
  is_patient = as.integer(stan_data$group == 2),
  alpha = alpha_summ$mean,
  beta = beta_summ$mean
)

params <- c("alpha", "beta")
params_string <- paste(params, collapse = " + ")
formula_string <- paste("is_patient ~", params_string)
formula_logit <- as.formula(formula_string)
formula_logit

# Fit the logistic regression (GLM)
fm <- glm(
  formula = formula_logit,
  family  = binomial(link = "logit"),
  data    = mydat
)

test_prob <- predict(fm, newdata = mydat, type = "response")
# Compute ROC curve and AUC using pROC
roc_obj <- roc(mydat$is_patient, test_prob, plot = TRUE, print.auc = TRUE)

t.test(beta ~ is_patient, data = mydat)



# Compute data frame for computing AUC ------------------------------------

# Summaries for alpha[i]
alpha_pos_summ <- fit1$summary("alpha_pos")
alpha_neg_summ <- fit1$summary("alpha_neg")

# Summaries for beta[i]
beta_summ <- fit1$summary("beta")

persev_summ <- fit1$summary("persev")

# Number of subjects
N <- stan_data$N

# Create a data frame with subject ID
mydat <- data.frame(
  subject_id = seq_len(N),

  # Convert group to a 0/1 indicator:
  # group=1 => is_patient=0, group=2 => is_patient=1
  is_patient = as.integer(stan_data$group == 2),

  # alpha: posterior mean (or median) for each subject
  alpha_pos = alpha_pos_summ$mean,
  alpha_neg = alpha_neg_summ$mean,

  # beta: posterior mean (or median)
  beta = beta_summ$mean,
  persev = persev_summ$mean
)

head(mydat)
summary(mydat$alpha_pos)
summary(mydat$alpha_neg)
summary(mydat$beta)
summary(mydat$persev)


# Compute AUC -------------------------------------------------------------

# Specify which parameters to include
params <- c("alpha_pos", "alpha_neg", "beta", "persev")

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

t.test(persev ~ is_patient, data = mydat)


# Rescorla-Wagner model plus inertia --------------------------------------

file <- here::here(
  "src", "wcst", "documentation", "stan", "02_rescorla_wagner_inertia.stan"
)

# Compile model
mod <- cmdstan_model(file)

fit2 <- mod$sample(
  data = stan_data,
  seed = 1234
)

# Save fit.
fit2$save_object(
  file = here::here("src", "wcst", "fits", "fit2.RDS")
)


# Models' comparisons -----------------------------------------------------

log_lik1 <- fit1$draws("log_lik", format = "matrix")
loo1 <- loo(log_lik1)

log_lik2 <- fit2$draws("log_lik", format = "matrix")
loo2 <- loo(log_lik2)

comp <- loo_compare(loo1, loo2)
print(comp, digits = 3)



# Summaries for alpha[i]
alpha_summ <- fit2$summary("alpha")
# This should have rows named alpha[1], alpha[2], ..., alpha[N]
# and columns including "mean", "median", "sd", etc.

# Summaries for beta[i]
beta_summ <- fit2$summary("beta")

# Summaries for beta[i]
inertia_summ <- fit2$summary("inertia")

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
  beta = beta_summ$mean,

  # inertia: posterior mean (or median)
  inertia = inertia_summ$mean
)

head(mydat)


# Compute AUC -------------------------------------------------------------

# Specify which parameters to include
params <- c("alpha", "beta", "inertia")

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
