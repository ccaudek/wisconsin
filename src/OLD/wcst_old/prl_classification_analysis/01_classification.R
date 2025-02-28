# Script name: 017_classification.R
# Project: Eating disorders Montecatini.
# Script purpose: Use the hDDM parameters to classify patients/controls.
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed Oct 27 06:40:19 2021
# Last Modified Date: Fri Jul 14 10:23:36 2023
# 
# Notes: 
# https://www.rebeccabarter.com/blog/2020-03-25_machine_learning/
# When the individual parameters of hDDMrl are computed without
# the knowledge of the group (NO split_by = is_patient), there is
# no shrinkage that makes more similar to each other the participants'
# individual parameters of each group. In these conditions, by using
# the parameters computed in the food condition, no ability for clinical 
# classification emerges: AUC = 0.55.


# Prelims -----------------------------------------------------------------

suppressPackageStartupMessages({
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


# Increase max print
options(max.print = .Machine$integer.max)

cores <- parallel::detectCores()
cores

dat <- rio::import(
  here::here(
    "src", "wcst", "prl_classification_analysis", "hddmrl_params.csv"
  )
)

dat$is_patient <- ifelse(dat$subj_idx < 51, 0, 1)
dat$is_patient <- factor(dat$is_patient)


dd <- dat
names(dd)

# Step 1: Remove the string "_subj" from the 'knode_name' column
dd <- dd %>%
  mutate(knode_name = gsub("_subj", "", knode_name))

dd$node <- NULL

dd <- dd %>%
  group_by(knode_name, subj_idx) %>%
  summarise(mean = mean(mean), .groups = 'drop') |> 
  dplyr::rename(param = knode_name)

# Step 2: Convert the data from long to wide format
dd_wide <- dd %>%
  pivot_wider(names_from = param, values_from = mean)

dd_wide$is_patient <- ifelse(dd_wide$subj_idx < 51, 0, 1)
dd_wide$is_patient |> table()

for_saving <- dd_wide |> 
  # dplyr::select(-c(is_patient)) |> 
  dplyr::rename(index = subj_idx)

rio::export(
  for_saving, "hddmrl_params_wide.csv"
)

model_glm = glm(
  is_patient ~ a + alpha + pos_alpha + t + v, 
  data = dd_wide, 
  family = "binomial"
)

test_prob = predict(model_glm, newdata = dd_wide, type = "response")
test_roc = roc(dd_wide$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)
# AUC: 0.822


# Bayesian AUC ------------------------------------------------------------

# Compile the model
model <- cmdstan_model(
  here::here(
    "src", "wcst", "prl_classification_analysis", "prl_logistic_model.stan")
  )

# Prepare data
stan_data <- list(
  N = nrow(dd_wide),
  P = 5,  # number of predictors
  X = as.matrix(dd_wide[, c("a", "alpha", "pos_alpha", "t", "v")]),
  y = dd_wide$is_patient
)

# Fit the model
fit <- model$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500
)

library(pROC)
library(tidybayes)

# Extract posterior samples
posterior_samples <- fit$draws(format = "df")

# Print the first few column names of posterior_samples
print(head(colnames(posterior_samples), 10))

# Print the first row of posterior_samples
print(posterior_samples[1, 1:10])

calculate_auc <- function(sample) {
  # Extract beta coefficients
  beta_cols <- paste0("beta[", 1:5, "]")  # Assuming 5 predictors
  betas <- as.numeric(sample[beta_cols])
  
  # Extract alpha (intercept)
  alpha <- as.numeric(sample["alpha"])
  
  # Ensure betas is a column vector
  betas <- matrix(betas, ncol = 1)
  
  # Calculate logits
  logits <- stan_data$X %*% betas + alpha
  
  # Convert to probabilities
  probs <- 1 / (1 + exp(-logits))
  
  # Ensure probs is numeric
  probs <- as.numeric(probs)
  
  # Calculate AUC
  auc <- roc(stan_data$y, probs)$auc
  return(auc)
}


# Try to calculate AUC for the first sample
test_auc <- calculate_auc(posterior_samples[1, ])
print(test_auc)


# If the above works, calculate AUC for all samples
aucs <- apply(posterior_samples, 1, calculate_auc)

# Calculate the 0.89 credibility interval for AUCs
ci_lower <- quantile(aucs, probs = 0.055)
median <- quantile(aucs, probs = 0.5)
ci_upper <- quantile(aucs, probs = 0.945)

# Print the credibility interval
print(paste("0.89 credibility interval for AUCs:", ci_lower, "-", median, "-", ci_upper))



# Comparisons between the two groups -------------------------------------------

dd_wide %>%
  group_by(is_patient) %>%
  summarise(
    MB_Arew_mean = mean(MB_Arew),
    MB_Arew_se = sd(MB_Arew) / sqrt(n()),
    MB_Apun_mean = mean(MB_Apun),
    MB_Apun_se = sd(MB_Apun) / sqrt(n()),
    MB_gamma_mean = mean(MB_gamma),
    MB_gamma_se = sd(MB_gamma) / sqrt(n()),
    MF_Arew_mean = mean(MF_Arew),
    MF_Arew_se = sd(MF_Arew) / sqrt(n()),
    MF_Apun_mean = mean(MF_Apun),
    MF_Apun_se = sd(MF_Apun) / sqrt(n()),
    MF_gamma_mean = mean(MF_gamma),
    MF_gamma_se = sd(MF_gamma) / sqrt(n()),
    temp_mean = mean(temp),
    temp_se = sd(temp) / sqrt(n()),
    w_mean = mean(w),
    w_se = sd(w) / sqrt(n())
  ) |> 
  as.data.frame()

hist(dd$MB_Arew)

m1 <- brm(
  MB_Arew ~ is_patient,
  family = student(),
  data = dd,
  chains = 4, 
  cores = 8,
  backend = "cmdstanr"
  #threads = threading(2)
)
pp_check(m1)
summary(m1, prob = 0.89)


hist(dd$MB_Apun)

m2 <- brm(
  MB_Apun ~ is_patient,
  family = student(),
  data = dd,
  chains = 4, 
  cores = 8,
  backend = "cmdstanr"
  #threads = threading(2)
)
pp_check(m2)
summary(m2, prob = 0.89)



# Split into train/test ---------------------------------------------------

d <- dd_wide |> 
  dplyr::select(-subj_idx)

# For a classification model, the outcome should be a `factor`
d$is_patient <- factor(d$is_patient)

set.seed(45823)
# Put 3/4 of the data into the training set.
data_split <- initial_split(d, strata = is_patient, prop = 0.7)
data_split

# The strata argument causes the random sampling to be conducted within the 
# stratification variable. This can help ensure that the number of data points 
# in the training data is equivalent to the proportions in the original data 
# set. 

# Create data frames for the two sets. 
# The training and testing sets can be extracted from the “split” object using 
# the training() and testing() functions.
train_data <- training(data_split)
test_data  <- testing(data_split)

dim(dat)
dim(train_data)
dim(test_data)

# Training set proportions by is_patient.
train_data %>% 
  count(is_patient) %>% 
  mutate(prop = n/sum(n))

# Test set proportions by is_patient.
test_data  %>% 
  count(is_patient) %>% 
  mutate(prop = n/sum(n))

# We’re going to want to do some parameter tuning, and to do that we’re going 
# to want to use cross-validation. So we can create a cross-validated version 
# of the training set.
prl_cv <- vfold_cv(train_data, v = 10, repeats = 500)


# Define a recipe ---------------------------------------------------------

prl_recipe <- 
  recipe(is_patient ~ ., data = d) %>%
  # and some pre-processing steps
  step_zv(all_predictors()) %>% 
  step_normalize(all_numeric()) 
prl_recipe


# Specify the model -------------------------------------------------------

lr_model <- 
  # specify that the model is a logistic regression
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") 


rf_model <- 
  # specify that the model is a random forest
  rand_forest() %>%
  # specify that the `mtry` parameter needs to be tuned
  set_args(mtry = tune()) %>%
  # select the engine/package that underlies the model
  set_engine("ranger", importance = "impurity") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") 




# Create a workflow -------------------------------------------------------

# set the workflow
lr_workflow <- workflow() %>%
  # add the recipe
  add_recipe(prl_recipe) %>%
  # add the model
  add_model(lr_model)


rf_workflow <- workflow() %>%
  # add the recipe
  add_recipe(prl_recipe) %>%
  # add the model
  add_model(rf_model)


# Tune the parameters -----------------------------------------------------

# Note that we still haven’t yet implemented the pre-processing steps in the 
# recipe nor have we fit the model. We’ve just written the framework. It is 
# only when we tune the parameters or fit the model that the recipe and model 
# frameworks are actually implemented.

lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 30))
lr_reg_grid %>% top_n(-5) # lowest penalty values
lr_reg_grid %>% top_n(5)  # highest penalty values

lr_res <- 
  lr_workflow %>% 
  tune_grid(
    resamples = prl_cv, #CV object
    grid = lr_reg_grid,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )

lr_res %>%
  collect_metrics()


# Finalize the workflow ---------------------------------------------------

# We want to add a layer to our workflow that corresponds to the tuned 
# parameter, i.e. sets mtry to be the value that yielded the best results. If 
# you didn’t tune any parameters, you can skip this step.

param_final <- lr_res %>%
  select_best(metric = "roc_auc")
param_final

lr_workflow <- lr_workflow %>%
  finalize_workflow(param_final)


# Evaluate the model on the test set --------------------------------------

lr_fit <- lr_workflow %>%
  # fit on the training set and evaluate on test set
  last_fit(data_split)
lr_fit

test_performance <- lr_fit %>% 
  collect_metrics()
test_performance 
# 1 accuracy binary         0.759 Preprocessor1_Model1
# 2 roc_auc  binary         0.781 Preprocessor1_Model1


# generate predictions from the test set
test_predictions <- lr_fit %>% 
  collect_predictions()
test_predictions

# generate a confusion matrix
test_predictions %>% 
  conf_mat(truth = is_patient, estimate = .pred_class)

test_predictions %>%
  ggplot() +
  geom_density(aes(x = .pred_0, fill = is_patient), 
               alpha = 0.5)

test_predictions <- lr_fit %>% 
  pull(.predictions)
test_predictions


lr_fit %>% 
  pluck(".workflow", 1) %>%   
  extract_fit_parsnip() %>% 
  vip(num_features = 20)

lr_fit %>% 
  collect_predictions() %>% 
  roc_curve(is_patient, .pred_0) %>% 
  autoplot()

#' We employed a logistic regression model to classify patients and controls 
#' based on the five parameters derived from the hDDMrl computational model. 
#' The dataset was split into training (70%) and testing (30%) subsets, 
#' maintaining the original distribution of patient and control groups. 
#' We utilized cross-validation with 10 folds and 500 repeats to optimize 
#' model parameters. Model tuning was conducted using a grid search over a 
#' range of penalty values, and the best-performing model was selected based 
#' on the highest ROC-AUC score. Finally, we evaluated the model's performance 
#' on the test set, achieving an accuracy of 75.9% and an ROC-AUC of 81.0%.

# eof ----

