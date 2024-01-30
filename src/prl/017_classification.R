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
  # Helper packages
  library("readr")       # for importing data
  library("vip") 
})


# Increase max print
options(max.print = .Machine$integer.max)

cores <- parallel::detectCores()
cores

source(here("code", "functions", "funs_prl.R"))
source(here("code", "functions", "funs_quest.R"))
source(here("code", "functions", "funs_rescorla_wagner.R"))
source(here("code", "functions", "funs_careless_resp.R"))
source(here("code", "functions", "funs_gen_data_for_hddm.R"))
source(here("code", "functions", "funs_param_analyses.R"))
source(here("code", "functions", "funs_param_classification_analyses.R"))


# Get parameter estimates


# if file does not exist:
if(
  !file.exists(
    here("data", "processed", "prl", "prl_params_cleaned", 
         "prl_params_no_split_by_group_and_subj_code.rds"
         # "prl_params_and_subj_code.rds"
    )
  )
) add_subj_code_to_hddm_params_and_save()
# add_subj_code_to_hddm_params_and_save()






params <- rio::import(
  here::here("data", "processed", "prl", "params", "subj_code_hddm_params.csv")
)

params_2grps <- params |> 
  dplyr::filter(group != "atrisk")

params_2grps$is_patient <- ifelse(
  params_2grps$group == "patients", 1, 0
)

# Get patients codes.
# patients_codes <- find_patients_codes()
# params$is_patient <- ifelse(params$subj_code %in% patients_codes, 1, 0)

#only_neutral <- params %>% 
#  dplyr::filter(stim == "neutral")
#
#only_food <- params %>% 
#  dplyr::filter(stim == "food")

# Select Dataframe.
dat <- params_2grps

dat$params <- NULL
dat$subj_code <- NULL
dat$diag_cat <- NULL
dat$group <- NULL
dat$is_patient <- factor(dat$is_patient)

dat_wide <- dat %>%
  pivot_wider(names_from = par, values_from = val)

glimpse(dat_wide)


library(pROC)

dd <- dat_wide |> 
  dplyr::select(-subj_idx)

model_glm = glm(is_patient ~ ., data = dd, family = "binomial")

test_prob = predict(model_glm, newdata = dd, type = "response")
test_roc = roc(dd$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)


# Split into train/test ---------------------------------------------------

set.seed(54645)
# Put 3/4 of the data into the training set.
data_split <- initial_split(dd, strata = is_patient, prop = 0.7)
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

dim(dat_wide)
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
prl_cv <- vfold_cv(train_data, v = 10, repeats = 100)


# Define a recipe ---------------------------------------------------------

prl_recipe <- 
  recipe(is_patient ~ ., data = dd) %>%
  # and some pre-processing steps
  step_zv(all_predictors()) %>% 
  step_normalize(all_numeric()) 
  # step_impute_knn(all_predictors()) %>% 
  # step_dummy(all_nominal_predictors()) 
prl_recipe

# prl_train_preprocessed <- prl_recipe %>%
#   # apply the recipe to the training data
#   prep(train_data) %>%
#   # extract the pre-processed training dataset
#   juice()
# prl_train_preprocessed


# Specify the model -------------------------------------------------------

lr_model <- 
  # specify that the model is a logistic regression
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") 

# rf_model <- 
#   # specify that the model is a random forest
#   rand_forest() %>%
#   # specify that the `mtry` parameter needs to be tuned
#   set_args(mtry = tune()) %>%
#   # select the engine/package that underlies the model
#   set_engine("ranger", importance = "impurity") %>%
#   # choose either the continuous regression or binary classification mode
#   set_mode("classification") 


# Create a workflow -------------------------------------------------------

# set the workflow
lr_workflow <- workflow() %>%
  # add the recipe
  add_recipe(prl_recipe) %>%
  # add the model
  add_model(lr_model)


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
# 1 accuracy binary         0.793 Preprocessor1_Model1
# 2 roc_auc  binary         0.871 Preprocessor1_Model1


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

