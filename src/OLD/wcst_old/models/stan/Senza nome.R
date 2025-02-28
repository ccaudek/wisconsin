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
    "src", "wcst", "traces", "wcst_params_classification.csv"
  )
)

dat$is_patient <- ifelse(params_2grps$index < 46, 1, 0)
dat$is_patient <- factor(dat$is_patient)

dat$index <- NULL

dd <- dat
names(dd)

model_glm = glm(is_patient ~ ., data = dd, family = "binomial")

test_prob = predict(model_glm, newdata = dd, type = "response")
test_roc = roc(dd$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)


# Split into train/test ---------------------------------------------------

set.seed(33625)
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
wcst_cv <- vfold_cv(train_data, v = 10, repeats = 200)


# Define a recipe ---------------------------------------------------------

wcst_recipe <- 
  recipe(is_patient ~ ., data = dd) %>%
  # and some pre-processing steps
  step_normalize(all_numeric()) %>%
  step_impute_knn(all_predictors())
prl_recipe


# Specify the model ------------------------------------------------------------

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

rf_workflow <- workflow() %>%
  # add the recipe
  add_recipe(wcst_recipe) %>%
  # add the model
  add_model(rf_model)


# Tune the parameters -----------------------------------------------------

# Note that we still haven’t yet implemented the pre-processing steps in the 
# recipe nor have we fit the model. We’ve just written the framework. It is 
# only when we tune the parameters or fit the model that the recipe and model 
# frameworks are actually implemented.

# specify which values eant to try
rf_grid <- expand.grid(mtry = c(3, 4, 5))
# extract results
rf_tune_results <- rf_workflow %>%
  tune_grid(resamples = wcst_cv, #CV object
            grid = rf_grid, # grid of values to try
            metrics = metric_set(accuracy, roc_auc) # metrics we care about
  )

# print results
rf_tune_results %>%
  collect_metrics()


# Finalize the workflow ---------------------------------------------------

# We want to add a layer to our workflow that corresponds to the tuned 
# parameter, i.e. sets mtry to be the value that yielded the best results. If 
# you didn’t tune any parameters, you can skip this step.

param_final <- rf_tune_results %>%
  select_best(metric = "accuracy")
param_final

rf_workflow <- rf_workflow %>%
  finalize_workflow(param_final)


# Evaluate the model on the test set --------------------------------------

rf_fit <- rf_workflow %>%
  # fit on the training set and evaluate on test set
  last_fit(data_split)
rf_fit


test_performance <- rf_fit %>% collect_metrics()
test_performance
# 1 accuracy binary         0.607 Preprocessor1_Model1
# 2 roc_auc  binary         0.673 Preprocessor1_Model1


# generate predictions from the test set
test_predictions <- rf_fit %>% collect_predictions()
test_predictions

# generate a confusion matrix
test_predictions %>% 
  conf_mat(truth = is_patient, estimate = .pred_class)

test_predictions %>%
  ggplot() +
  geom_density(aes(x = .pred_1, fill = is_patient), 
               alpha = 0.5)

test_predictions <- rf_fit %>% pull(.predictions)
test_predictions

final_model <- fit(rf_workflow, dd)
final_model

ranger_obj <- pull_workflow_fit(final_model)$fit
ranger_obj$variable.importance



# eof ----

