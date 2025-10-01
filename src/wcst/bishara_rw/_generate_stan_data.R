# Overview ----------------------------------------------------------------
# Associated project: WCST, taks-switching, PRL
# Script purpose: save stan_data for further processing (models' comparison
# for the WCST task)
#
# Written by: Corrado Caudek (corrado.caudek@unifi.it)
# Version: 2025-10-01
# Last update:
# Status: final
# Notes:

# Load necessary libraries ------------------------------------------------

if (!requireNamespace("pacman")) install.packages("pacman")

pacman::p_load()

# Load necessary libraries ------------------------------------------------

if (!requireNamespace("pacman")) install.packages("pacman")

pacman::p_load(
  here,
  tidyverse,
  cmdstanr,
  posterior,
  bayesplot,
  insight,
  pROC,
  stringr,
  loo
)


# Source functions --------------------------------------------------------

source(here::here(
  "src",
  "wcst",
  "documentation",
  "functions",
  "funs_wcst.R"
))

source(here::here(
  "src",
  "wcst",
  "documentation",
  "functions",
  "funs_model_selection_wcst.R"
))

# process_and_prepare_stan_data()
source(
  here::here(
    "src",
    "wcst",
    "documentation",
    "functions",
    "funs_input_for_stan_wcst.R"
  )
)

# Data preparation ----------------------------------------------------------

generate_csv("controls")
generate_csv("patients")

# Generate the input for the Stan models
stan_data <- process_and_prepare_stan_data()
str(stan_data)

# Save data -----------------------------------------------------------------

saveRDS(
  stan_data,
  here::here(
    "src",
    "wcst",
    "data",
    "wcst_stan_list.RDS"
  )
)

# eof -----
