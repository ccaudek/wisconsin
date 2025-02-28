# Script name: 001_wcst_descriptive_stats.R
# Project: WCST project
# Script purpose: Compare descriptve stats between patients and controls
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Wed Jan 24 14:35:42 2024
# Last Modified Date: Wed Jan 24 14:35:42 2024
#
# 👉 

# Prelims
suppressPackageStartupMessages({
  library("here")
  library("tidyverse")
  library("brms")
  library("cmdstanr")
})

# Version 2.34.1 does not work
set_cmdstan_path("/Users/corrado/.cmdstan/cmdstan-2.33.1")
cmdstan_path()
cmdstan_version()


d <- readRDS(
  here("data", "processed", "wcst", "behavioral_stats.rds")
)

d %>% 
  group_by(group) %>% 
  summarise(
    prop_pers_err = mean(prop_pers_err),
    prop_non_pers_err = mean(prop_non_pers_err),
    prop_err = mean(prop_err),
    prop_cor = mean(prop_cor),
    n = n_distinct(subj_name)
  )

hist(d$prop_pers_err)

bf1 <- bf(prop_pers_err ~ group)

fit1 <- brm(
  bf1, 
  data = d, 
  family = asym_laplace(),
  chains = 4, cores = 4
)

pp_check(fit1)
loo1 <- loo(fit1)
plot(loo1)

summary(fit1)
conditional_effects(fit1, "group")



# prop_non_pers_err ----

hist(d$prop_non_pers_err)

d |> 
  group_by(group) |> 
  summarize(
    e = mean(prop_non_pers_err)
  )

bf2 <- bf(prop_non_pers_err ~ group)

fit2 <- brm(
  bf2, 
  data = d, 
  family = asym_laplace(),
  chains = 4, cores = 4
)

pp_check(fit2)

summary(fit2)
conditional_effects(fit2, "group")


# prop_err -----

hist(d$prop_err)

bf3 <- bf(prop_err ~ group)

fit3 <- brm(
  bf3, 
  data = d, 
  family = asym_laplace(),
  chains = 4, cores = 4
)

pp_check(fit3)

summary(fit3)
conditional_effects(fit3, "group")


# prop_cor ---


hist(d$prop_cor)

bf4 <- bf(prop_cor ~ group)

fit4 <- brm(
  bf4, 
  data = d, 
  family = asym_laplace(),
  chains = 4, cores = 4
)

pp_check(fit4)

summary(fit4)
conditional_effects(fit4, "group")


# Perseverative responses ------------------------------------------------------

# Import raw data for all three groups included in the project
raw_df <- rio::import(
  here::here(
    "data", "processed", "wcst", "raw_data_project.csv"
  )
)
length(unique(raw_df$subj_name))

raw_df$group <- factor(raw_df$group)
raw_df$group <- relevel(raw_df$group, ref = "hc")

raw_df$is_correct <- ifelse(
  raw_df$is_correct == 2, 0, 
  ifelse(raw_df$is_correct == 1, 1, NA)
)



# Calculate perseverative responses 
# perseverative_df <- raw_df %>%
#   group_by(subj_name) %>%
#   mutate(
#     prev_chosen_card = lag(chosen_card),
#     prev_correct_card = lag(correct_card),
#     # A perseverative response occurs when the chosen card matches the previous chosen card
#     # and the current correct card is different from the previous correct card
#     is_pers_resp = if_else(
#       chosen_card == prev_chosen_card & correct_card != prev_correct_card, 1, 0),
#     is_pers_err = if_else(
#       (chosen_card == prev_chosen_card) & (correct_card != prev_correct_card) & (is_correct == 2), 1, 0)
#   ) %>%
#   ungroup()

# perseverative_df <- raw_df %>%
#   group_by(subj_name) %>%
#   mutate(
#     prev_name_of_task = lag(name_of_task),
#     rule_changed = if_else(name_of_task != prev_name_of_task & !is.na(prev_name_of_task), TRUE, FALSE),
#     post_rule_change_period = cumsum(rule_changed)
#   ) %>%
#   mutate(
#     prev_chosen_card = lag(chosen_card),
#     is_pers_resp = if_else(
#       chosen_card == prev_chosen_card & post_rule_change_period == lag(post_rule_change_period) & !is.na(prev_name_of_task), 1, 0, missing = 0),
#     is_pers_err = if_else(
#       is_pers_resp == 1 & is_correct == 0, 1, 0, missing = 0)
#   ) %>%
#   ungroup()

perseverative_df <- raw_df %>%
  group_by(subj_name) %>%
  mutate(
    prev_name_of_task = lag(name_of_task),
    rule_changed = if_else(name_of_task != prev_name_of_task & !is.na(prev_name_of_task), TRUE, FALSE),
    post_rule_change_period = cumsum(rule_changed),
    trial_num = row_number()  # Add a column for trial number within each subject
  ) %>%
  mutate(
    prev_chosen_card = lag(chosen_card),
    is_pers_resp = if_else(
      chosen_card == prev_chosen_card & post_rule_change_period == lag(post_rule_change_period) & 
        !is.na(prev_chosen_card), 1, 0, missing = 0),
    is_pers_err = if_else(
      is_pers_resp == 1 & is_correct == 0 & trial_num > 10, 1, 0, missing = 0)
  ) %>%
  ungroup()

# Calculate proportion of perseverative responses for each participant
proportion_perseverative <- perseverative_df %>%
  group_by(group, subj_name) %>%
  summarize(
    total_responses = n(),
    pers_resp = sum(is_pers_resp, na.rm = TRUE),
    pers_err  = sum(is_pers_err, na.rm = TRUE),
    prop_pers_resp = pers_resp / total_responses,
    prop_pers_err = pers_err / 50,
  ) %>%
  ungroup()

# Display the result
proportion_perseverative %>%
  group_by(group) %>%
  summarise(
    prop_pers_resp = mean(prop_pers_resp),
    prop_pers_err = mean(prop_pers_err)
  )

# Model perseverative responses
mod1 <- brm(
  prop_pers_resp ~ group, 
  data = proportion_perseverative, 
  family = asym_laplace(),
  chains = 4, 
  cores = 4, 
  backend = "cmdstanr"
)
pp_check(mod1)
summary(mod1)
conditional_effects(mod1, "group")

# Model perseverative errors
mod2 <- brm(
  prop_pers_err ~ group, 
  data = proportion_perseverative, 
  family = asym_laplace(),
  chains = 4, 
  cores = 4, 
  iter = 10000,
  backend = "cmdstanr"
)
pp_check(mod2)
summary(mod2)
conditional_effects(mod2, "group")


# Model errors
mod3 <- brm(
  is_error ~ group + (1 | subj_name), 
  data = raw_df, 
  family = bernoulli(),
  chains = 4, 
  cores = 4, 
  backend = "cmdstanr"
)
pp_check(mod3)
summary(mod3)
conditional_effects(mod3, "group")


