# Script name: 003_descript_stats.R
# Project: Eating disorders Montecatini
# Script purpose: Descriptive stats on task-switching RT data
# @author: Corrado Caudek <corrado.caudek@unifi.it>
# Date Created: Tue Jun  7 09:20:09 2022
# Last Modified Date: Tue Jun  7 09:20:09 2022
#
# 👉 

suppressPackageStartupMessages({
  library("here")
  library("tidyverse")
  library("brms")
  library("tidybayes")        # Manipulate Stan objects in a tidy way
  library("broom")            # Convert model objects to data frames
  library("broom.mixed")      # Convert brms model objects to data frames
  library("ggdist")           # Special geoms for posterior distributions
  library("ggrepel")          # Automatically position labels
  library("patchwork")        # Combine ggplot objects
  library("emmeans")          # Calculate marginal effects in even fancier ways
  library("cmdstanr")
  library("tidyr")
})

# Version 2.34.1 does not work
set_cmdstan_path("/Users/corrado/.cmdstan/cmdstan-2.33.1")
cmdstan_path()
cmdstan_version()

# Source functions 
source(here::here("src", "functions", "funs_task_switching.R"))


# Read data ---------------------------------------------------------------

# Read raw switching data
d <- read_clean_task_switch_data()

# Add response variable
d1 <- add_response_var(d)

# Wrangle data
wrangle_task_switch_data(d1)

task_switch_df <- readRDS(
  here::here(
    "data", "processed", "task_switching", 
    "task_switch_data_for_descript_stats.rds"
  )
)

task_switch_df %>% 
  group_by(group) %>% 
  summarise(
    n = n_distinct(subj_id)
  )

# fm <- lme4::lmer(log(rt/1000) ~ diag_cat + (1 | subj_code), two_groups)
# car::Anova(fm)


# RT plot 

plot_df <- task_switch_df %>%
  dplyr::filter(is_correct == 1) |> 
  group_by(resp_transition, task_transition, group) %>%
  summarise(
    sd = sd(rt) / sqrt(n()),
    y = median(rt, na.rm = TRUE)
  ) |> 
  ungroup()

plot_df <- plot_df %>%
  mutate_at(
    c('resp_transition', 'task_transition', 'group'), 
    as.factor
  )

plot_df %>% 
  ggplot(aes(resp_transition, y)) +
  geom_line(
    aes(#linetype = task_transition, 
        group = task_transition,
        color = task_transition)
    ) +
  geom_point() +
  geom_pointrange(
    aes(ymin = y-sd, ymax = y+sd, color = task_transition)
  ) +
  facet_grid(~ group) +
  labs(
    x = "Response transition",
    y = "Reaction times (ms)",
    color = "Task transition"
  ) +
  theme(legend.position = "bottom")

hist(task_switch_df$rt)


# Error rates

plot_df <- task_switch_df %>%
  group_by(group, subj_id, resp_transition, task_transition) %>%
  summarise(
    er = 100 * (1 - mean(is_correct, na.rm = TRUE)),
  ) %>% 
  ungroup() %>% 
  group_by(resp_transition, task_transition, group) %>%
  summarise(
    sd = sd(er) / sqrt(n()),
    y = mean(er, trim = 0.1, na.rm = TRUE)
  )

plot_df <- plot_df %>%
  mutate_at(
    c('resp_transition', 'task_transition', 'group'), 
    as.factor
  )


plot_df %>% 
  ggplot(aes(resp_transition, y)) +
  geom_line(
    aes(#linetype = task_transition, 
      group = task_transition,
      color = task_transition)
  ) +
  geom_point() +
  geom_pointrange(
    aes(ymin = y-sd, ymax = y+sd, color = task_transition)
  ) +
  facet_grid(~ group) +
  labs(
    x = "Response transition",
    y = "Error rate (%)",
    color = "Task transition"
  ) +
  theme(legend.position = "bottom")

# Model for RTs

task_switch_df$RT <- task_switch_df$rt / 1000

correct_data_df <- task_switch_df |> 
  dplyr::filter(is_correct == 1 & RT < 3.5)

hist(correct_data_df$RT)

# Define a function to identify outliers based on the Tukey criterion
identify_outliers <- function(x) {
  q1 <- quantile(x, .25)
  q3 <- quantile(x, .75)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  return(x >= lower_bound & x <= upper_bound)
}

# Apply the function to each subj_id group and filter out outliers
filtered_data_df <- correct_data_df %>%
  group_by(subj_id) %>%
  mutate(outlier = identify_outliers(rt)) %>%
  filter(outlier) %>%
  select(-outlier)

hist(filtered_data_df$RT)


fit_rt <- brm(
  RT ~ block_type + (task_transition * resp_transition) * group +
    (block_type + task_transition * resp_transition | subj_id) +
    (1 | food_img) + (1 | plant_img), 
  data = filtered_data_df, 
  # algorithm = "meanfield",
  family = shifted_lognormal(),
  chains = 4, 
  cores = 4,
  #iter = 5000,
  # control = list(
  #   adapt_delta = 0.95,
  #   max_treedepth = 15
  # ),
  threads = threading(4),
  backend = "cmdstanr"
)
pp_check(fit_rt)

fit1_rt <- brm(
  RT ~ block_type + (task_transition * resp_transition) + group +
    (block_type + task_transition * resp_transition | subj_id) +
    (1 | food_img) + (1 | plant_img), 
  data = filtered_data_df, 
  # algorithm = "meanfield",
  family = shifted_lognormal(),
  chains = 4, 
  cores = 4,
  #iter = 5000,
  # control = list(
  #   adapt_delta = 0.95,
  #   max_treedepth = 15
  # ),
  threads = threading(4),
  backend = "cmdstanr"
)


conditional_effects(fit_rt, "task_transition:group")
summary(fit_rt)

loo1 <- loo(fit_rt)
plot(loo1)

loo2 <- loo(fit2_rt)
loo_compare(loo1, loo2)


# eof ------------







error_rates_df <- two_groups %>%
  group_by(diag_cat, subj_code, resp_transition, task_transition) %>%
  summarise(
    er = (1 - mean(is_correct, na.rm = TRUE)),
  ) %>% 
  ungroup()







hist(error_rates_df$er)


error_rates_df$error_rate <- error_rates_df$er / 100

error_rates_df$error_rate <- ifelse(
  error_rates_df$error_rate == 0, 1/120, error_rates_df$error_rate
)


only_hc <- error_rates_df[error_rates_df$diag_cat == "AN", ]


bf_er <- bf(
  error_rate ~ resp_transition * task_transition + 
    (resp_transition + task_transition | subj_code)
)

fit_er <- brm(
  bf_er, 
  data = only_hc, 
  # family = exponential(),
  family = Beta(),
  chains = 4, 
  cores = 4,
  iter = 5000,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 15
  ),
  threads = threading(4),
  backend = "cmdstan"
)

pp_check(fit_er)

summary(fit_er)
broom.mixed::tidy(fit_er, effects = "fixed")

conditions <- make_conditions(mod_rt, "diag_cat")
conditional_effects(
  mod_rt, 
  "resp_transition:task_transition",
  conditions = conditions)

bayes_R2(fit_er)



# Classification ----------------------------------------------------------


two_groups <- two_groups %>%  
  mutate_if(is.character, as.factor)

two_groups$task_transition <- dplyr::recode_factor(
  two_groups$task_transition,
  "repeat" = "repetition" 
)
summary(two_groups$task_transition)


all_df <- two_groups %>%
  group_by(diag_cat, subj_code, resp_transition, task_transition) %>%
  summarise(
    er = 100 * (1 - mean(is_correct, na.rm = TRUE)),
  ) %>% 
  ungroup() 
summary(all_df)


hc_df <- all_df %>% 
  dplyr::filter(diag_cat == "HC")
hc_df$is_patient <- 0

an_df <- all_df %>% 
  dplyr::filter(diag_cat == "AN")
an_df$is_patient <- 1






# select random 29 rows of the dataframe 
# set.seed(2)
subj_names_hc <- unique(hc_df$subj_code)
chosen_subjects <- sample(subj_names_hc, 29)

hc29_df <- hc_df[hc_df$subj_code %in% chosen_subjects, ]

two_groups29_df <- rbind(an_df, hc29_df)
summary(two_groups29_df)


two_groups29_wider_df <- two_groups29_df %>%
  pivot_wider(
    #!c("diag_cat", "subj_code", "resp_transition", "is_patient"),
    names_from = task_transition, values_from = er)
summary(two_groups29_wider_df)

two_groups29_wider_df <- two_groups29_wider_df %>% 
  group_by(subj_code, resp_transition) %>% 
  mutate(
    tt_effect = switch - repetition
  )
summary(two_groups29_wider_df)

two_groups29_wider_df$repetition <- NULL
two_groups29_wider_df$switch <- NULL
two_groups29_wider_df$diag_cat <- NULL

two_groups29_wider2_df <- two_groups29_wider_df %>%
  pivot_wider(names_from = resp_transition, values_from = tt_effect)

two_groups29_wider2_df <- two_groups29_wider2_df %>% 
  group_by(subj_code) %>% 
  mutate(
    task_switch_eff = switch - repetition
  )

two_groups29_wider2_df %>% 
  group_by(is_patient) %>% 
  summarise(
    m = mean(task_switch_eff, na.rm = TRUE), 
    sd = sd(task_switch_eff, na.rm = TRUE) / sqrt(n())
  )


fm <- glm(
  is_patient ~ task_switch_eff,
  family = binomial(),
  data = two_groups29_wider2_df
)

car::Anova(fm)


library("pROC")
test_prob = predict(fm, newdata = two_groups29_wider2_df, type = "response")
test_roc = roc(two_groups29_wider2_df$is_patient ~ test_prob, plot = TRUE, print.auc = TRUE)


# End of file -------------------------------------------------------------

nrep <- 1000
res_auc <- rep(NA, nrep)
for (i in 1:nrep) {
  
  # select random 29 rows of the dataframe 
  subj_names_hc <- unique(hc_df$subj_code)
  chosen_subjects <- sample(subj_names_hc, 29)
  
  hc29_df <- hc_df[hc_df$subj_code %in% chosen_subjects, ]
  
  two_groups29_df <- rbind(an_df, hc29_df)
  
  temp <- two_groups29_df %>%
    pivot_wider(names_from = c(resp_transition, task_transition), values_from = er)
  
  temp$diag_cat <- NULL
  temp$subj_code <- NULL
  
  fm <- glm(
    is_patient ~ .,
    family = binomial(),
    data = temp
  )
  
  test_prob = predict(fm, newdata = temp, type = "response")
  test_roc = roc(temp$is_patient ~ test_prob, plot = FALSE, print.auc = FALSE)
  
  res_auc[i] <- readr::parse_number(as.character(test_roc$auc))

}

hist(res_auc)
mean(res_auc)
sd(res_auc)

#-------------------------------------------------------------------------------








