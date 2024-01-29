
library("here")
library("tidyverse")
library("brms")


prl_dat <- readRDS(
  here::here(
    "data", "processed", "prl", "raw_prl_data", "prl_tot_raw_data.rds")
)

length(unique(prl_dat$subj_name))
# [1] 338

# Recode
prl_dat$feedback <- ifelse(
  prl_dat$feedback == 2, 0, ifelse(prl_dat$feedback == 1, 1, NA)
)

source(here::here("data", "list_participants.R"))

d_an <- prl_dat[prl_dat$subj_name %in% participants_list[[1]], ]
d_an$group <- "AN"
length(unique(d_an$subj_name))

d_hc <- prl_dat[prl_dat$subj_name %in% participants_list[[2]], ]
d_hc$group <- "HC"
length(unique(d_hc$subj_name))

d_ri <- prl_dat[prl_dat$subj_name %in% participants_list[[3]], ]
d_ri$group <- "RI"
length(unique(d_ri$subj_name))

d <- rbind(d_an, d_hc, d_ri)
d$group <- factor(d$group)

d |> 
  group_by(group) |> 
  summarize(
    n = n_distinct(subj_name)
  )

d1 <- d |>
  dplyr::filter(stimulus_type == "food") |> 
  dplyr::select(subj_name, trial, group, feedback, is_target_img_chosen) |> 
  arrange(subj_name, trial)

d1$choice <- d1$is_target_img_chosen
# 1: food
# 0: no food

# Create a shifted column for comparison
d1 <- d1 %>%
  group_by(group, subj_name) %>%
  mutate(
    previous_choice = lag(choice),
    previous_feedback = lag(feedback),
    stay = ifelse(choice == previous_choice, 1, 0),
    perseverative_error = ifelse(
      previous_feedback == 0 & choice == previous_choice, 1, 0)
  ) |> 
  ungroup()

set.seed(1)
imp <- mice(d1)
d1_imp <- mice::complete(imp, 1)

# Compute Win-Stay and Lose-Shift
win_stay_lose_shift <- d1_imp |>
  group_by(group, subj_name) |>
  summarize(
    win_stay = mean(stay[previous_feedback == 1], na.rm = TRUE),
    lose_shift = mean(1 - stay[previous_feedback == 0], na.rm = TRUE),
    e_pers = mean(perseverative_error)
  ) |> 
  ungroup()


hist(win_stay_lose_shift$win_stay)
hist(win_stay_lose_shift$lose_shift)
hist(win_stay_lose_shift$e_pers)


mod_ws <- brm(
  win_stay ~ group,
  family = zero_one_inflated_beta(),
  data = win_stay_lose_shift,
  backend = "cmdstanr"
)
pp_check(mod_ws)
loo_ws <- loo(mod_ws)
plot(loo_ws)

summary(mod_ws, prob = 0.89)
#           Estimate Est.Error l-89% CI u-89% CI Rhat Bulk_ESS Tail_ESS
# Intercept     1.11      0.12     0.92     1.31 1.00     4569     2517
# groupHC       0.31      0.17     0.03     0.58 1.00     4349     2883
# groupRI      -0.02      0.20    -0.33     0.31 1.00     4347     3003

win_stay_lose_shift |> 
  group_by(group) |> 
  summarize(
    ws = median(win_stay),
    ls = median(lose_shift)
  )

bayes_R2(mod_ws, probs = 0.89)
#      Estimate  Est.Error        Q89
# R2 0.04367431 0.03043061 0.08351171


mod_ls <- brm(
  lose_shift ~ group,
  family = asym_laplace(),
  data = win_stay_lose_shift,
  backend = "cmdstanr"
)
pp_check(mod_ls)
loo_ls <- loo(mod_ls)
plot(loo_ls)

summary(mod_ls, prob = 0.89)
#           Estimate Est.Error l-89% CI u-89% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.54      0.02     0.50     0.58 1.00     2187     2277
# groupHC       0.04      0.03    -0.00     0.08 1.00     2459     2265
# groupRI       0.04      0.03    -0.01     0.09 1.00     2525     2489

bayes_R2(mod_ls, probs = 0.89)
#      Estimate  Est.Error        Q89
# R2 0.02453708 0.01885976 0.04930178



mod_epers <- brm(
  e_pers ~ group,
  family = asym_laplace(),
  data = win_stay_lose_shift,
  backend = "cmdstanr"
)
pp_check(mod_epers)
loo_ep <- loo(mod_epers)
plot(loo_ep)

summary(mod_epers, prob = 0.89)
#           Estimate Est.Error l-89% CI u-89% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.23      0.01     0.21     0.25 1.00     2069     2496
# groupHC      -0.04      0.01    -0.06    -0.02 1.00     1793     2618
# groupRI      -0.03      0.01    -0.05    -0.00 1.00     2004     2276

win_stay_lose_shift |> 
  group_by(group) |> 
  summarize(
    p = mean(e_pers)
  )
# group     p
# 1 AN    0.220
# 2 HC    0.178
# 3 RI    0.191
