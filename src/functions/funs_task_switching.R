

# logistic_reg_classification() -------------------------------------------

logistic_reg_classification <- function() {
  
  task_switch_df <- readRDS(
    here::here(
      "data", "processed", "task_switching", "data_for_descript_stats", 
      "task_switch_data_for_descript_stats_clean.rds"
    )
  )
  
  two_groups <- task_switch_df %>% 
    dplyr::filter(diag_cat == "AN" | diag_cat == "HC")
  
  two_groups$diag_cat <- factor(two_groups$diag_cat)
  two_groups$is_error <- ifelse(two_groups$is_correct == 1, 0, 1)
  
  
  out <- two_groups %>% 
    group_by(diag_cat) %>% 
    summarise(
      n = n_distinct(subj_code)
    ) %>% 
    ungroup()
  # 31 patients completed the task-switching experiment.
  N_PATIENTS <- out[2, 2] %>% as.numeric()
  
  hc_df <- two_groups %>% 
    dplyr::filter(diag_cat == "HC")
  hc_df$is_patient <- 0
  
  an_df <- two_groups %>% 
    dplyr::filter(diag_cat == "AN")
  an_df$is_patient <- 1
  
  nrep <- 1000 # repetitions of the AUC classification calculation, after
               # randomly selecting 31 HC controls from the 254 available.
  
  res_auc <- rep(NA, nrep)
  
  for (i in 1:nrep) {
    
    # Select random 31 HCs. 
    subj_names_hc <- unique(hc_df$subj_code)
    chosen_subjects <- sample(subj_names_hc, N_PATIENTS)
    
    hc31_df <- hc_df[hc_df$subj_code %in% chosen_subjects, ]
    
    two_groups31_df <- rbind(an_df, hc31_df)
    
    temp <- two_groups31_df %>%
      pivot_wider(
        names_from = c(resp_transition, task_transition), 
        values_from = is_error
      ) %>% 
      dplyr::select(-c(diag_cat, subj_code)) %>% 
      ungroup()
    
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
  
}



# save_fitted_rt_model_task_switch() --------------------------------------

save_fitted_rt_model_task_switch <- function() {
  
  task_switch_df <- readRDS(
    here::here(
      "data", "processed", "task_switching", "data_for_descript_stats", 
      "task_switch_data_for_descript_stats_clean.rds"
    )
  )
  
  two_groups <- task_switch_df %>% 
    dplyr::filter(diag_cat == "AN" | diag_cat == "HC")
  
  two_groups$diag_cat <- factor(two_groups$diag_cat)
  
  df <- two_groups %>% 
    dplyr::select(
      rt, diag_cat, resp_transition, task_transition, subj_code
    )
  
  bysubj_df <- df %>% 
    group_by(diag_cat, subj_code, resp_transition, task_transition) %>% 
    summarise(
      mrt = median(rt)
    ) %>% 
    mutate(rt = mrt / 1000) %>% 
    dplyr::select(-mrt)

  bf_rt <- bf(
    rt ~ diag_cat * resp_transition * task_transition + 
      (resp_transition + task_transition | subj_code)
  )
  
  mod_rt <- brm(
    bf_rt, 
    data = bysubj_df, 
    family = lognormal(),
    chains = 4, 
    cores = 4,
    iter = 3000,
    file = here::here("scripts", "R", "scripts_task_switching", "brm_files",
                      "mod_rt"),
    backend = "cmdstan"
  )
  
  # pp_check(mod_rt)
  # 
  # summary(mod_rt)
  # broom.mixed::tidy(mod_rt, effects = "fixed")
  # 
  # conditions <- make_conditions(mod_rt, "diag_cat")
  # conditional_effects(
  #   mod_rt, 
  #   "resp_transition:task_transition",
  #   conditions = conditions)
  # 
  # bayes_R2(fit_er)
}




# save_fitted_err_rate_model_task_switch() --------------------------------


save_fitted_err_rate_model_task_switch <- function() {
  
  task_switch_df <- readRDS(
    here::here(
      "data", "processed", "task_switching", "data_for_descript_stats", 
      "task_switch_data_for_descript_stats_clean.rds"
    )
  )
  
  two_groups <- task_switch_df %>% 
    dplyr::filter(diag_cat == "AN" | diag_cat == "HC")
  
  two_groups$diag_cat <- factor(two_groups$diag_cat)
  
  two_groups$is_error <- ifelse(two_groups$is_correct == 1, 0, 1)
  
  df <- two_groups %>% 
    dplyr::select(
      is_error, diag_cat, resp_transition, task_transition, subj_code
    )
  
  bf_er <- bf(
    is_error ~ diag_cat * resp_transition * task_transition + 
      (resp_transition + task_transition | subj_code)
  )
  
  mod_error_rate <- brm(
    bf_er, 
    data = df, 
    family = bernoulli(),
    chains = 4, 
    cores = 4,
    iter = 4000,
    # control = list(
    #   adapt_delta = 0.95,
    #   max_treedepth = 15
    # ),
    # threads = threading(4),
    algorithm = "meanfield",
    tol_rel_obj = 0.00000001,
    file = here::here("scripts", "R", "scripts_task_switching", "brm_files", 
                      "mod_error_rate"),
    backend = "cmdstan"
  )
  
}



# gen_task_switch_performance_fig() ---------------------------------------

gen_task_switch_performance_fig <- function() {
  
  fig_rt <- gen_task_switch_rt_fig()
  fig_acc <- gen_task_switch_accuracy_fig()
  
  p <- fig_rt / fig_acc
  
  ggsave(
    here::here("reports", "figures", "task_switch_performance.pdf"), 
    p, 
    width = 6, height = 7)
  
}


gen_task_switch_rt_fig <- function() {
  
  task_switch_df <- readRDS(
    here::here(
      "data", "processed", "task_switching", "data_for_descript_stats", 
      "task_switch_data_for_descript_stats_clean.rds"
    )
  )
  
  two_groups <- task_switch_df %>% 
    dplyr::filter(diag_cat == "AN" | diag_cat == "HC")
  
  two_groups %>% 
    group_by(diag_cat) %>% 
    summarise(
      n = n_distinct(subj_code)
    )
  
  two_groups$diag_cat <- factor(two_groups$diag_cat, levels=c('AN', 'HC'))
  
  # RT plot 
  
  plot_df <- two_groups %>%
    group_by(resp_transition, task_transition, diag_cat) %>%
    summarise(
      sd = sd(rt) / sqrt(n()),
      y = mean(rt, trim = 0.1, na.rm = TRUE)
    )
  
  plot_df <- plot_df %>%
    mutate_at(
      c('resp_transition', 'task_transition', 'diag_cat'), 
      as.factor
    )
  
  p <- plot_df %>% 
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
    facet_grid(~ diag_cat) +
    labs(
      x = "Response transition",
      y = "Mean reaction time (ms)",
      color = "Task transition"
    ) +
    theme_apa() +
    scale_colour_grey(start = 0.4, end = 0.8) +
    theme(legend.position = "none") 
    
  p
}


# gen_task_switch_accuracy_fig() ------------------------------------------

gen_task_switch_accuracy_fig <- function() {
  
  task_switch_df <- readRDS(
    here::here(
      "data", "processed", "task_switching", "data_for_descript_stats", 
      "task_switch_data_for_descript_stats_clean.rds"
    )
  )
  
  two_groups <- task_switch_df %>% 
    dplyr::filter(diag_cat == "AN" | diag_cat == "HC")
  
  two_groups %>% 
    group_by(diag_cat) %>% 
    summarise(
      n = n_distinct(subj_code)
    )
  
  two_groups$diag_cat <- factor(two_groups$diag_cat, levels=c('AN', 'HC'))
  
  plot_df <- two_groups %>%
    group_by(diag_cat, subj_code, resp_transition, task_transition) %>%
    summarise(
      er = 100 * (1 - mean(is_correct, na.rm = TRUE)),
    ) %>% 
    ungroup() %>% 
    group_by(resp_transition, task_transition, diag_cat) %>%
    summarise(
      sd = sd(er) / sqrt(n()),
      y = mean(er, trim = 0.1, na.rm = TRUE)
    )
  
  plot_df <- plot_df %>%
    mutate_at(
      c('resp_transition', 'task_transition', 'diag_cat'), 
      as.factor
    )
  
  p <- plot_df %>% 
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
    facet_grid(~ diag_cat) +
    labs(
      x = "Response transition",
      y = "Error rate (%)",
      color = "Task transition"
    ) +
    theme_apa() +
    scale_colour_grey(start = 0.4, end = 0.8) +
    theme(legend.position = "bottom") 
  
  p
}


# read_and_tidy_switch_data() ---------------------------------------------

read_and_tidy_switch_data <- function() {
  
  dat_1 <- read_clean_task_switch_data()
  dat_2 <- add_response_var(dat_1)
  dat_3 <- wrangle_task_switch_data(dat_2)
  
}


# read_clean_task_switch_data() -------------------------------------------

read_clean_task_switch_data <- function() {
  
  # Read raw data
  task_switching <- rio::import(
    here::here(
      "data", "processed", "task_switching", "task_switching_data.csv"
    )
  )
  length(unique(task_switching$subj_id))
  
  # Remove subjects with accuracy smaller than 0.85. There are two patients
  # "ch_na_2007_06_23_908_f" "ch_pi_2004_02_25_126_f", with accuracty equal to
  # 0.842 and 0.8. I keep them.
  bad_subjects <- c(
    "an_to_2000_11_30_575_f", "ba_pa_2000_01_08_543_f", "bi_sa_2001_03_01_675_f",
    # "ch_na_2007_06_23_908_f", "ch_pi_2004_02_25_126_f", 
    "cr_ci_1999_08_21_931_f",
    "da_sc_1993_01_08_813_f", "do_za_2002_03_14_283_f", "el_la_1999_06_15_464_f",
    "em_sa_2001_09_21_707_f", "fr_as_1997_02_10_127_f", "fr_ba_1997_10_29_663_f",
    "fr_da_1994_04_19_591_f", "gi_me_2001_03_31_627_f", "gr_po_2002_02_12_110_f",
    "pa_pe_1994_12_31_482_f", "pa_pi_2001_09_28_281_f", "va_ve_1991_11_08_609_f",
    "vi_bi_1996_01_17_478_f", "vi_te_2001_06_08_644_f"
  )
  
  task_switching_clean <- task_switching[
    !(task_switching$subj_id %in% bad_subjects), 
  ]
  
  task_switching_clean$is_correct <- ifelse(
    task_switching_clean$is_correct == 2, 0, 
    ifelse(
      task_switching_clean$is_correct == 3, NA, 
      task_switching_clean$is_correct
    )
  )
  
  task_switching_clean
}


# add_response_var() ------------------------------------------------------

#' @description The purpose of this function is to add the `response` variable
#' which encoded the subject's response, which was  not recorded in the raw
#' data. By `response` I mean the answser a subject gives to the task request.
#' In each of the three blocks of trials (40 trials: food; 40 trials: plants; 
#' 40 trials: mixed) twere were two task requests. 
#' (1) Is the image a hyper-caloric food? Key "B": yes; key "N": no.
#' (2) Is the image a flower? Key "B": yes; key "N": no.
#' In the first two blocks, it is irrelevant whether the images were presented
#' in the upper or lower quadrant. In the "mixed" block, the subject was asked 
#' to perform task (1) when the images were presented in the upper quadrant; the
#' task (2) when the images were presented in the lower quadrant.
#' Sometimes, the subject was answering to different questions by pressing the 
#' same key. The variable `response` is codiding the kind of the subject's
#' response, in order to understand whether the response kind has changed or 
#' not with respect to the previous trial. I could define 4 different values for
#' the subject's response: 2 values x two questions. Or else I can reverse the 
#' key values when the subject perform different tasks (even though, in reality,
#' she/he pressed the same key when performing two different tasks). 
add_response_var <- function(df) {
  
  df <- df |> 
    dplyr::rename(
      stim = stim_category
    )
  
  df <- df %>%
    mutate(
      response = case_when(
        stim == "food" & food_img < 5 & is_correct == 1 ~ "B",
        stim == "food" & food_img > 4 & is_correct == 1 ~ "N",
        stim == "food" & food_img < 5 & is_correct == 0 ~ "N",
        stim == "food" & food_img > 4 & is_correct == 0 ~ "B",
        #
        stim == "plants" & plant_img < 5 & is_correct == 1 ~ "N",
        stim == "plants" & plant_img > 4 & is_correct == 1 ~ "B",
        stim == "plants" & plant_img < 5 & is_correct == 0 ~ "B",
        stim == "plants" & plant_img > 4 & is_correct == 0 ~ "N",
        #
        stim == "mixed" & (stim_position == 1 | stim_position == 2) & 
          food_img < 5 & is_correct == 1 ~ "B",
        stim == "mixed" & (stim_position == 1 | stim_position == 2) & 
          food_img > 4 & is_correct == 1 ~ "N",
        stim == "mixed" & (stim_position == 1 | stim_position == 2) & 
          food_img < 5 & is_correct == 0 ~ "N",
        stim == "mixed" & (stim_position == 1 | stim_position == 2) & 
          food_img > 4 & is_correct == 0 ~ "B",
        # 
        stim == "mixed" & (stim_position == 3 | stim_position == 4) & 
          plant_img < 5 & is_correct == 1 ~ "N",
        stim == "mixed" & (stim_position == 3 | stim_position == 4) & 
          plant_img > 4 & is_correct == 1 ~ "B",
        stim == "mixed" & (stim_position == 3 | stim_position == 4) & 
          plant_img < 5 & is_correct == 0 ~ "B",
        stim == "mixed" & (stim_position == 3 | stim_position == 4) & 
          plant_img > 4 & is_correct == 0 ~ "N",
      )
    )
  
  df$response <- factor(df$response)
  
  df
}


# wrangle_task_switch_data() ----------------------------------------------

wrangle_task_switch_data <- function(df) {
  
  # Define response in the previous trial 
  df <- df %>% 
    group_by(subj_id) %>% 
    mutate(prev_resp = lag(response))
  
  # Remove NAs (the first trial for each participant)
  df$resp_transition <- rep(NA, nrow(df))
  
  # Define response transition
  for (i in 1:nrow(df)) {
    if (is.na(df$prev_resp[i])) {
      df$resp_transition[i] = NA
    } else {
      df$resp_transition[i] = ifelse(
        df$response[i] == df$prev_resp[i], 
        "repetition", "switch"
      )
    }
  }
  
  # Remove NAs
  df <- df[!is.na(df$resp_transition), ]
  
  df <- df %>%
    dplyr::rename(
      task_transition = is_task_switch
    )
  
  df <- df %>% 
    mutate(task_transition = case_when(
      task_transition == 0 ~ 'repetition',
      task_transition == 1 ~ 'switch')
    )
  
  df <- df %>%
    mutate_at(
      c('resp_transition', 'task_transition', 'group'), 
    as.factor
  )
  # levels(task_switching_clean[['resp_transition']])
  
  saveRDS(
    df,
    here::here(
      "data", "processed", "task_switching", 
      "task_switch_data_for_descript_stats.rds"
    )
  )
}
