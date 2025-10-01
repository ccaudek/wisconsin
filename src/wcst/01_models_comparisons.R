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

num_cores <- parallel::detectCores() # Automatically detect the number of cores
# Ensure we don't exceed the number of available cores
parallel_chains <- min(4, num_cores) # Limit to 4 chains or available cores


# Bishara model -----------------------------------------------------------
#
# file <- here::here(
#   "src",
#   "wcst",
#   "documentation",
#   "stan",
#   "bishara_1.stan"
# )
#
# # Compile model
# mod <- cmdstan_model(file)
#
# # Fit con sampling (più lento ma stime complete)
# fit <- mod$sample(
#   data = stan_data,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 2000,
#   refresh = 100,
#   max_treedepth = 12,
#   adapt_delta = 0.95
# )
#
# # OPPURE: Fit con optimization (molto più veloce, stime MAP)
# fit_opt <- mod$optimize(
#   data = stan_data,
#   algorithm = "lbfgs",
#   init = function() list(
#     r = runif(stan_data$N, 0.3, 0.9),
#     p = runif(stan_data$N, 0.3, 0.9),
#     d = runif(stan_data$N, 0.5, 2)
#   )
# )
#
# # Estrai parametri
# params_summary <- fit$summary(c("r", "p", "d"))
# print(params_summary)
#
# params_df <- fit$summary(c("r", "p", "d")) %>%
#   filter(grepl("\\[", variable)) %>%
#   mutate(
#     param = gsub("\\[.*", "", variable),
#     subject = as.numeric(gsub(".*\\[|\\]", "", variable))
#   ) %>%
#   left_join(
#     data.frame(subject = 1:88, group = stan_data$group),
#     by = "subject"
#   )
#
# # Plot r vs p colorato per gruppo
# ggplot(params_df %>% pivot_wider(id_cols = c(subject, group),
#                                  names_from = param,
#                                  values_from = mean),
#        aes(x = r, y = p, color = factor(group))) +
#   geom_point(size = 3, alpha = 0.7) +
#   labs(title = "Bishara Model: Reward vs Punishment Learning",
#        x = "r (reward learning rate)",
#        y = "p (punishment learning rate)",
#        color = "Group") +
#   theme_minimal()
#
# # Confronto con parametri precedenti
# cat("\nVariabilità tra soggetti:\n")
# cat("r  SD:", sd(params_df$mean[params_df$param == "r"]), "\n")
# cat("p  SD:", sd(params_df$mean[params_df$param == "p"]), "\n")
# cat("d  SD:", sd(params_df$mean[params_df$param == "d"]), "\n")
#
#
# # Prepara i dati aggiungendo etichette gruppo
# params_df <- params_df %>%
#   mutate(
#     group_label = factor(group,
#                          levels = c(0, 1),
#                          labels = c("Controls", "Patients"))
#   )
#
# # Boxplot per i tre parametri
# ggplot(params_df, aes(x = group_label, y = mean, fill = group_label)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = 16) +
#   geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
#   facet_wrap(~ param, scales = "free_y", ncol = 3,
#              labeller = labeller(param = c(
#                "r" = "r (Reward Learning)",
#                "p" = "p (Punishment Learning)",
#                "d" = "d (Decision Consistency)"
#              ))) +
#   labs(title = "Bishara WCST Model Parameters by Group",
#        x = "Group",
#        y = "Parameter Estimate",
#        fill = "Group") +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     strip.text = element_text(size = 11, face = "bold")
#   ) +
#   scale_fill_manual(values = c("Controls" = "#4CAF50", "Patients" = "#F44336"))
#
# # Test statistici
# cat("\n=== STATISTICAL TESTS ===\n\n")
#
# for (par in c("r", "p", "d")) {
#   cat(sprintf("Parameter: %s\n", par))
#
#   param_data <- params_df %>% filter(param == par)
#
#   # t-test
#   test_result <- t.test(mean ~ group, data = param_data)
#
#   # Cohen's d
#   controls <- param_data$mean[param_data$group == 1]
#   patients <- param_data$mean[param_data$group == 0]
#
#   pooled_sd <- sqrt(((length(controls)-1)*sd(controls)^2 +
#                        (length(patients)-1)*sd(patients)^2) /
#                       (length(controls) + length(patients) - 2))
#   cohens_d <- (mean(controls) - mean(patients)) / pooled_sd
#
#   cat(sprintf("  Controls: M = %.3f, SD = %.3f\n",
#               mean(controls), sd(controls)))
#   cat(sprintf("  Patients: M = %.3f, SD = %.3f\n",
#               mean(patients), sd(patients)))
#   cat(sprintf("  t(%.0f) = %.3f, p = %.4f\n",
#               test_result$parameter, test_result$statistic, test_result$p.value))
#   cat(sprintf("  Cohen's d = %.3f\n\n", cohens_d))
# }
#
# # Summary statistics per gruppo
# cat("\n=== SUMMARY BY GROUP ===\n\n")
# params_summary <- params_df %>%
#   group_by(param, group_label) %>%
#   summarise(
#     n = n(),
#     mean = mean(mean),
#     median = median(mean),
#     sd = sd(mean),
#     min = min(mean),
#     max = max(mean),
#     .groups = "drop"
#   ) %>%
#   arrange(param, group_label)
#
# print(params_summary, n = 100)

######### ----------------------------------------------------------------------

file <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "bishara_single_subject.stan"
)

model_single <- cmdstan_model(file)

# Prepara contenitore per i risultati
all_params <- data.frame(
  subject = integer(),
  r = numeric(),
  p = numeric(),
  d = numeric(),
  log_lik = numeric(),
  converged = logical()
)

# Fit per ogni soggetto
cat("Fitting individual subjects...\n")

for (subj in 1:stan_data$N) {
  # Prepara dati per questo soggetto
  n_trials <- stan_data$Tsubj[subj]

  subj_data <- list(
    T = n_trials,
    rew = stan_data$rew[subj, 1:n_trials],
    los = stan_data$los[subj, 1:n_trials],
    resp_choice = stan_data$resp_choice[subj, 1:n_trials],
    resp_color = stan_data$resp_color[subj, 1:n_trials],
    resp_shape = stan_data$resp_shape[subj, 1:n_trials],
    resp_number = stan_data$resp_number[subj, 1:n_trials]
  )

  # Fit con optimization (MAP estimation, come nel paper)
  fit <- model_single$optimize(
    data = subj_data,
    algorithm = "lbfgs",
    init = function()
      list(
        r = runif(1, 0.3, 0.7),
        p = runif(1, 0.3, 0.7),
        d = runif(1, 0.5, 2)
      ),
    jacobian = TRUE,
    refresh = 0
  )

  # Estrai parametri
  params <- fit$mle()

  all_params <- rbind(
    all_params,
    data.frame(
      subject = subj,
      r = params["r"],
      p = params["p"],
      d = params["d"],
      log_lik = fit$summary()$estimate[fit$summary()$variable == "log_lik"],
      converged = fit$return_codes() == 0
    )
  )

  if (subj %% 10 == 0) {
    cat(sprintf("  Fitted %d/%d subjects\n", subj, stan_data$N))
  }
}

# Aggiungi informazione sul gruppo
all_params <- all_params %>%
  mutate(
    group = stan_data$group[subject],
    group_label = factor(
      group,
      levels = c(1, 2),
      labels = c("Controls", "Patients")
    )
  )

# Verifica convergenza
cat("\nConvergence check:\n")
cat(sprintf(
  "  Converged: %d/%d subjects\n",
  sum(all_params$converged),
  nrow(all_params)
))

# Statistiche descrittive
cat("\n=== VARIABILITY ACROSS SUBJECTS ===\n")
cat(sprintf(
  "r  - Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
  mean(all_params$r),
  sd(all_params$r),
  min(all_params$r),
  max(all_params$r)
))
cat(sprintf(
  "p  - Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
  mean(all_params$p),
  sd(all_params$p),
  min(all_params$p),
  max(all_params$p)
))
cat(sprintf(
  "d  - Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
  mean(all_params$d),
  sd(all_params$d),
  min(all_params$d),
  max(all_params$d)
))

# Test statistici tra gruppi
cat("\n=== GROUP COMPARISONS ===\n\n")

for (par in c("r", "p", "d")) {
  param_data <- all_params[[par]]
  controls <- param_data[all_params$group == 1]
  patients <- param_data[all_params$group == 2]

  test <- t.test(param_data ~ all_params$group)

  pooled_sd <- sqrt(
    ((length(controls) - 1) *
      sd(controls)^2 +
      (length(patients) - 1) * sd(patients)^2) /
      (length(controls) + length(patients) - 2)
  )
  cohens_d <- (mean(controls) - mean(patients)) / pooled_sd

  cat(sprintf("Parameter: %s\n", par))
  cat(sprintf(
    "  Controls: M = %.3f (SD = %.3f)\n",
    mean(controls),
    sd(controls)
  ))
  cat(sprintf(
    "  Patients: M = %.3f (SD = %.3f)\n",
    mean(patients),
    sd(patients)
  ))
  cat(sprintf(
    "  t(%.0f) = %.3f, p = %.4f\n",
    test$parameter,
    test$statistic,
    test$p.value
  ))
  cat(sprintf("  Cohen's d = %.3f\n", cohens_d))

  if (test$p.value < 0.05) {
    cat(sprintf("  *** SIGNIFICANT at p < .05 ***\n"))
  }
  cat("\n")
}

# Visualizzazioni
# 1. Boxplot per parametro
params_long <- all_params %>%
  tidyr::pivot_longer(
    cols = c(r, p, d),
    names_to = "parameter",
    values_to = "estimate"
  )

p1 <- ggplot(
  params_long,
  aes(x = group_label, y = estimate, fill = group_label)
) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(
    ~parameter,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(
      parameter = c(
        "r" = "r (Reward Learning)",
        "p" = "p (Punishment Learning)",
        "d" = "d (Decision Consistency)"
      )
    )
  ) +
  labs(
    title = "Bishara WCST Model - Individual Subject Fits",
    subtitle = "MAP estimates via optimization",
    x = "Group",
    y = "Parameter Estimate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold")
  ) +
  scale_fill_manual(values = c("Patients" = "#F44336", "Controls" = "#4CAF50"))

print(p1)

# 2. Scatter plot r vs p
p2 <- ggplot(all_params, aes(x = r, y = p, color = group_label)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  labs(
    title = "Reward vs Punishment Learning",
    x = "r (Reward Learning Rate)",
    y = "p (Punishment Learning Rate)",
    color = "Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Patients" = "#F44336", "Controls" = "#4CAF50"))

print(p2)

# 3. Salva risultati
write.csv(all_params, "bishara_individual_params.csv", row.names = FALSE)
cat("\nResults saved to 'bishara_individual_params.csv'\n")


# Rescorla-Wagner dual learning model ------------------------------------------

file <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "08_rw_dual_learning.stan"
)

# Compile model
mod <- cmdstan_model(file)


# Fit per ogni soggetto
all_params <- data.frame(
  subject = integer(),
  alpha_pos = numeric(),
  alpha_neg = numeric(),
  beta = numeric(),
  log_lik = numeric(),
  converged = logical()
)

cat("Fitting RW model with dual learning rates...\n")

for (subj in 1:stan_data$N) {
  n_trials <- stan_data$Tsubj[subj]

  subj_data <- list(
    T = n_trials,
    rew = stan_data$rew[subj, 1:n_trials],
    los = stan_data$los[subj, 1:n_trials],
    resp_choice = stan_data$resp_choice[subj, 1:n_trials],
    resp_color = stan_data$resp_color[subj, 1:n_trials],
    resp_shape = stan_data$resp_shape[subj, 1:n_trials],
    resp_number = stan_data$resp_number[subj, 1:n_trials]
  )

  fit <- mod$optimize(
    data = subj_data,
    algorithm = "lbfgs",
    # init = list(
    #   alpha_pos = runif(1, 0.3, 0.7),
    #   alpha_neg = runif(1, 0.3, 0.7),
    #   beta = runif(1, 0.8, 1.5)
    # ),
    refresh = 0
  )

  params <- fit$mle()

  all_params <- rbind(
    all_params,
    data.frame(
      subject = subj,
      alpha_pos = params["alpha_pos"],
      alpha_neg = params["alpha_neg"],
      beta = params["beta"],
      log_lik = fit$summary()$estimate[fit$summary()$variable == "log_lik"],
      converged = (fit$return_codes() == 0)
    )
  )

  if (subj %% 10 == 0) {
    cat(sprintf("  %d/%d subjects\n", subj, stan_data$N))
  }
}

# Aggiungi gruppo
all_params <- all_params %>%
  mutate(
    group = stan_data$group[subject],
    group_label = factor(
      group,
      levels = c(1, 2),
      labels = c("Controls", "Patients")
    )
  )

# Statistiche
cat("\n=== VARIABILITY ===\n")
cat(sprintf(
  "alpha_pos: M=%.3f, SD=%.3f\n",
  mean(all_params$alpha_pos),
  sd(all_params$alpha_pos)
))
cat(sprintf(
  "alpha_neg: M=%.3f, SD=%.3f\n",
  mean(all_params$alpha_neg),
  sd(all_params$alpha_neg)
))
cat(sprintf(
  "beta:      M=%.3f, SD=%.3f\n",
  mean(all_params$beta),
  sd(all_params$beta)
))

# Test tra gruppi
cat("\n=== GROUP COMPARISONS ===\n\n")

for (par in c("alpha_pos", "alpha_neg", "beta")) {
  vals <- all_params[[par]]
  controls <- vals[all_params$group == 1]
  patients <- vals[all_params$group == 2]

  test <- t.test(controls, patients)
  cohens_d <- (mean(controls) - mean(patients)) /
    sqrt(
      ((length(controls) - 1) *
        sd(controls)^2 +
        (length(patients) - 1) * sd(patients)^2) /
        (length(controls) + length(patients) - 2)
    )

  cat(sprintf("%s:\n", par))
  cat(sprintf("  Controls: %.3f (%.3f)\n", mean(controls), sd(controls)))
  cat(sprintf("  Patients: %.3f (%.3f)\n", mean(patients), sd(patients)))
  cat(sprintf(
    "  t=%.2f, p=%.4f, d=%.2f",
    test$statistic,
    test$p.value,
    cohens_d
  ))
  if (test$p.value < 0.05) cat(" ***")
  cat("\n\n")
}

# Plot
params_long <- all_params %>%
  pivot_longer(
    cols = c(alpha_pos, alpha_neg, beta),
    names_to = "parameter",
    values_to = "estimate"
  )

ggplot(params_long, aes(x = group_label, y = estimate, fill = group_label)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~parameter, scales = "free_y") +
  labs(title = "RW Model - Dual Learning Rates") +
  theme_minimal() +
  scale_fill_manual(values = c("Controls" = "#4CAF50", "Patients" = "#F44336"))


#########

file <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "09_rw_dual_learning.stan"
)

# Compile model
mod <- cmdstan_model(file)

fit <- mod$sample(
  data = subj_data,
  seed = 123,
  chains = 4,
  parallel_chains = 2,
  iter_warmup = 1000,
  iter_sampling = 2000,
  refresh = 0
)

summ <- fit$summary(c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse"))
summ

# options(mc.cores = parallel::detect_cores())

# --- INPUT ---
# - 'mod' è il modello già compilato (cmdstan_model)
# - 'stan_data' contiene le matrici/vettori completi
#   campi minimi richiesti per per-soggetto: Tsubj, rew, los, resp_choice, resp_color, resp_shape, resp_number
# - 'stan_data$group' è un vettore (1=Controls, 2=Patients)

# Scegli la statistica riassuntiva da usare dal posterior
stat_col <- "median" # "median" oppure "mean"

# helper: estrae una riga (tibble) con i riassunti posterior per un soggetto
summarize_subject_fit <- function(
  fit,
  subj_id,
  stat_col = c("median", "mean")
) {
  stat_col <- match.arg(stat_col)
  pars <- c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")
  S <- fit$summary(pars)
  # riga unica: portiamo in wide con le colonne scelte (median o mean)
  val_col <- if (stat_col == "median") "median" else "mean"
  out <- S %>%
    select(variable, all_of(val_col), rhat, ess_bulk, ess_tail) %>%
    tidyr::pivot_wider(names_from = variable, values_from = all_of(val_col)) %>%
    mutate(
      subject = subj_id,
      converged = ifelse(all(is.finite(rhat) & rhat < 1.01), TRUE, NA)
    ) # flag semplice
  out
}

# loop sui soggetti
fit_one_subject <- function(
  subj,
  mod,
  stan_data,
  iter_warmup = 800,
  iter_sampling = 800,
  seed = 123
) {
  n_trials <- stan_data$Tsubj[subj]
  subj_data <- list(
    T = n_trials,
    rew = stan_data$rew[subj, 1:n_trials],
    los = stan_data$los[subj, 1:n_trials],
    resp_choice = stan_data$resp_choice[subj, 1:n_trials],
    resp_color = stan_data$resp_color[subj, 1:n_trials],
    resp_shape = stan_data$resp_shape[subj, 1:n_trials],
    resp_number = stan_data$resp_number[subj, 1:n_trials]
  )
  fit <- mod$sample(
    data = subj_data,
    seed = seed + subj,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  fit
}

message("Fitting RW dual-alpha + sticky + lapse per subject...")
N <- stan_data$N

all_params <- map_dfr(1:N, function(subj) {
  fit_s <- fit_one_subject(subj, mod, stan_data)
  summ <- summarize_subject_fit(fit_s, subj, stat_col = stat_col)
  # aggiungo anche log_lik totale se vuoi (da generated quantities)
  # NB: per come è definito, il modello ha 'log_lik' scalare: usiamo summary
  ll <- fit_s$summary("log_lik")$median
  summ$log_lik <- ll
  summ
})

# Aggiungi gruppo e label
all_params <- all_params %>%
  mutate(
    group = stan_data$group[subject],
    group_label = factor(
      group,
      levels = c(1, 2),
      labels = c("Controls", "Patients")
    )
  )

# Stampa statistiche globali
cat("\n=== VARIABILITY (", stat_col, ") ===\n", sep = "")
cat(sprintf(
  "alpha_pos: M=%.3f, SD=%.3f\n",
  mean(all_params$alpha_pos, na.rm = TRUE),
  sd(all_params$alpha_pos, na.rm = TRUE)
))
cat(sprintf(
  "alpha_neg: M=%.3f, SD=%.3f\n",
  mean(all_params$alpha_neg, na.rm = TRUE),
  sd(all_params$alpha_neg, na.rm = TRUE)
))
cat(sprintf(
  "beta:      M=%.3f, SD=%.3f\n",
  mean(all_params$beta, na.rm = TRUE),
  sd(all_params$beta, na.rm = TRUE)
))
cat(sprintf(
  "kappa:     M=%.3f, SD=%.3f\n",
  mean(all_params$kappa, na.rm = TRUE),
  sd(all_params$kappa, na.rm = TRUE)
))
cat(sprintf(
  "lapse:     M=%.3f, SD=%.3f\n",
  mean(all_params$lapse, na.rm = TRUE),
  sd(all_params$lapse, na.rm = TRUE)
))

# Long per il plot
params_long <- all_params %>%
  select(subject, group_label, alpha_pos, alpha_neg, beta, kappa, lapse) %>%
  pivot_longer(
    cols = c(alpha_pos, alpha_neg, beta, kappa, lapse),
    names_to = "parameter",
    values_to = "estimate"
  )

# Boxplot + jitter per gruppo, faccettato per parametro
p <- ggplot(
  params_long,
  aes(x = group_label, y = estimate, fill = group_label)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.6) +
  facet_wrap(~parameter, scales = "free_y") +
  labs(
    title = paste(
      "RW dual-α + sticky + lapse — per-subject",
      toupper(stat_col)
    ),
    x = NULL,
    y = "Posterior summary"
  ) +
  scale_fill_manual(
    values = c("Controls" = "#4CAF50", "Patients" = "#F44336")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p)

# (opzionale) salva CSV e figura
# readr::write_csv(all_params, "all_params_per_subject.csv")
# ggsave("params_boxplot_by_group.png", p, width = 10, height = 6, dpi = 300)

###### eof #################

fit <- model$optimize(
  data = subj_data,
  algorithm = "lbfgs",
  init = list(
    alpha_pos = runif(1, 0.3, 0.7),
    alpha_neg = runif(1, 0.3, 0.7),
    beta = runif(1, 0.8, 1.5)
  ),
  refresh = 0
)


# Rescorla-Wagner model ---------------------------------------------------

file <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "05_model_no_hierarchy.stan"
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

fit_indep <- mod$variational(
  data = stan_data,
  seed = 1234
)

# Save fit.
fit1$save_object(
  file = here::here("src", "wcst", "fits", "fit1.RDS")
)

alpha_summ <- fit1$summary("alpha")
beta_summ <- fit1$summary("beta")
inertia_summ <- fit1$summary("inertia")

alpha_pos_summ <- fit1$summary("alpha_pos")
alpha_neg_summ <- fit1$summary("alpha_neg")
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
  family = binomial(link = "logit"),
  data = mydat
)

test_prob <- predict(fm, newdata = mydat, type = "response")
# Compute ROC curve and AUC using pROC
roc_obj <- roc(mydat$is_patient, test_prob, plot = TRUE, print.auc = TRUE)

t.test(beta ~ is_patient, data = mydat)


# Empirical Bayes ----------------------

file_no_hier <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "05_model_no_hierarchy.stan"
)

mod_no_hier <- cmdstan_model(file_no_hier)

# stan_data: deve contenere i campi attesi dal tuo .stan (N, T, Tsubj, holdout, ecc.)
fit_indep <- mod_no_hier$sample(
  data = stan_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)


# (Facoltativo) alternativa VI, sconsigliata per EB perché può introdurre bias:
# fit_indep <- mod_no_hier$variational(data = stan_data, seed = 1234)

# ── 2) Estrai i posterior e calcola statistiche EB per i prior ────────────────
# Helper per estrarre matrice draws (iterazioni x soggetti) di un parametro vettoriale
get_param_matrix <- function(fit, par_name) {
  # estrae solo le colonne tipo "alpha[1]", "alpha[2]"...
  da <- as_draws_df(fit$draws(par_name))
  # converti in matrix ignorando eventuali colonne 'lp__' ecc.
  keep <- grep(sprintf("^%s\\[", par_name), names(da))
  as.matrix(da[, keep, drop = FALSE])
}

alpha_mat <- get_param_matrix(fit_indep, "alpha")
beta_mat <- get_param_matrix(fit_indep, "beta")
inert_mat <- get_param_matrix(fit_indep, "inertia")

# Posterior mean per SOGGETTO (colMeans su ogni matrice di draws)
alpha_bar_by_subj <- colMeans(alpha_mat)
beta_bar_by_subj <- colMeans(beta_mat)
inert_bar_by_subj <- colMeans(inert_mat)

# ── 2a) Hyperparam Beta(a,b) per alpha in [0,1] via method-of-moments ─────────
# MOM per Beta: dato mean m e sd s (s sulla distribuzione tra soggetti),
#   a = m * ((m*(1-m)/s^2) - 1)
#   b = (1-m) * ((m*(1-m)/s^2) - 1)
# Attenzione: serve s^2 < m*(1-m); in caso contrario clamp s.

m_alpha <- mean(alpha_bar_by_subj)
s_alpha <- sd(alpha_bar_by_subj)
# clamp s per stabilità numerica
s2_max <- m_alpha * (1 - m_alpha) - 1e-6
s_alpha <- min(s_alpha, sqrt(max(s2_max, 1e-8)))

mm_term <- (m_alpha * (1 - m_alpha) / (s_alpha^2)) - 1
a_alpha <- m_alpha * mm_term
b_alpha <- (1 - m_alpha) * mm_term

# fallback di sicurezza (se dati estremi)
if (
  !is.finite(a_alpha) || a_alpha <= 0 || !is.finite(b_alpha) || b_alpha <= 0
) {
  # usa un Beta che concentra intorno a m_alpha con moderata certezza
  a_alpha <- max(2 * m_alpha, 0.5)
  b_alpha <- max(2 * (1 - m_alpha), 0.5)
}

# ── 2b) Hyperparam Lognormal per beta >0 e inertia >0 ─────────────────────────
# Fit su log dei posterior means per soggetto
mu_log_beta <- mean(log(beta_bar_by_subj))
sigma_log_beta <- sd(log(beta_bar_by_subj))
mu_log_inertia <- mean(log(inert_bar_by_subj))
sigma_log_inertia <- sd(log(inert_bar_by_subj))

# clamp minimi per stabilità
sigma_log_beta <- max(sigma_log_beta, 1e-3)
sigma_log_inertia <- max(sigma_log_inertia, 1e-3)

empirical_stats <- list(
  m_alpha = m_alpha,
  s_alpha = s_alpha,
  a_alpha = a_alpha,
  b_alpha = b_alpha,
  mu_log_beta = mu_log_beta,
  sigma_log_beta = sigma_log_beta,
  mu_log_inertia = mu_log_inertia,
  sigma_log_inertia = sigma_log_inertia
)
print(empirical_stats)

# ── 3) Prepara i dati per il MODELLO GERARCHICO con prior EB ──────────────────
# Assumiamo che il modello gerarchico prenda questi iper-parametri come `data`
stan_data_hier <- stan_data
stan_data_hier$a_alpha <- a_alpha
stan_data_hier$b_alpha <- b_alpha
stan_data_hier$mu_log_beta <- mu_log_beta
stan_data_hier$sigma_log_beta <- sigma_log_beta
stan_data_hier$mu_log_inertia <- mu_log_inertia
stan_data_hier$sigma_log_inertia <- sigma_log_inertia

# ── 4) Compila e fitta il modello GERARCHICO ──────────────────────────────────
# Sostituisci con il percorso reale del tuo modello gerarchico:
file_ml <- here::here(
  "src",
  "wcst",
  "documentation",
  "stan",
  "07_indep_flat_priors.stan"
)
mod_ml <- cmdstan_model(file_ml)

fit_ml <- mod_hier$sample(
  data = stan_data_hier,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500
)

# ── 5) Verifiche post-fitting richieste ───────────────────────────────────────
# 5.1) Confronto variabilità tra soggetti (posterior mean per soggetto)
alpha_mat_hier <- get_param_matrix(fit_ml, "alpha")
alpha_bar_by_subj_hier <- colMeans(alpha_mat_hier)

cat("SD between subjects (independent): ", sd(alpha_bar_by_subj), "\n")
cat("SD between subjects (hierarchical):", sd(alpha_bar_by_subj_hier), "\n")

# 5.2) Stampa degli iper-sigma del modello gerarchico
# Adatta i nomi dei parametri agli effettivi nel tuo .stan:
#   es. "tau_alpha", "tau_beta", "tau_inertia" oppure "sigma_alpha", ...
pars_sigma <- c(
  "sigma_alpha",
  "sigma_beta",
  "sigma_inertia",
  "tau_alpha",
  "tau_beta",
  "tau_inertia"
)
pars_sigma <- pars_sigma[pars_sigma %in% fit_hier$metadata()$parameter_names]
if (length(pars_sigma)) {
  print(fit_hier$summary(pars_sigma), n = Inf)
} else {
  message(
    "Nessun parametro di scala (sigma/tau) trovato con quei nomi: aggiorna i nomi."
  )
}

# ── 6) (Facoltativo) QC grafico: shrinkage ────────────────────────────────────
df_shrink <- tibble(
  subj = seq_along(alpha_bar_by_subj),
  alpha_indep = alpha_bar_by_subj,
  alpha_hier = alpha_bar_by_subj_hier
)

ggplot(df_shrink, aes(x = alpha_indep, y = alpha_hier)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.6) +
  labs(
    x = "Posterior mean alpha (independent)",
    y = "Posterior mean alpha (hierarchical)",
    title = "Shrinkage check: alpha (per subject)"
  ) +
  theme_minimal()


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
  family = binomial(link = "logit"),
  data = mydat
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
  "src",
  "wcst",
  "documentation",
  "stan",
  "02_rescorla_wagner_inertia.stan"
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
  family = binomial(link = "logit"),
  data = mydat
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
