# ============================================================
# RW dual-α + sticky + lapse (empirical WCST) — full pipeline
# - Usa le tue funzioni per preparare i dati
# - Scrive/compila il modello Stan migliorato
# - Fit per-soggetto (optimize o sample)
# - Estrazione parametri SENZA NA + grafici
# ============================================================

# 0) Packages ---------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
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

set.seed(20251001)

# 1) Sorgenti del progetto --------------------------------------------------
source(here::here("src", "wcst", "documentation", "functions", "funs_wcst.R"))
source(here::here(
  "src",
  "wcst",
  "documentation",
  "functions",
  "funs_model_selection_wcst.R"
))
source(here::here(
  "src",
  "wcst",
  "documentation",
  "functions",
  "funs_input_for_stan_wcst.R"
))

# 2) Dati empirici -> stan_data --------------------------------------------
generate_csv("controls")
generate_csv("patients")

stan_data <- process_and_prepare_stan_data()
str(stan_data)

# Aspettative:
# stan_data$N, $T, $Tsubj[1:N], $group[1:N],
# $rew[N,T], $los[N,T] (NB: los in {-1,0}! -> convertiamo a {1,0}),
# $resp_choice[N,T] in {1..4}, $resp_color/shape/number (extra)
stopifnot(
  is.integer(stan_data$resp_choice),
  is.integer(stan_data$rew),
  is.integer(stan_data$los)
)

# 3) Scrivi/compila il modello Stan ----------------------------------------
stan_dir <- here::here("src", "wcst", "documentation", "stan")
if (!dir.exists(stan_dir)) dir.create(stan_dir, recursive = TRUE)
stan_file <- file.path(stan_dir, "10_rw_dual_sticky_lapse.stan")

stan_code <- '
data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> choice;
  array[T] int<lower=0, upper=1> rew;
  array[T] int<lower=0, upper=1> los;
  array[T] int<lower=0, upper=1> resp_color;
  array[T] int<lower=0, upper=1> resp_shape;
  array[T] int<lower=0, upper=1> resp_number;
}
parameters {
  real alpha_pos_raw;   // (0,1) via inv_logit
  real alpha_neg_raw;   // (0,1) via inv_logit
  real beta_raw;        // (0,+inf) via exp
  real kappa;           // sticky bias
  real lapse_raw;       // (0,1) via inv_logit
}
transformed parameters {
  real<lower=0, upper=1> alpha_pos = inv_logit(alpha_pos_raw);
  real<lower=0, upper=1> alpha_neg = inv_logit(alpha_neg_raw);
  real<lower=0>          beta      = exp(beta_raw);
  real<lower=0, upper=1> lapse     = inv_logit(lapse_raw);
}
model {
  // Priors debolmente informativi
  alpha_pos_raw ~ normal(0, 1);
  alpha_neg_raw ~ normal(0, 1);
  beta_raw      ~ normal(0.0, 0.7);
  kappa         ~ student_t(3, 0, 1);
  lapse_raw     ~ normal(-2.2, 1.0);  // media ~ 0.10 su scala naturale

  {
    vector[4] Q = rep_vector(0.0, 4);
    int prev_choice = 0;

    for (t in 1:T) {
      vector[4] logits = beta * Q;
      if (prev_choice > 0) logits[prev_choice] += kappa;

      // Mixture con lapse (uniforme su 4 azioni)
      target += log_mix(
        1.0 - lapse,
        categorical_logit_lpmf(choice[t] | logits),
        -log(4.0)
      );

      // Outcome in {-1, +1} da rew/los binari
      {
        int a = choice[t];
        real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
        real pe = r - Q[a];
        if (pe >= 0) Q[a] += alpha_pos * pe;
        else         Q[a] += alpha_neg * pe;
        prev_choice = a;
      }
    }
  }
}
generated quantities {
  real log_lik;
  {
    vector[4] Q = rep_vector(0.0, 4);
    int prev_choice = 0;
    real lp = 0;
    for (t in 1:T) {
      vector[4] logits = exp(beta_raw) * Q;
      if (prev_choice > 0) logits[prev_choice] += kappa;
      lp += log_mix(
        1.0 - inv_logit(lapse_raw),
        categorical_logit_lpmf(choice[t] | logits),
        -log(4.0)
      );
      {
        int a = choice[t];
        real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
        real pe = r - Q[a];
        if (pe >= 0) Q[a] += inv_logit(alpha_pos_raw) * pe;
        else         Q[a] += inv_logit(alpha_neg_raw) * pe;
        prev_choice = a;
      }
    }
    log_lik = lp;
  }
}
'
writeLines(stan_code, con = stan_file)

if (is.null(cmdstanr::cmdstan_path())) {
  message("CmdStan non trovato: avvio installazione (richiede toolchain C++).")
  cmdstanr::install_cmdstan()
}
mod <- cmdstanr::cmdstan_model(stan_file)

# 4) Scelte di fitting ------------------------------------------------------
fit_mode <- "sample" # "optimize" (MAP) oppure "sample" (MCMC)
iter_warmup <- 1000
iter_sampling <- 2000
chains <- 4
parallel_chains <- min(chains, max(1, parallel::detectCores() - 1))

# 5) Helper: prepara dati per SOGGETTO i -----------------------------------
prep_subject_data <- function(i) {
  Ti <- stan_data$Tsubj[i]
  stopifnot(Ti <= ncol(stan_data$resp_choice))

  # converti los {-1,0} -> {1,0} e rew in {0,1}
  los_bin <- as.integer(stan_data$los[i, 1:Ti] == -1L)
  rew_bin <- as.integer(stan_data$rew[i, 1:Ti] == 1L)

  list(
    T = Ti,
    choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    rew = rew_bin,
    los = los_bin,
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti] == 1L),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti] == 1L),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti] == 1L)
  )
}

# 6) Fit per soggetto (due versioni) ---------------------------------------
fit_one_subject_opt <- function(i) {
  dat <- prep_subject_data(i)

  fit <- mod$optimize(
    data = dat,
    algorithm = "lbfgs",
    init = function()
      list(
        alpha_pos_raw = rnorm(1, 0.4, 0.6),
        alpha_neg_raw = rnorm(1, -0.6, 0.6),
        beta_raw = rnorm(1, 0.2, 0.5),
        kappa = rnorm(1, 0.0, 0.7),
        lapse_raw = rnorm(1, -2.0, 0.7)
      ),
    jacobian = TRUE,
    refresh = 0
  )

  pars <- fit$mle() # named numeric vector
  sm <- fit$summary("log_lik")
  ll <- sm$estimate[sm$variable == "log_lik"]

  tibble(
    subject = i,
    alpha_pos = unname(pars["alpha_pos"]),
    alpha_neg = unname(pars["alpha_neg"]),
    beta = unname(pars["beta"]),
    kappa = unname(pars["kappa"]),
    lapse = unname(pars["lapse"]),
    log_lik = ll,
    converged = (fit$return_codes() == 0)
  )
}

fit_one_subject_mcmc <- function(
  i,
  chains,
  parallel_chains,
  iter_warmup,
  iter_sampling
) {
  dat <- prep_subject_data(i)

  fit <- mod$sample(
    data = dat,
    seed = 20251001 + i,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )

  S <- fit$summary(c(
    "alpha_pos",
    "alpha_neg",
    "beta",
    "kappa",
    "lapse",
    "log_lik"
  ))

  tibble(
    subject = i,
    alpha_pos = S$median[S$variable == "alpha_pos"],
    alpha_neg = S$median[S$variable == "alpha_neg"],
    beta = S$median[S$variable == "beta"],
    kappa = S$median[S$variable == "kappa"],
    lapse = S$median[S$variable == "lapse"],
    log_lik = S$median[S$variable == "log_lik"],
    converged = all(
      S$rhat[
        S$variable %in% c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")
      ] <
        1.05
    )
  )
}

# 7) Runner + esecuzione sicura --------------------------------------------
# evita eventuali omonimie pregresse
if (
  exists("fit_one_subject", inherits = FALSE) &&
    !is.function(get("fit_one_subject"))
) {
  rm(fit_one_subject)
}

message("Fitting per-soggetto: mode = ", fit_mode)
N <- stan_data$N

runner <- if (identical(fit_mode, "optimize")) {
  fit_one_subject_opt
} else {
  \(i)
    fit_one_subject_mcmc(
      i,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling
    )
}

# (facoltativo) dry-run su primo soggetto
print(runner(1))

# esecuzione con handling degli errori (continua sui soggetti validi)
safe_runner <- purrr::safely(runner, otherwise = NULL, quiet = FALSE)
res_list <- purrr::map(seq_len(N), safe_runner)

fail_idx <- which(purrr::map_lgl(res_list, ~ !is.null(.x$error)))
if (length(fail_idx)) {
  message("Soggetti falliti: ", paste(fail_idx, collapse = ", "))
  message("Primo errore: ", res_list[[fail_idx[1]]]$error$message)
}

all_params <- res_list |>
  purrr::map("result") |>
  purrr::compact() |>
  dplyr::bind_rows()

stopifnot(nrow(all_params) > 0L)

# 8) Aggiungi gruppo e controlli NA ----------------------------------------
all_params <- all_params %>%
  mutate(
    group = stan_data$group[subject],
    group_label = factor(
      group,
      levels = c(1, 2),
      labels = c("Controls", "Patients")
    )
  )

stopifnot(all(is.finite(all_params$alpha_pos)))
stopifnot(all(is.finite(all_params$alpha_neg)))
stopifnot(all(is.finite(all_params$beta)))
stopifnot(all(is.finite(all_params$kappa)))
stopifnot(all(is.finite(all_params$lapse)))

# 9) Sommari e plot ---------------------------------------------------------
cat(
  "\n=== SUMMARY (",
  toupper(ifelse(fit_mode == "sample", "posterior median", "MAP")),
  ") ===\n",
  sep = ""
)
for (par in c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")) {
  cat(sprintf(
    "%-9s M=%.3f  SD=%.3f  min=%.3f  max=%.3f\n",
    par,
    mean(all_params[[par]]),
    sd(all_params[[par]]),
    min(all_params[[par]]),
    max(all_params[[par]])
  ))
}

params_long <- all_params %>%
  select(subject, group_label, alpha_pos, alpha_neg, beta, kappa, lapse) %>%
  pivot_longer(
    cols = c(alpha_pos, alpha_neg, beta, kappa, lapse),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")
    )
  )

p <- ggplot(
  params_long,
  aes(x = group_label, y = estimate, fill = group_label)
) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.6) +
  facet_wrap(~parameter, scales = "free_y") +
  labs(
    title = paste0("RW dual-α + sticky + lapse — per-subject (", fit_mode, ")"),
    x = NULL,
    y = "Estimate"
  ) +
  scale_fill_manual(
    values = c("Controls" = "#4CAF50", "Patients" = "#F44336")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p)

# 10) Salvataggi ------------------------------------------------------------
out_dir <- here::here("src", "wcst", "documentation", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
readr::write_csv(
  all_params,
  file.path(out_dir, "rw_dual_sticky_lapse_params.csv")
)
ggsave(
  file.path(out_dir, "rw_dual_sticky_lapse_boxplots.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

cat(
  "\nSaved:\n - ",
  file.path(out_dir, "rw_dual_sticky_lapse_params.csv"),
  "\n - ",
  file.path(out_dir, "rw_dual_sticky_lapse_boxplots.png"),
  "\n - ",
  stan_file,
  "\n\n",
  sep = ""
)
