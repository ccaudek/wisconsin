# ============================================================
# RW dual-α + sticky + lapse — modello gerarchico unico
# - legge i dati empirici con le tue funzioni
# - costruisce e compila il modello gerarchico
# - MCMC, estrazione parametri per soggetto e iper-parametri
# - LOO/ELPD complessivo per confronto modelli
# ============================================================

# 0) Packages ---------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here,
  tidyverse,
  cmdstanr,
  posterior,
  bayesplot,
  loo
)

set.seed(20251001)

# 1) Sorgenti progetto e dati ----------------------------------------------
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

generate_csv("controls")
generate_csv("patients")
stan_data_raw <- process_and_prepare_stan_data()
str(stan_data_raw)

# 2) Costruisci data list per il modello gerarchico -------------------------
to_int_mat <- function(M) {
  M2 <- as.matrix(M)
  storage.mode(M2) <- "integer"
  unname(M2) # rimuove dimnames, ma mantiene dim [N,T]
}

stan_data_hier <- list(
  N = as.integer(stan_data_raw$N),
  T = as.integer(stan_data_raw$T),
  Tsubj = as.integer(stan_data_raw$Tsubj),
  choice = to_int_mat(stan_data_raw$resp_choice), # 1..4
  rew = to_int_mat(stan_data_raw$rew), # 0/1 già ok
  los = to_int_mat(stan_data_raw$los) # -1/0 (il modello usa check == -1)
)

# sanity-check veloci:
stopifnot(all(stan_data_hier$choice >= 1 & stan_data_hier$choice <= 4))
stopifnot(all(stan_data_hier$rew %in% c(0L, 1L)))
stopifnot(all(stan_data_hier$los %in% c(-1L, 0L)))

# 3) Scrivi/compila il modello ---------------------------------------------
stan_dir <- here::here("src", "wcst", "documentation", "stan")
if (!dir.exists(stan_dir)) dir.create(stan_dir, recursive = TRUE)
stan_file <- file.path(stan_dir, "11_rw_dual_sticky_lapse_hier.stan")

# Se non è già presente sul disco, scrivilo (incolla il codice Stan qui sopra se serve)
# writeLines(stan_code_hier, con = stan_file) # <- opzionale se vuoi scriverlo da R

mod <- cmdstan_model(stan_file)

# 4) Sampling ---------------------------------------------------------------
chains <- 4
parallel_chains <- min(chains, max(1, parallel::detectCores() - 1))
iter_warmup <- 1000
iter_sampling <- 1000
adapt_delta <- 0.9
max_treedepth <- 12

fit <- mod$sample(
  data = stan_data_hier,
  seed = 20251001,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  refresh = 100,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth
)

print(fit$summary(c(
  "mu_ap",
  "mu_an",
  "mu_b",
  "mu_k",
  "mu_l",
  "sigma_ap",
  "sigma_an",
  "sigma_b",
  "sigma_k",
  "sigma_l"
)))

# 5) Estrai parametri per soggetto ------------------------------------------
# alpha_pos[N], alpha_neg[N], beta[N], kappa[N], lapse[N]
S <- fit$summary(
  variables = c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")
)

rx <- "^(alpha_pos|alpha_neg|beta|kappa|lapse)\\[(\\d+)\\]$"
M <- str_match(S$variable, rx)
# M[,2] = param, M[,3] = subject come stringa

params_subj <- S %>%
  mutate(
    param = M[, 2],
    subject = suppressWarnings(as.integer(M[, 3]))
  ) %>%
  filter(!is.na(subject), !is.na(param)) %>%
  # se (param, subject) appare più di una volta, fai il collapse:
  group_by(subject, param) %>%
  summarise(median = median(median), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = param, values_from = median) %>%
  arrange(subject) %>%
  mutate(
    group = stan_data_raw$group[subject],
    group_label = factor(
      group,
      levels = c(1, 2),
      labels = c("Controls", "Patients")
    )
  )

# controllo che adesso siano tutte colonne numeriche base
stopifnot(all(sapply(
  params_subj[, c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")],
  is.numeric
)))
stopifnot(all(is.finite(as.matrix(params_subj[, c(
  "alpha_pos",
  "alpha_neg",
  "beta",
  "kappa",
  "lapse"
)]))))

# 6) Plot per parametro ------------------------------------------------------
params_long <- params_subj %>%
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
    title = "RW dual-α + sticky + lapse — gerarchico (posterior medians)",
    x = NULL,
    y = "Estimate"
  ) +
  scale_fill_manual(
    values = c("Controls" = "#4CAF50", "Patients" = "#F44336")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p)

# 7) LOO/ELPD complessivo ---------------------------------------------------
# Estrai i log_lik per osservazione: draws x (N*T)
ll_draws <- fit$draws("log_lik", format = "matrix") # S_draws x (N*T)
N <- stan_data_hier$N
T <- stan_data_hier$T
stopifnot(ncol(ll_draws) == N * T)

# Crea maschera delle osservazioni reali (TRUE se t <= Tsubj[s])
mask <- rep(FALSE, N * T)
idx <- 1
for (s in 1:N) {
  for (t in 1:T) {
    mask[idx] <- (t <= stan_data_hier$Tsubj[s])
    idx <- idx + 1
  }
}

ll_obs <- ll_draws[, mask, drop = FALSE] # solo colonne osservate
loo_fit <- loo::loo(ll_obs, r_eff = NA)

print(loo_fit)

# 8) Salvataggi -------------------------------------------------------------
out_dir <- here::here("src", "wcst", "documentation", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

readr::write_csv(params_subj, file.path(out_dir, "hier_params_per_subject.csv"))
ggsave(
  file.path(out_dir, "hier_params_boxplots.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300
)

saveRDS(loo_fit, file.path(out_dir, "hier_model_loo.rds"))
cat(
  "\nSaved:\n - hier_params_per_subject.csv\n - hier_params_boxplots.png\n - hier_model_loo.rds\n"
)
