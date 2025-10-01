# ============================================================
# Confronto modelli RW per-soggetto via DIC
# - Fit MCMC per ogni soggetto e modello (full, no_kappa, no_lapse, basic)
# - DIC_i per soggetto; somma DIC_totale per modello
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, tidyverse, cmdstanr, posterior, stringr, matrixStats)

set.seed(20251001)

# --- Dati: tua pipeline ----------------------------------------------------

stan_data <- readRDS(
  here::here(
    "src",
    "wcst",
    "data",
    "wcst_stan_list.RDS"
  )
)

# --- Stan model ------------------------------------------------------------
stan_dir <- here::here("src", "wcst", "bishara_rw", "4_rw_dual_sticky_lapse")
stan_file <- file.path(stan_dir, "rw_dual_sticky_lapse.stan")
stopifnot(file.exists(stan_file))
mod <- cmdstan_model(stan_file)

# --- Helper: prepara dati soggetto + flag modello -------------------------
prep_subject_data <- function(i, use_kappa = 1L, use_lapse = 1L) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    rew = as.integer(stan_data$rew[i, 1:Ti] == 1L),
    los = as.integer(stan_data$los[i, 1:Ti] == -1L),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti] == 1L),
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti] == 1L),
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti] == 1L),
    use_kappa = as.integer(use_kappa),
    use_lapse = as.integer(use_lapse)
  )
}

# --- Verosimiglianza in R per D(θbar) -------------------------------------
loglik_rw_R <- function(
  choice,
  rew01,
  los01,
  pars,
  use_kappa = 1L,
  use_lapse = 1L
) {
  T <- length(choice)
  Q <- rep(0, 4L)
  prev <- 0L
  ll <- 0
  for (t in seq_len(T)) {
    logits <- pars$beta * Q
    if (use_kappa == 1L && prev > 0L) logits[prev] <- logits[prev] + pars$kappa
    # log-softmax: lp(choice) = logit(choice) - log(sum exp(logits))
    lp_model <- logits[choice[t]] - matrixStats::logSumExp(logits)
    if (use_lapse == 1L) {
      ll <- ll + log((1 - pars$lapse) * exp(lp_model) + pars$lapse * 0.25)
    } else {
      ll <- ll + lp_model
    }
    a <- choice[t]
    r <- (rew01[t] == 1L) - (los01[t] == 1L) # in {-1,+1}
    pe <- r - Q[a]
    if (pe >= 0) Q[a] <- Q[a] + pars$alpha_pos * pe else
      Q[a] <- Q[a] + pars$alpha_neg * pe
    prev <- a
  }
  ll # log-likelihood (costanti incluse)
}

# --- DIC per soggetto a partire dal fit MCMC -------------------------------
compute_dic_subject <- function(fit, data_i, use_kappa, use_lapse) {
  # Dbar = -2 * mean_draws( sum_t log_lik[t] )
  ll_mat <- fit$draws("log_lik", format = "matrix") # draws x T
  Dbar <- -2 * mean(rowSums(ll_mat))

  # θbar = posterior mean su scala naturale
  summ <- fit$summary(c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse"))
  theta_bar <- list(
    alpha_pos = summ$mean[summ$variable == "alpha_pos"],
    alpha_neg = summ$mean[summ$variable == "alpha_neg"],
    beta = summ$mean[summ$variable == "beta"],
    kappa = summ$mean[summ$variable == "kappa"],
    lapse = summ$mean[summ$variable == "lapse"]
  )

  # Dhat = -2 * loglik( y | θbar )
  Dhat <- -2 *
    loglik_rw_R(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar,
      use_kappa = use_kappa,
      use_lapse = use_lapse
    )

  pD <- Dbar - Dhat
  DIC <- Dhat + 2 * pD # == 2*Dbar - Dhat
  list(DIC = DIC, Dbar = Dbar, Dhat = Dhat, pD = pD, theta_bar = theta_bar)
}

# --- Fit MCMC per soggetto (ritorna DIC e parametri) -----------------------
fit_subject_mcmc_dic <- function(
  i,
  use_kappa,
  use_lapse,
  chains = 2,
  iter_warmup = 800,
  iter_sampling = 800,
  parallel_chains = min(chains, max(1, parallel::detectCores() - 1))
) {
  data_i <- prep_subject_data(i, use_kappa, use_lapse)
  fit <- mod$sample(
    data = data_i,
    seed = 20251001 + i,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )
  dic <- compute_dic_subject(fit, data_i, use_kappa, use_lapse)
  pars <- fit$summary(c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse")) %>%
    select(variable, median) %>%
    pivot_wider(names_from = variable, values_from = median) %>%
    mutate(subject = i)
  list(dic = dic, pars = pars)
}

# --- Runner per un modello su tutti i soggetti -----------------------------
fit_model_dic_all <- function(
  model_id,
  use_kappa,
  use_lapse,
  chains = 2,
  iter_warmup = 800,
  iter_sampling = 800
) {
  N <- stan_data$N
  message(
    "Model ",
    model_id,
    "  (use_kappa=",
    use_kappa,
    ", use_lapse=",
    use_lapse,
    ")"
  )
  safe <- purrr::safely(fit_subject_mcmc_dic, otherwise = NULL, quiet = FALSE)
  res <- purrr::map(
    seq_len(N),
    ~ safe(
      .x,
      use_kappa,
      use_lapse,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling
    )
  )
  fail <- which(purrr::map_lgl(res, ~ !is.null(.x$error)))
  if (length(fail)) {
    message("Falliti soggetti: ", paste(fail, collapse = ", "))
    message("Primo errore: ", res[[fail[1]]]$error$message)
  }
  ok <- purrr::compact(purrr::map(res, "result"))
  stopifnot(length(ok) > 0)

  dic_tbl <- purrr::map_dfr(seq_along(ok), function(j) {
    tibble(
      subject = ok[[j]]$pars$subject,
      DIC = ok[[j]]$dic$DIC,
      Dbar = ok[[j]]$dic$Dbar,
      Dhat = ok[[j]]$dic$Dhat,
      pD = ok[[j]]$dic$pD
    )
  }) %>%
    arrange(subject)

  pars_tbl <- bind_rows(purrr::map(ok, "pars")) %>%
    mutate(
      group = stan_data$group[subject],
      group_label = factor(
        group,
        levels = c(1, 2),
        labels = c("Controls", "Patients")
      )
    )

  list(
    model_id = model_id,
    use_kappa = use_kappa,
    use_lapse = use_lapse,
    dic_by_subject = dic_tbl,
    dic_total = dic_tbl %>%
      summarise(
        DIC_total = sum(DIC),
        Dbar_total = sum(Dbar),
        Dhat_total = sum(Dhat),
        pD_total = sum(pD)
      ),
    params = pars_tbl
  )
}

# --- Modelli da confrontare ------------------------------------------------
models <- tribble(
  ~model_id,
  ~use_kappa,
  ~use_lapse,
  "full",
  1L,
  1L,
  "no_kappa",
  0L,
  1L,
  "no_lapse",
  1L,
  0L,
  "basic_rw",
  0L,
  0L
)

# --- Esegui (riduci iter_* per run veloce) --------------------------------
fits_dic <- purrr::pmap(
  models,
  ~ fit_model_dic_all(
    ..1,
    ..2,
    ..3,
    chains = 2,
    iter_warmup = 600,
    iter_sampling = 600
  )
)

# --- Tabella confronto finale ----------------------------------------------
dic_tot <- bind_rows(lapply(fits_dic, `[[`, "dic_total")) %>%
  mutate(model_id = models$model_id) %>%
  select(model_id, everything()) %>%
  arrange(DIC_total)

print(dic_tot, n = Inf)

# (facoltativo) ΔDIC rispetto al migliore
best <- min(dic_tot$DIC_total)
dic_tot %>%
  mutate(delta_DIC = DIC_total - best) %>%
  arrange(delta_DIC) %>%
  print(n = Inf)

# --- Salvataggi -------------------------------------------------------------
out_dir <- here::here("src", "wcst", "documentation", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(
  list(models = models, fits = fits_dic, dic_table = dic_tot),
  file.path(out_dir, "subjectwise_model_comparison_dic.rds")
)
readr::write_csv(
  dic_tot,
  file.path(out_dir, "subjectwise_model_comparison_dic.csv")
)
