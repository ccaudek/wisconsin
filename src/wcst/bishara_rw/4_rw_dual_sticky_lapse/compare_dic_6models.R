# ============================================================
# Confronto modelli (per-soggetto) via DIC — 6 modelli
# 1) full                : RW dual-alpha + sticky + lapse
# 2) no_kappa            : RW dual-alpha + (no sticky) + lapse
# 3) no_lapse            : RW dual-alpha + sticky + (no lapse)
# 4) basic_rw            : RW dual-alpha (no sticky, no lapse)  [storico]
# 5) rw_single_alpha_dic : RW single-alpha + beta                [nuovo]
# 6) bishara             : Attention model (Bishara et al.)      [nuovo]
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, tidyverse, cmdstanr, posterior, matrixStats, purrr)

set.seed(20251001)

# ------------------------------------------------------------
# Dati (pipeline esistente)
# ------------------------------------------------------------
stan_data <- readRDS(here::here("src", "wcst", "data", "wcst_stan_list.RDS"))

# ------------------------------------------------------------
# Percorso ai file .stan (spostati nella cartella `stan`)
# ------------------------------------------------------------
stan_dir <- here::here(
  "src",
  "wcst",
  "bishara_rw",
  "4_rw_dual_sticky_lapse",
  "stan"
)
req <- c(
  "rw_dual_sticky_lapse.stan",
  "rw_single_alpha_dic.stan",
  "bishara_single_subject_aligned.stan"
)
stopifnot(all(file.exists(file.path(stan_dir, req))))

# ------------------------------------------------------------
# Prep dati per soggetto — RW vs Bishara (diversa codifica resp_*)
# ------------------------------------------------------------
# RW models expect resp_* come indicatori 0/1 (come nei .stan RW)
prep_subject_data_rw <- function(i, use_kappa = 1L, use_lapse = 1L) {
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

# Bishara model expects resp_* come INDICI 1..4 del mazzo che matcha color/shape/number
# e rew/los in {-1,0,1} ma nel codice si usa solo (rew==1) per il feedback corretto.
prep_subject_data_bishara <- function(i) {
  Ti <- stan_data$Tsubj[i]
  list(
    T = as.integer(Ti),
    rew = as.integer(stan_data$rew[i, 1:Ti]), # valori originari (-1/0/1) o 0/1: ok
    los = as.integer(stan_data$los[i, 1:Ti]), # non utilizzato nel .stan, ma richiesto dal data block
    resp_choice = as.integer(stan_data$resp_choice[i, 1:Ti]),
    resp_color = as.integer(stan_data$resp_color[i, 1:Ti]), # *** INDICI 1..4 ***
    resp_shape = as.integer(stan_data$resp_shape[i, 1:Ti]), # *** INDICI 1..4 ***
    resp_number = as.integer(stan_data$resp_number[i, 1:Ti]) # *** INDICI 1..4 ***
  )
}

# ------------------------------------------------------------
# Verosimiglianze in R per Dhat
# ------------------------------------------------------------

# RW dual-alpha (+sticky/lapse opzionali) — parametri: alpha_pos, alpha_neg, beta, (kappa), (lapse)
loglik_rw_dual_R <- function(
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
    if (use_kappa == 1L && prev > 0L && !is.null(pars$kappa))
      logits[prev] <- logits[prev] + pars$kappa
    lp <- logits[choice[t]] - matrixStats::logSumExp(logits)
    if (use_lapse == 1L && !is.null(pars$lapse)) {
      ll <- ll + log((1 - pars$lapse) * exp(lp) + pars$lapse * 0.25)
    } else {
      ll <- ll + lp
    }
    a <- choice[t]
    r <- (rew01[t] == 1L) - (los01[t] == 1L) # +1/-1
    pe <- r - Q[a]
    if (pe >= 0) Q[a] <- Q[a] + pars$alpha_pos * pe else
      Q[a] <- Q[a] + pars$alpha_neg * pe
    prev <- a
  }
  ll
}

# RW single-alpha — parametri: alpha, beta
loglik_rw_singlealpha_R <- function(choice, rew01, los01, pars) {
  T <- length(choice)
  Q <- rep(0, 4L)
  ll <- 0
  for (t in seq_len(T)) {
    logits <- pars$beta * Q
    lp <- logits[choice[t]] - matrixStats::logSumExp(logits)
    ll <- ll + lp
    a <- choice[t]
    r <- (rew01[t] == 1L) - (los01[t] == 1L)
    Q[a] <- Q[a] + pars$alpha * (r - Q[a])
  }
  ll
}

# Bishara et al. — parametri: r (reward LR), p (punishment LR), d (decision consistency)
# Likelihood esatta come nel .stan: Eq.9 per scelta; update Eq.6/8; feedback signal Eq.13-14 con f=1.
loglik_bishara_R <- function(
  resp_choice,
  rew,
  resp_color,
  resp_shape,
  resp_number,
  pars
) {
  T <- length(resp_choice)
  a <- rep(1 / 3, 3) # attenzione iniziale [1/3,1/3,1/3]
  ll <- 0
  for (t in seq_len(T)) {
    if (resp_choice[t] > 0) {
      # match vectors m[k, 3]
      m <- matrix(0, nrow = 4, ncol = 3)
      if (resp_color[t] >= 1 && resp_color[t] <= 4) m[resp_color[t], 1] <- 1
      if (resp_shape[t] >= 1 && resp_shape[t] <= 4) m[resp_shape[t], 2] <- 1
      if (resp_number[t] >= 1 && resp_number[t] <= 4) m[resp_number[t], 3] <- 1

      a_pow <- a^pars$d
      pile_values <- as.vector(m %*% a_pow) # dot_product(a_pow, m[k])

      # log p(choice_t)
      denom <- matrixStats::logSumExp(log(pile_values))
      ll <- ll + (log(pile_values[resp_choice[t]]) - denom)

      # scelto quale dimensione?
      chosen_dim <- 0L
      if (resp_choice[t] == resp_color[t]) chosen_dim <- 1L else if (
        resp_choice[t] == resp_shape[t]
      )
        chosen_dim <- 2L else if (resp_choice[t] == resp_number[t])
        chosen_dim <- 3L

      if (chosen_dim > 0L) {
        feedback_correct <- ifelse(rew[t] == 1L, 1L, 0L)

        # dimensioni "consistenti" col feedback
        consistent_dims <- rep(0L, 3)
        if (feedback_correct == 1L) {
          # corretto: consistenti se matchano il mazzo scelto
          for (dim in 1:3)
            if (m[resp_choice[t], dim] == 1) consistent_dims[dim] <- 1L
        } else {
          # errato: consistenti se NON matchano il mazzo scelto
          for (dim in 1:3)
            if (m[resp_choice[t], dim] == 0) consistent_dims[dim] <- 1L
        }

        # segnale s con f=1 (proporzionale alle a correnti sulle dimensioni consistenti)
        denom_s <- sum((a[consistent_dims == 1])^1.0)
        s <- rep(0, 3)
        for (dim in 1:3) {
          if (consistent_dims[dim] == 1L) s[dim] <- (a[dim]^1.0) / denom_s
        }

        # update attenzione: corretto -> r ; errato -> p
        if (feedback_correct == 1L) {
          a <- (1 - pars$r) * a + pars$r * s
        } else {
          a <- (1 - pars$p) * a + pars$p * s
        }
      }
    }
  }
  ll
}

# ------------------------------------------------------------
# Calcolo DIC generico (accetta sia log_lik vettoriale che scalare)
# ------------------------------------------------------------
compute_dic_subject_generic <- function(fit, Dhat) {
  # Estrae log_lik per draw: matrice draws x K (K=T oppure 1)
  ll_mat <- fit$draws("log_lik", format = "matrix")
  Dbar <- -2 * mean(rowSums(ll_mat))
  Dhat_val <- -2 * Dhat
  pD <- Dbar - Dhat_val
  DIC <- Dhat_val + 2 * pD
  list(DIC = DIC, Dbar = Dbar, Dhat = Dhat_val, pD = pD)
}

# ------------------------------------------------------------
# Registry dei modelli
# ------------------------------------------------------------
registry <- tibble::tribble(
  ~model_id,
  ~stan_file,
  ~prep_fun,
  ~theta_vars,
  ~loglik_fun,
  ~flags,
  "full",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse"),
  "loglik_rw_dual_R",
  list(use_kappa = 1L, use_lapse = 1L),
  "no_kappa",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse"),
  "loglik_rw_dual_R",
  list(use_kappa = 0L, use_lapse = 1L),
  "no_lapse",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  c("alpha_pos", "alpha_neg", "beta", "kappa", "lapse"),
  "loglik_rw_dual_R",
  list(use_kappa = 1L, use_lapse = 0L),
  "basic_rw",
  "rw_dual_sticky_lapse.stan",
  prep_subject_data_rw,
  c("alpha_pos", "alpha_neg", "beta"),
  "loglik_rw_dual_R",
  list(use_kappa = 0L, use_lapse = 0L),
  "rw_single_alpha_dic",
  "rw_single_alpha_dic.stan",
  prep_subject_data_rw,
  c("alpha", "beta"),
  "loglik_rw_singlealpha_R",
  list(),
  "bishara",
  "bishara_single_subject_aligned.stan",
  prep_subject_data_bishara,
  c("r", "p", "d"),
  "loglik_bishara_R",
  list()
)

# Pre-compilazione
registry <- registry %>%
  mutate(
    stan_path = file.path(stan_dir, stan_file),
    mod_obj = purrr::map(stan_path, cmdstan_model)
  )

# ------------------------------------------------------------
# Fit per soggetto per un modello
# ------------------------------------------------------------
fit_subject_one <- function(
  row,
  i,
  chains = 2,
  iter_warmup = 600,
  iter_sampling = 600,
  parallel_chains = min(chains, max(1, parallel::detectCores() - 1))
) {
  # prepara dati
  data_i <- do.call(row$prep_fun[[1]], c(list(i = i), row$flags[[1]]))

  # sample
  fit <- row$mod_obj[[1]]$sample(
    data = data_i,
    seed = 20251001 + i,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 0
  )

  # posterior mean dei parametri richiesti
  summ <- fit$summary(variables = row$theta_vars[[1]])
  theta_bar <- setNames(as.list(summ$mean), summ$variable)

  # costruisci Dhat con la loglik R coerente
  ll_fun <- get(row$loglik_fun[[1]])

  Dhat_arglist <- switch(
    row$model_id,
    "full" = list(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar,
      use_kappa = 1L,
      use_lapse = 1L
    ),
    "no_kappa" = list(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar,
      use_kappa = 0L,
      use_lapse = 1L
    ),
    "no_lapse" = list(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar,
      use_kappa = 1L,
      use_lapse = 0L
    ),
    "basic_rw" = list(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar,
      use_kappa = 0L,
      use_lapse = 0L
    ),
    "rw_single_alpha_dic" = list(
      choice = data_i$choice,
      rew01 = data_i$rew,
      los01 = data_i$los,
      pars = theta_bar
    ),
    "bishara" = list(
      resp_choice = data_i$resp_choice,
      rew = data_i$rew,
      resp_color = data_i$resp_color,
      resp_shape = data_i$resp_shape,
      resp_number = data_i$resp_number,
      pars = theta_bar
    )
  )

  Dhat <- do.call(ll_fun, Dhat_arglist)
  dic <- compute_dic_subject_generic(fit, Dhat)

  # parametri (median) per traccia
  pars <- fit$summary(variables = row$theta_vars[[1]]) %>%
    select(variable, median) %>%
    pivot_wider(names_from = variable, values_from = median) %>%
    mutate(subject = i)

  list(dic = dic, pars = pars)
}

# ------------------------------------------------------------
# Runner su tutti i soggetti per un modello
# ------------------------------------------------------------
fit_model_all <- function(
  row,
  chains = 2,
  iter_warmup = 600,
  iter_sampling = 600
) {
  N <- stan_data$N
  message(
    "Model ",
    row$model_id,
    " (use_kappa=",
    row$flags[[1]]$use_kappa %||% NA_integer_,
    ", use_lapse=",
    row$flags[[1]]$use_lapse %||% NA_integer_,
    ")"
  )

  safe <- safely(fit_subject_one, otherwise = NULL, quiet = FALSE)

  res <- map(
    seq_len(N),
    ~ safe(
      row = row,
      i = .x,
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling
    )
  )

  fail <- which(map_lgl(res, ~ !is.null(.x$error)))
  if (length(fail)) {
    message("Falliti soggetti: ", paste(fail, collapse = ", "))
    message("Primo errore: ", res[[fail[1]]]$error$message)
  }
  ok <- compact(map(res, "result"))
  stopifnot(length(ok) > 0)

  dic_tbl <- map_dfr(seq_along(ok), function(j) {
    tibble(
      subject = ok[[j]]$pars$subject,
      DIC = ok[[j]]$dic$DIC,
      Dbar = ok[[j]]$dic$Dbar,
      Dhat = ok[[j]]$dic$Dhat,
      pD = ok[[j]]$dic$pD
    )
  }) %>%
    arrange(subject)

  pars_tbl <- bind_rows(map(ok, "pars")) %>%
    mutate(
      group = stan_data$group[subject],
      group_label = factor(
        group,
        levels = c(1, 2),
        labels = c("Controls", "Patients")
      )
    )

  list(
    model_id = row$model_id,
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

# ------------------------------------------------------------
# Esecuzione per tutti i modelli
# ------------------------------------------------------------
fits_dic <- purrr::map(
  seq_len(nrow(registry)),
  ~ fit_model_all(
    registry %>% dplyr::slice(.x),
    chains = 2,
    iter_warmup = 600,
    iter_sampling = 600
  )
)

# Tabella confronto finale
dic_tot <- bind_rows(lapply(fits_dic, `[[`, "dic_total")) %>%
  mutate(model_id = registry$model_id) %>%
  select(model_id, everything()) %>%
  arrange(DIC_total)

print(dic_tot, n = Inf)

# Delta DIC rispetto al migliore
best <- min(dic_tot$DIC_total)
dic_tot %>%
  mutate(delta_DIC = DIC_total - best) %>%
  arrange(delta_DIC) %>%
  print(n = Inf)

# Salvataggi
out_dir <- here::here("src", "wcst", "documentation", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(
  list(registry = registry, fits = fits_dic, dic_table = dic_tot),
  file.path(out_dir, "subjectwise_model_comparison_dic.rds")
)
readr::write_csv(
  dic_tot,
  file.path(out_dir, "subjectwise_model_comparison_dic.csv")
)


# --- PATCH APPLICATO IL 2025-10-01 --------------------------------------

# ===== PATCH: robust iteration over registry rows =====
# Instead of group_split(row_number()) which returns a list of 1-col tibbles
# and may pass atomic vectors, iterate by explicit row index and pass a 1-row tibble.

fits_dic <- purrr::map(
  seq_len(nrow(registry)),
  ~ fit_model_all(
    registry %>% dplyr::slice(.x),
    chains = 2,
    iter_warmup = 600,
    iter_sampling = 600
  )
)
