data {
  // ---- struttura dati task ----
  int<lower=1> N;                        // soggetti
  int<lower=1> T;                        // max trials
  array[N] int<lower=1, upper=T> Tsubj;  // n trials per soggetto
  array[N] int<lower=0, upper=1> holdout;

  array[N, T] real rew;                  // rinforzo (+)
  array[N, T] real los;                  // punizione (tipicamente <= 0 o 0/1 negativo)
  array[N, T] int rule_choice;           // (non usato qui ma lasciato per compat)
  array[N, T] int resp_choice;           // risposta effettiva (codifica che mappa su colore/forma/numero)
  array[N, T] int resp_color;            // codice opzione se scelta colore
  array[N, T] int resp_shape;            // codice opzione se scelta forma
  array[N, T] int resp_number;           // codice opzione se scelta numero

  // ---- iper-parametri "empirical Bayes" fissati (da Step-1) ----
  // alpha in (0,1) -> Beta(a_alpha, b_alpha)
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;

  // beta, inertia > 0 -> LogNormal(mu, sigma) su scala log
  real mu_log_beta;
  real<lower=1e-9> sigma_log_beta;

  real mu_log_inertia;
  real<lower=1e-9> sigma_log_inertia;
}

parameters {
  // Parametri per-soggetto (parzialmente "poolati" solo tramite i prior EB fissi)
  vector<lower=0, upper=1>[N] alpha;   // learning rate
  vector<lower=0>[N]          beta;    // inverse temperature
  vector<lower=0>[N]          inertia; // inerzia/perseverazione
}

model {
  // ---- Priors (Empirical Bayes): fissati a iper-parametri passati in 'data' ----
  alpha   ~ beta(a_alpha, b_alpha);
  beta    ~ lognormal(mu_log_beta,    sigma_log_beta);
  inertia ~ lognormal(mu_log_inertia, sigma_log_inertia);

  // ---- Likelihood ----
  for (i in 1:N) {
    if (holdout[i] == 0) {
      vector[3] Q = rep_vector(0.0, 3);  // valori RW per 3 dimensioni (colore/forma/numero)
      int last_dim = 0;

      for (t in 1:Tsubj[i]) {
        if (resp_choice[i, t] > 0) {
          // utilità (logits) per le 3 dimensioni
          vector[3] dim_logits = beta[i] * Q;

          // inerzia: bonus alla dimensione scelta nello step precedente
          if (last_dim > 0) {
            dim_logits[last_dim] += inertia[i];
          }

          vector[3] dim_probs = softmax(dim_logits);

          // ricava quale dimensione è stata effettivamente scelta in questo trial
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) {
            chosen_dim = 1;
          } else if (resp_choice[i, t] == resp_shape[i, t]) {
            chosen_dim = 2;
          } else if (resp_choice[i, t] == resp_number[i, t]) {
            chosen_dim = 3;
          }

          if (chosen_dim > 0) {
            // log-likelihood della scelta
            target += log(dim_probs[chosen_dim]);

            // outcome = reward + loss (assunto già scalato in [0,1] o compatibile con RW)
            real outcome = rew[i, t] + los[i, t];

            // aggiornamento RW per la dimensione scelta
            Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);

            // traccia della dimensione scelta per l'inerzia al prossimo trial
            last_dim = chosen_dim;
          }
        }
      } // end for t
    }   // end if holdout
  }     // end for i
}

generated quantities {
  // log-likelihood per-soggetto (utile per LOO)
  array[N] real log_lik;

  for (i in 1:N) {
    log_lik[i] = 0;
    vector[3] Q = rep_vector(0.0, 3);
    int last_dim = 0;

    for (t in 1:Tsubj[i]) {
      if (resp_choice[i, t] > 0) {
        vector[3] dim_logits = beta[i] * Q;
        if (last_dim > 0) {
          dim_logits[last_dim] += inertia[i];
        }
        vector[3] dim_probs = softmax(dim_logits);

        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) {
          chosen_dim = 1;
        } else if (resp_choice[i, t] == resp_shape[i, t]) {
          chosen_dim = 2;
        } else if (resp_choice[i, t] == resp_number[i, t]) {
          chosen_dim = 3;
        }

        if (chosen_dim > 0) {
          log_lik[i] += log(dim_probs[chosen_dim]);
          real outcome = rew[i, t] + los[i, t];
          Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
          last_dim = chosen_dim;
        }
      }
    }
  }
}
