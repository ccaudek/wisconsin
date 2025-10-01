data {
  int<lower=1> N;                        // soggetti
  int<lower=1> T;                        // trials max
  array[N] int<lower=1> Tsubj;           // trials osservati per soggetto
  array[N, T] int<lower=1, upper=4> choice; // scelta 1..4
  array[N, T] int<lower=0, upper=1> rew;    // 1 reward, 0 altrimenti
  array[N, T] int<lower=-1, upper=0> los;   // -1 loss, 0 altrimenti
}

parameters {
  // Non-centered: param_raw[s] = mu + sigma * z[s]
  vector[N] z_ap;     // alpha_pos_raw
  vector[N] z_an;     // alpha_neg_raw
  vector[N] z_b;      // beta_raw
  vector[N] z_k;      // kappa
  vector[N] z_l;      // lapse_raw

  real mu_ap;
  real mu_an;
  real mu_b;
  real mu_k;
  real mu_l;

  real<lower=0> sigma_ap;
  real<lower=0> sigma_an;
  real<lower=0> sigma_b;
  real<lower=0> sigma_k;
  real<lower=0> sigma_l;
}

transformed parameters {
  // Parametri per soggetto (s) su scale utili
  vector[N] alpha_pos_raw = mu_ap + sigma_ap * z_ap;
  vector[N] alpha_neg_raw = mu_an + sigma_an * z_an;
  vector[N] beta_raw      = mu_b  + sigma_b  * z_b;
  vector[N] kappa         = mu_k  + sigma_k  * z_k;
  vector[N] lapse_raw     = mu_l  + sigma_l  * z_l;

  vector<lower=0, upper=1>[N] alpha_pos = inv_logit(alpha_pos_raw);
  vector<lower=0, upper=1>[N] alpha_neg = inv_logit(alpha_neg_raw);
  vector<lower=0>[N]          beta      = exp(beta_raw);
  vector<lower=0, upper=1>[N] lapse     = inv_logit(lapse_raw);
}

model {
  // Hyper-priors molto deboli (pooling minimo)
  mu_ap ~ normal(0, 2.5);
  mu_an ~ normal(0, 2.5);
  mu_b  ~ normal(0, 2.5);
  mu_k  ~ normal(0, 2.5);
  mu_l  ~ normal(-2.2, 2.5);     // media ~0.10 su scala naturale

  sigma_ap ~ student_t(3, 0, 5); // half-Student-t(3,5)
  sigma_an ~ student_t(3, 0, 5);
  sigma_b  ~ student_t(3, 0, 5);
  sigma_k  ~ student_t(3, 0, 5);
  sigma_l  ~ student_t(3, 0, 5);

  z_ap ~ std_normal();
  z_an ~ std_normal();
  z_b  ~ std_normal();
  z_k  ~ std_normal();
  z_l  ~ std_normal();

  // Likelihood
  for (s in 1:N) {
    vector[4] Q = rep_vector(0.0, 4);
    int prev = 0;

    for (t in 1:Tsubj[s]) {
      vector[4] logits = beta[s] * Q;
      if (prev > 0) logits[prev] += kappa[s];

      // mixture con lapse (uniforme su 4 azioni)
      target += log_mix(
        1.0 - lapse[s],
        categorical_logit_lpmf(choice[s, t] | logits),
        -log(4.0)
      );

      {
        int a = choice[s, t];
        real r = (rew[s, t] == 1 ? 1.0 : 0.0) - (los[s, t] == -1 ? 1.0 : 0.0);
        real pe = r - Q[a];
        if (pe >= 0) Q[a] += alpha_pos[s] * pe;
        else         Q[a] += alpha_neg[s] * pe;
        prev = a;
      }
    }
  }
}

generated quantities {
  // log-lik per osservazione (N x T), con zeri fuori dai trials osservati
  array[N, T] real log_lik;
  {
    for (s in 1:N) {
      vector[4] Q = rep_vector(0.0, 4);
      int prev = 0;
      for (t in 1:T) {
        if (t <= Tsubj[s]) {
          vector[4] logits = exp(beta_raw[s]) * Q;
          if (prev > 0) logits[prev] += kappa[s];
          log_lik[s, t] = log_mix(
            1.0 - inv_logit(lapse_raw[s]),
            categorical_logit_lpmf(choice[s, t] | logits),
            -log(4.0)
          );
          {
            int a = choice[s, t];
            real r = (rew[s, t] == 1 ? 1.0 : 0.0) - (los[s, t] == -1 ? 1.0 : 0.0);
            real pe = r - Q[a];
            if (pe >= 0) Q[a] += inv_logit(alpha_pos_raw[s]) * pe;
            else         Q[a] += inv_logit(alpha_neg_raw[s]) * pe;
            prev = a;
          }
        } else {
          log_lik[s, t] = 0; // verrà filtrato in R
        }
      }
    }
  }
}
