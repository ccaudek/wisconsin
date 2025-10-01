
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

