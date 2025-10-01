data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> choice;
  array[T] int<lower=0, upper=1> rew;
  array[T] int<lower=0, upper=1> los;     // 1 se loss (da los==-1), 0 altrimenti
  array[T] int<lower=0, upper=1> resp_color;
  array[T] int<lower=0, upper=1> resp_shape;
  array[T] int<lower=0, upper=1> resp_number;

  int<lower=0, upper=1> use_kappa;        // 1 = usa sticky
  int<lower=0, upper=1> use_lapse;        // 1 = usa lapse
}
parameters {
  real alpha_pos_raw;
  real alpha_neg_raw;
  real beta_raw;
  real kappa;
  real lapse_raw;
}
transformed parameters {
  real<lower=0, upper=1> alpha_pos = inv_logit(alpha_pos_raw);
  real<lower=0, upper=1> alpha_neg = inv_logit(alpha_neg_raw);
  real<lower=0>          beta      = exp(beta_raw);
  real<lower=0, upper=1> lapse     = inv_logit(lapse_raw);
}
model {
  alpha_pos_raw ~ normal(0, 1);
  alpha_neg_raw ~ normal(0, 1);
  beta_raw      ~ normal(0.0, 0.7);
  kappa         ~ student_t(3, 0, 1);
  lapse_raw     ~ normal(-2.2, 1.0);

  {
    vector[4] Q = rep_vector(0.0, 4);
    int prev = 0;
    for (t in 1:T) {
      vector[4] logits = beta * Q;
      if (use_kappa == 1 && prev > 0) logits[prev] += kappa;

      real lp_model = categorical_logit_lpmf(choice[t] | logits);
      if (use_lapse == 1)
        target += log_mix(1.0 - lapse, lp_model, -log(4.0));
      else
        target += lp_model;

      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
      real pe = r - Q[a];
      if (pe >= 0) Q[a] += alpha_pos * pe;
      else         Q[a] += alpha_neg * pe;
      prev = a;
    }
  }
}
generated quantities {
  array[T] real log_lik;
  {
    vector[4] Q = rep_vector(0.0, 4);
    int prev = 0;
    for (t in 1:T) {
      vector[4] logits = exp(beta_raw) * Q;
      if (use_kappa == 1 && prev > 0) logits[prev] += kappa;

      real lp_model = categorical_logit_lpmf(choice[t] | logits);
      if (use_lapse == 1)
        log_lik[t] = log_mix(1.0 - inv_logit(lapse_raw), lp_model, -log(4.0));
      else
        log_lik[t] = lp_model;

      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
      real pe = r - Q[a];
      if (pe >= 0) Q[a] += inv_logit(alpha_pos_raw) * pe;
      else         Q[a] += inv_logit(alpha_neg_raw) * pe;
      prev = a;
    }
  }
}
