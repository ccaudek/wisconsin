data {
  int<lower=1> T;                          // trials
  array[T] int<lower=1, upper=4> choice;   // 1..4
  array[T] int<lower=0, upper=1> rew;      // 1 reward, 0 altrimenti
  array[T] int<lower=0, upper=1> los;      // 1 loss, 0 altrimenti
  // compatibilità con la tua pipeline (non usati in likelihood)
  array[T] int<lower=0, upper=1> resp_color;
  array[T] int<lower=0, upper=1> resp_shape;
  array[T] int<lower=0, upper=1> resp_number;
}
parameters {
  real alpha_raw;     // -> (0,1) via inv_logit
  real beta_raw;      // -> (0,+inf) via exp
}
transformed parameters {
  real<lower=0, upper=1> alpha = inv_logit(alpha_raw);
  real<lower=0>          beta  = exp(beta_raw);
}
model {
  // priors deboli
  alpha_raw ~ normal(0, 1);
  beta_raw  ~ normal(0, 0.7);

  {
    vector[4] Q = rep_vector(0.0, 4);
    int prev = 0;
    for (t in 1:T) {
      vector[4] logits = beta * Q;
      target += categorical_logit_lpmf(choice[t] | logits);

      // update RW single-alpha (positivo o negativo sempre α)
      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
      real pe = r - Q[a];
      Q[a] += alpha * pe;
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
      log_lik[t] = categorical_logit_lpmf(choice[t] | logits);

      int a = choice[t];
      real r = (rew[t] == 1 ? 1.0 : 0.0) - (los[t] == 1 ? 1.0 : 0.0);
      real pe = r - Q[a];
      Q[a] += inv_logit(alpha_raw) * pe;
      prev = a;
    }
  }
}
