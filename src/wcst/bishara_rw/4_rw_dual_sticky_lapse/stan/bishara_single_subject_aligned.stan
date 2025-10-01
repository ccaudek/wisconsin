
// bishara_single_subject.stan (patched, modern array syntax)
// Attention-based WCST model (Bishara et al.)
// Emits per-trial log_lik_t and total log_lik for WAIC/LOO/DIC.

data {
  int<lower=1> T;
  array[T] int<lower=1, upper=4> resp_choice;  // chosen pile 1..4
  array[T] int<lower=1, upper=4> resp_color;   // pile matching color at t
  array[T] int<lower=1, upper=4> resp_shape;   // pile matching shape at t
  array[T] int<lower=1, upper=4> resp_number;  // pile matching number at t
  array[T] int rew;                            // feedback: 1=correct, else!=1
}

parameters {
  real<lower=0, upper=1> r;   // reward learning rate
  real<lower=0, upper=1> p;   // punishment learning rate
  real<lower=0> d;            // decision consistency
}

model {
  // Priors (tune to match your original spec if needed)
  r ~ beta(2, 2);
  p ~ beta(2, 2);
  d ~ gamma(2, 1);

  // Likelihood
  {
    vector[3] a = rep_vector(1.0/3.0, 3);   // initial attention
    for (t in 1:T) {
      // Match matrix m (4 x 3)
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      // Decision rule: values per pile = m * a^d
      vector[3] a_pow = pow(a, d);
      vector[4] pile_values = m * a_pow;

      // Normalize to probabilities
      vector[4] p_choice = pile_values / sum(pile_values);

      // Trial log-likelihood
      target += categorical_lpmf(resp_choice[t] | p_choice);

      // Feedback-consistent dimensions (f = 1)
      int correct = (rew[t] == 1);
      vector[3] s = rep_vector(0, 3);
      vector[3] mask;
      for (j in 1:3) {
        mask[j] = correct ? m[resp_choice[t], j] : 1 - m[resp_choice[t], j];
      }
      real denom_s = sum( a .* mask );
      if (denom_s > 0) s = (a .* mask) / denom_s;

      // Attention update
      if (correct == 1) a = (1 - r) * a + r * s;
      else              a = (1 - p) * a + p * s;
    }
  }
}

generated quantities {
  vector[T] log_lik_t;
  real log_lik;

  {
    vector[3] a = rep_vector(1.0/3.0, 3);
    log_lik = 0;
    for (t in 1:T) {
      matrix[4,3] m = rep_matrix(0, 4, 3);
      m[resp_color[t],  1] = 1;
      m[resp_shape[t],  2] = 1;
      m[resp_number[t], 3] = 1;

      vector[3] a_pow = pow(a, d);
      vector[4] pile_values = m * a_pow;
      vector[4] p_choice = pile_values / sum(pile_values);

      log_lik_t[t] = categorical_lpmf(resp_choice[t] | p_choice);
      log_lik += log_lik_t[t];

      int correct = (rew[t] == 1);
      vector[3] s = rep_vector(0, 3);
      vector[3] mask;
      for (j in 1:3) {
        mask[j] = correct ? m[resp_choice[t], j] : 1 - m[resp_choice[t], j];
      }
      real denom_s = sum( a .* mask );
      if (denom_s > 0) s = (a .* mask) / denom_s;
      if (correct == 1) a = (1 - r) * a + r * s;
      else              a = (1 - p) * a + p * s;
    }
  }
}
