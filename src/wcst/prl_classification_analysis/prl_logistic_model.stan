data {
  int<lower=0> N; // number of observations
  int<lower=0> P; // number of predictors
  matrix[N, P] X; // predictor matrix
  array[N] int<lower=0, upper=1> y; // binary outcome
}
parameters {
  vector[P] beta; // regression coefficients
  real alpha; // intercept
}
model {
  // Priors
  beta ~ normal(0, 5);
  alpha ~ normal(0, 5);
  
  // Likelihood
  y ~ bernoulli_logit(alpha + X * beta);
}
generated quantities {
  vector[N] y_pred;
  
  for (n in 1 : N) {
    y_pred[n] = bernoulli_logit_rng(alpha + dot_product(X[n], beta));
  }
}
