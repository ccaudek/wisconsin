data {
  int<lower=0> N1; // number of observations for task 1
  int<lower=0> N2; // number of observations for task 2
  int<lower=0> K1; // number of predictors for task 1
  int<lower=0> K2; // number of predictors for task 2
  matrix[N1, K1] X1; // predictor matrix for task 1
  matrix[N2, K2] X2; // predictor matrix for task 2
  array[N1] int<lower=0, upper=1> y1; // binary outcome for task 1
  array[N2] int<lower=0, upper=1> y2; // binary outcome for task 2
}
parameters {
  vector[K1] beta1; // coefficients for task 1
  vector[K2] beta2; // coefficients for task 2
  real<lower=0, upper=1> auc1; // AUC for task 1
  real<lower=0, upper=1> auc2; // AUC for task 2
}
model {
  // Priors
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  auc1 ~ beta(2, 2);
  auc2 ~ beta(2, 2);
  
  // Likelihood
  y1 ~ bernoulli_logit(X1 * beta1);
  y2 ~ bernoulli_logit(X2 * beta2);
  
  // AUC computation (approximation using logistic regression coefficients)
  target += normal_lpdf(inv_Phi(auc1) | sqrt(sum(square(X1 * beta1)))
                                        / sqrt(8), 1);
  target += normal_lpdf(inv_Phi(auc2) | sqrt(sum(square(X2 * beta2)))
                                        / sqrt(8), 1);
}
generated quantities {
  real auc_diff = auc2 - auc1; // Difference in AUC (task2 - task1)
}
