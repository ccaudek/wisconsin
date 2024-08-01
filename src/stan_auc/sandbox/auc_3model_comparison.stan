functions {
  real calculate_auc(vector y_pred, array[] int y) {
    int N = num_elements(y);
    int n_pos = sum(y);
    int n_neg = N - n_pos;
    vector[n_pos] pos_scores;
    vector[n_neg] neg_scores;
    int pos_index = 1;
    int neg_index = 1;
    
    for (i in 1 : N) {
      if (y[i] == 1) {
        pos_scores[pos_index] = y_pred[i];
        pos_index += 1;
      } else {
        neg_scores[neg_index] = y_pred[i];
        neg_index += 1;
      }
    }
    
    int concordant = 0;
    for (i in 1 : n_pos) {
      for (j in 1 : n_neg) {
        concordant += (pos_scores[i] > neg_scores[j]);
      }
    }
    
    return concordant * 1.0 / (n_pos * n_neg);
  }
}
data {
  int<lower=0> N1; // number of observations for model 1
  int<lower=0> N2; // number of observations for model 2
  int<lower=0> N3; // number of observations for model 3
  int<lower=0> K1; // number of predictors for model 1
  int<lower=0> K2; // number of predictors for model 2
  int<lower=0> K3; // number of predictors for model 3
  matrix[N1, K1] X1; // predictor matrix for model 1
  matrix[N2, K2] X2; // predictor matrix for model 2
  matrix[N3, K3] X3; // predictor matrix for model 3
  array[N1] int<lower=0, upper=1> y1; // binary outcome for model 1
  array[N2] int<lower=0, upper=1> y2; // binary outcome for model 2
  array[N3] int<lower=0, upper=1> y3; // binary outcome for model 3
}
parameters {
  vector[K1] beta1; // coefficients for model 1
  vector[K2] beta2; // coefficients for model 2
  vector[K3] beta3; // coefficients for model 3
}
model {
  // Priors
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  beta3 ~ normal(0, 5);
  
  // Likelihood
  y1 ~ bernoulli_logit(X1 * beta1);
  y2 ~ bernoulli_logit(X2 * beta2);
  y3 ~ bernoulli_logit(X3 * beta3);
}
generated quantities {
  real auc1;
  real auc2;
  real auc3;
  real auc_diff1;
  real auc_diff2;
  
  vector[N1] y_pred1 = inv_logit(X1 * beta1);
  vector[N2] y_pred2 = inv_logit(X2 * beta2);
  vector[N3] y_pred3 = inv_logit(X3 * beta3);
  
  auc1 = calculate_auc(y_pred1, y1);
  auc2 = calculate_auc(y_pred2, y2);
  auc3 = calculate_auc(y_pred3, y3);
  
  auc_diff1 = auc2 - auc1; // Difference in AUC (model2 - model1)
  auc_diff2 = auc1 - auc3; // Difference in AUC (model2 - model3)
}
