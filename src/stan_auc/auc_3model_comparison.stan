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
  int<lower=0> N4; // number of observations for model 4
  int<lower=0> N5; // number of observations for model 5
  int<lower=0> N6; // number of observations for model 6
  //
  int<lower=0> K1; // number of predictors for model 1
  int<lower=0> K2; // number of predictors for model 2
  int<lower=0> K3; // number of predictors for model 3
  int<lower=0> K4; // number of predictors for model 4
  int<lower=0> K5; // number of predictors for model 5
  int<lower=0> K6; // number of predictors for model 6
  //
  matrix[N1, K1] X1; // predictor matrix for model 1
  matrix[N2, K2] X2; // predictor matrix for model 2
  matrix[N3, K3] X3; // predictor matrix for model 3
  matrix[N4, K4] X4; // predictor matrix for model 4
  matrix[N5, K5] X5; // predictor matrix for model 5
  matrix[N6, K6] X6; // predictor matrix for model 6
  //
  array[N1] int<lower=0, upper=1> y1; // binary outcome for model 1
  array[N2] int<lower=0, upper=1> y2; // binary outcome for model 2
  array[N3] int<lower=0, upper=1> y3; // binary outcome for model 3
  array[N4] int<lower=0, upper=1> y4; // binary outcome for model 4
  array[N5] int<lower=0, upper=1> y5; // binary outcome for model 5
  array[N6] int<lower=0, upper=1> y6; // binary outcome for model 6
}
parameters {
  vector[K1] beta1; // coefficients for model 1
  vector[K2] beta2; // coefficients for model 2
  vector[K3] beta3; // coefficients for model 3
  vector[K4] beta4; // coefficients for model 4
  vector[K5] beta5; // coefficients for model 5
  vector[K6] beta6; // coefficients for model 6
}
model {
  // Priors
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  beta3 ~ normal(0, 5);
  beta4 ~ normal(0, 5);
  beta5 ~ normal(0, 5);
  beta6 ~ normal(0, 5);
  
  // Likelihood
  y1 ~ bernoulli_logit(X1 * beta1);
  y2 ~ bernoulli_logit(X2 * beta2);
  y3 ~ bernoulli_logit(X3 * beta3);
  y4 ~ bernoulli_logit(X4 * beta4);
  y5 ~ bernoulli_logit(X5 * beta5);
  y6 ~ bernoulli_logit(X6 * beta6);
}
generated quantities {
  real auc1;
  real auc2;
  real auc3;
  real auc4;
  real auc5;
  real auc6;
  real auc_diff_params_prl_wcst;
  real auc_diff_params_prl_ts;
  real auc_diff_params_ts_wcst;
  real auc_diff_behav_indices_prl_wcst;
  real auc_diff_behav_indices_prl_ts;
  real auc_diff_behav_indices_ts_wcst;
  real auc_diff_params_indices_prl;
  real auc_diff_params_indices_wcst;
  real auc_diff_params_indices_ts;
  
  vector[N1] y_pred1 = inv_logit(X1 * beta1);
  vector[N2] y_pred2 = inv_logit(X2 * beta2);
  vector[N3] y_pred3 = inv_logit(X3 * beta3);
  vector[N4] y_pred4 = inv_logit(X4 * beta4);
  vector[N5] y_pred5 = inv_logit(X5 * beta5);
  vector[N6] y_pred6 = inv_logit(X6 * beta6);
  
  auc1 = calculate_auc(y_pred1, y1);
  auc2 = calculate_auc(y_pred2, y2);
  auc3 = calculate_auc(y_pred3, y3);
  auc4 = calculate_auc(y_pred4, y4);
  auc5 = calculate_auc(y_pred5, y5);
  auc6 = calculate_auc(y_pred6, y6);
  
  // M1: WCST Steinke params
  // M2: PRL params
  // M3: WCST behav indices
  // M4: Task Switching DDM params
  // M5: Task Switching behav indices
  // M6: PRL behav indices
  //
  auc_diff_params_prl_wcst = auc2 - auc1; // AUC Difference params PRL - WCST
  auc_diff_params_prl_ts = auc2 - auc4; // AUC Difference params PRL - TS
  auc_diff_params_ts_wcst = auc4 - auc1; // AUC Difference params TS - WCST
  //
  auc_diff_behav_indices_prl_wcst = auc6 - auc3; // AUC Difference behav indices PRL - WCST
  auc_diff_behav_indices_prl_ts = auc6 - auc5; // AUC Difference behav indices PRL - TS
  auc_diff_behav_indices_ts_wcst = auc5 - auc3; // AUC Difference behav indices TS - WCST
  //
  auc_diff_params_indices_prl = auc2 - auc6; // AUC Difference params and behav indices PRL
  auc_diff_params_indices_wcst = auc1 - auc3; // AUC Difference params and behav indices WCST
  auc_diff_params_indices_ts = auc4 - auc5; // AUC Difference params and behav indices TS
}
