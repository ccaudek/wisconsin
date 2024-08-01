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
}
model {
  // Priors (less informative)
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  
  // Likelihood
  y1 ~ bernoulli_logit(X1 * beta1);
  y2 ~ bernoulli_logit(X2 * beta2);
}
generated quantities {
  real auc1;
  real auc2;
  real auc_diff;
  
  // AUC computation using ROC curve method
  {
    vector[N1] p1 = inv_logit(X1 * beta1);
    vector[N2] p2 = inv_logit(X2 * beta2);
    real tp1 = 0;
    real fp1 = 0;
    real tn1 = 0;
    real fn1 = 0;
    real tp2 = 0;
    real fp2 = 0;
    real tn2 = 0;
    real fn2 = 0;
    
    for (i in 1 : N1) {
      if (y1[i] == 1 && p1[i] >= 0.5) 
        tp1 += 1;
      if (y1[i] == 0 && p1[i] >= 0.5) 
        fp1 += 1;
      if (y1[i] == 0 && p1[i] < 0.5) 
        tn1 += 1;
      if (y1[i] == 1 && p1[i] < 0.5) 
        fn1 += 1;
    }
    
    for (i in 1 : N2) {
      if (y2[i] == 1 && p2[i] >= 0.5) 
        tp2 += 1;
      if (y2[i] == 0 && p2[i] >= 0.5) 
        fp2 += 1;
      if (y2[i] == 0 && p2[i] < 0.5) 
        tn2 += 1;
      if (y2[i] == 1 && p2[i] < 0.5) 
        fn2 += 1;
    }
    
    real tpr1 = tp1 / (tp1 + fn1);
    real fpr1 = fp1 / (fp1 + tn1);
    real tpr2 = tp2 / (tp2 + fn2);
    real fpr2 = fp2 / (fp2 + tn2);
    
    auc1 = (tpr1 + (1 - fpr1)) / 2;
    auc2 = (tpr2 + (1 - fpr2)) / 2;
  }
  
  auc_diff = auc2 - auc1; // Difference in AUC (task2 - task1)
}
