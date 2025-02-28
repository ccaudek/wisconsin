data {
  int<lower=1> N; // #subjects
  int<lower=1> T; // #trials
  array[N] int<lower=1, upper=T> Tsubj;
  array[N] int<lower=0, upper=1> holdout;
  
  array[N, T] real rew;
  array[N, T] real los;
  
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}
parameters {
  // Group-level means
  vector[8] mu_p;
  
  // Group-level standard deviations
  // Use a half-normal (or normal(0,0.5) truncated to positive) for more constrained priors
  vector<lower=0>[8] sigma;
  
  // Raw subject-level parameters
  vector[N] MB_Arew_pr;
  vector[N] MB_Apun_pr;
  vector[N] MB_gamma_pr;
  
  vector[N] MF_Arew_pr;
  vector[N] MF_Apun_pr;
  vector[N] MF_gamma_pr;
  
  vector[N] temp_pr;
  vector[N] w_pr;
}
transformed parameters {
  // Subject-level parameters
  // Use inv_logit for stable transformations from R -> (0,1)
  vector<lower=0, upper=1>[N] MB_Arew = inv_logit(mu_p[1]
                                                  + sigma[1] * MB_Arew_pr);
  vector<lower=0, upper=1>[N] MB_Apun = inv_logit(mu_p[2]
                                                  + sigma[2] * MB_Apun_pr);
  vector<lower=0, upper=1>[N] MB_gamma = inv_logit(mu_p[3]
                                                   + sigma[3] * MB_gamma_pr);
  
  vector<lower=0, upper=1>[N] MF_Arew = inv_logit(mu_p[4]
                                                  + sigma[4] * MF_Arew_pr);
  vector<lower=0, upper=1>[N] MF_Apun = inv_logit(mu_p[5]
                                                  + sigma[5] * MF_Apun_pr);
  vector<lower=0, upper=1>[N] MF_gamma = inv_logit(mu_p[6]
                                                   + sigma[6] * MF_gamma_pr);
  
  // Temperature mapped to something > 0
  // If you want more moderate temperatures, use a smaller multiplier than 5.
  vector<lower=0>[N] temp = exp(mu_p[7] + sigma[7] * temp_pr);
  // Using exp ensures strictly positive. If you want it strictly between 0 and 5:
  // temp = 5 * inv_logit(mu_p[7] + sigma[7] * temp_pr);
  
  vector<lower=0, upper=1>[N] w = inv_logit(mu_p[8] + sigma[8] * w_pr);
}
model {
  // Priors on group-level means
  mu_p ~ normal(0, 0.5); // more concentrated prior to keep parameters away from extremes
  
  // Priors on group-level sds 
  // Encourage small variability (reducing individual differences)
  sigma ~ normal(0, 0.5); // If needed, use sigma ~ normal(0,0.5)T[0,] for truncation
  // This keeps sigma small, reducing the scatter and hence reducing extreme values.
  
  // Subject-level effects
  MB_Arew_pr ~ normal(0, 1);
  MB_Apun_pr ~ normal(0, 1);
  MB_gamma_pr ~ normal(0, 1);
  
  MF_Arew_pr ~ normal(0, 1);
  MF_Apun_pr ~ normal(0, 1);
  MF_gamma_pr ~ normal(0, 1);
  
  temp_pr ~ normal(0, 1);
  w_pr ~ normal(0, 1);
  
  // Likelihood
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      vector[4] Qeff;
      vector[3] QMB = rep_vector(0.0, 3);
      vector[4] QMF = rep_vector(0.0, 4);
      
      for (t in 1 : Tsubj[i]) {
        // Compute Qeff
        Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        
        Qeff[resp_color[i, t]] += w[i] * QMB[1];
        Qeff[resp_shape[i, t]] += w[i] * QMB[2];
        Qeff[resp_number[i, t]] += w[i] * QMB[3];
        Qeff += (1 - w[i]) * QMF;
        
        // Softmax choice
        vector[4] p = softmax(Qeff / temp[i]);
        resp_choice[i, t] ~ categorical(p);
        
        // Decay
        QMB *= MB_gamma[i];
        QMF *= MF_gamma[i];
        
        // Compute PEs
        vector[3] PEMB = rep_vector(0.0, 3);
        vector[4] PEMF = rep_vector(0.0, 4);
        
        PEMB[rule_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMB[rule_choice[i, t]];
        PEMF[resp_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMF[resp_choice[i, t]];
        
        // Update Q-values
        if (rew[i, t] == 1) {
          QMB += MB_Arew[i] * PEMB;
          QMF += MF_Arew[i] * PEMF;
        } else {
          QMB += MB_Apun[i] * PEMB;
          QMF += MF_Apun[i] * PEMF;
        }
      }
    }
  }
}
generated quantities {
  real<lower=0, upper=1> mu_MB_Arew = inv_logit(mu_p[1]);
  real<lower=0, upper=1> mu_MB_Apun = inv_logit(mu_p[2]);
  real<lower=0, upper=1> mu_MB_gamma = inv_logit(mu_p[3]);
  
  real<lower=0, upper=1> mu_MF_Arew = inv_logit(mu_p[4]);
  real<lower=0, upper=1> mu_MF_Apun = inv_logit(mu_p[5]);
  real<lower=0, upper=1> mu_MF_gamma = inv_logit(mu_p[6]);
  
  // Using exp transform for temp:
  // If you changed temp parameterization above, reflect that here:
  real<lower=0> mu_temp = exp(mu_p[7]);
  // Or if using the logistic transform for temp:
  // real<lower=0> mu_temp = 5 * inv_logit(mu_p[7]);
  
  real<lower=0, upper=1> mu_w = inv_logit(mu_p[8]);
  
  // Log likelihood and posterior predictive
  array[N] real log_lik;
  array[N, T] int y_pred;
  
  for (i in 1 : N) {
    log_lik[i] = 0;
    for (t in 1 : T) {
      y_pred[i, t] = -1;
    }
  }
  
  for (i in 1 : N) {
    vector[4] QMF = rep_vector(0.0, 4);
    vector[3] QMB = rep_vector(0.0, 3);
    
    if (holdout[i] == 0) {
      for (t in 1 : Tsubj[i]) {
        vector[4] Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        
        Qeff[resp_color[i, t]] += w[i] * QMB[1];
        Qeff[resp_shape[i, t]] += w[i] * QMB[2];
        Qeff[resp_number[i, t]] += w[i] * QMB[3];
        Qeff += (1 - w[i]) * QMF;
        
        vector[4] p = softmax(Qeff / temp[i]);
        log_lik[i] += categorical_lpmf(resp_choice[i, t] | p);
        
        y_pred[i, t] = categorical_rng(p);
        
        QMB *= MB_gamma[i];
        QMF *= MF_gamma[i];
        
        vector[3] PEMB = rep_vector(0.0, 3);
        vector[4] PEMF = rep_vector(0.0, 4);
        
        PEMB[rule_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMB[rule_choice[i, t]];
        PEMF[resp_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMF[resp_choice[i, t]];
        
        if (rew[i, t] == 1) {
          QMB += MB_Arew[i] * PEMB;
          QMF += MF_Arew[i] * PEMF;
        } else {
          QMB += MB_Apun[i] * PEMB;
          QMF += MF_Apun[i] * PEMF;
        }
      }
    }
  }
}
