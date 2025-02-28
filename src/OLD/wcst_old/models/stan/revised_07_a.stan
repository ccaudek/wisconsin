// Enhanced Stan model for RL fitting with improved hierarchical variability and noise handling
data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Number of trials
  array[N] int<lower=1, upper=T> Tsubj; // Trials per subject
  array[N] int<lower=0, upper=1> holdout; // Holdout indicator
  
  array[N, T] real rew;
  array[N, T] real los;
  
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}
parameters {
  // Group-level parameters
  vector[8] mu_p; // Group-level means
  vector<lower=0>[8] sigma; // Group-level standard deviations
  
  // Raw subject-level deviations (non-centered parameterization)
  array[8] vector[N] raw_subject_offsets; // Subject-level deviations
  
  // Noise scaling for likelihood
  real<lower=0> epsilon; // Additional noise parameter
  
  // Hyperparameters for mixture modeling
  simplex[2] group_weights; // Mixture weights for latent subgroups
  array[8] vector[2] group_mu_p; // Means for latent subgroups
  array[8] vector<lower=0>[2] group_sigma; // Standard deviations for latent subgroups
}
transformed parameters {
  // Transformed subject-level parameters
  vector[N] MB_Arew;
  vector[N] MB_Apun;
  vector[N] MB_gamma;
  vector[N] MF_Arew;
  vector[N] MF_Apun;
  vector[N] MF_gamma;
  vector[N] temp;
  vector[N] w;
  
  for (i in 1 : N) {
    MB_Arew[i] = inv_logit(mu_p[1] + sigma[1] * raw_subject_offsets[1][i]);
    MB_Apun[i] = inv_logit(mu_p[2] + sigma[2] * raw_subject_offsets[2][i]);
    MB_gamma[i] = inv_logit(mu_p[3] + sigma[3] * raw_subject_offsets[3][i]);
    
    MF_Arew[i] = inv_logit(mu_p[4] + sigma[4] * raw_subject_offsets[4][i]);
    MF_Apun[i] = inv_logit(mu_p[5] + sigma[5] * raw_subject_offsets[5][i]);
    MF_gamma[i] = inv_logit(mu_p[6] + sigma[6] * raw_subject_offsets[6][i]);
    
    temp[i] = exp(mu_p[7] + sigma[7] * raw_subject_offsets[7][i]);
    w[i] = inv_logit(mu_p[8] + sigma[8] * raw_subject_offsets[8][i]);
  }
}
model {
  // Priors on group-level parameters
  for (g in 1 : 2) {
    group_mu_p[ : , g] ~ normal(0, 1); // Priors for group means
    group_sigma[ : , g] ~ cauchy(0, 2); // Priors for group-level variability
  }
  group_weights ~ dirichlet(rep_vector(1.0, 2)); // Dirichlet prior for mixture weights
  
  mu_p ~ normal(0, 1); // Priors for overall means
  sigma ~ cauchy(0, 2); // Priors for overall variability
  
  // Priors on subject-level deviations
  for (k in 1 : 8) {
    raw_subject_offsets[k] ~ normal(0, 1);
  }
  
  // Prior for additional noise
  epsilon ~ normal(0, 1); // Regularizes additional noise in likelihood
  
  // Likelihood
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      vector[4] Qeff;
      vector[3] QMB;
      vector[4] QMF;
      vector[3] PEMB;
      vector[4] PEMF;
      
      Qeff = rep_vector(0.0, 4);
      QMB = rep_vector(0.0, 3);
      QMF = rep_vector(0.0, 4);
      
      for (t in 1 : Tsubj[i]) {
        Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = w[i] * QMB[1];
        Qeff[resp_shape[i, t]] += w[i] * QMB[2];
        Qeff[resp_number[i, t]] += w[i] * QMB[3];
        Qeff += (1 - w[i]) * QMF;
        
        Qeff = softmax(Qeff / (temp[i] + epsilon)); // Apply softmax with additional noise
        resp_choice[i, t] ~ categorical(Qeff);
        
        QMB *= MB_gamma[i];
        QMF *= MF_gamma[i];
        
        PEMB = rep_vector(0.0, 3);
        PEMF = rep_vector(0.0, 4);
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
generated quantities {
  // Group-level means for reporting
  real<lower=0, upper=1> mu_MB_Arew = inv_logit(mu_p[1]);
  real<lower=0, upper=1> mu_MB_Apun = inv_logit(mu_p[2]);
  real<lower=0, upper=1> mu_MB_gamma = inv_logit(mu_p[3]);
  
  real<lower=0, upper=1> mu_MF_Arew = inv_logit(mu_p[4]);
  real<lower=0, upper=1> mu_MF_Apun = inv_logit(mu_p[5]);
  real<lower=0, upper=1> mu_MF_gamma = inv_logit(mu_p[6]);
  
  real<lower=0> mu_temp = exp(mu_p[7]);
  real<lower=0, upper=1> mu_w = inv_logit(mu_p[8]);
  
  // Log-likelihood for each subject
  array[N] real log_lik;
  array[N, T] real y_pred;
  
  for (i in 1 : N) {
    log_lik[i] = 0;
    for (t in 1 : Tsubj[i]) {
      vector[4] Qeff;
      vector[3] QMB;
      vector[4] QMF;
      
      Qeff = rep_vector(-1.0, 4);
      QMB = rep_vector(0.0, 3);
      QMF = rep_vector(0.0, 4);
      
      Qeff[resp_color[i, t]] = w[i] * QMB[1];
      Qeff[resp_shape[i, t]] += w[i] * QMB[2];
      Qeff[resp_number[i, t]] += w[i] * QMB[3];
      Qeff += (1 - w[i]) * QMF;
      
      Qeff = softmax(Qeff / (temp[i] + epsilon));
      log_lik[i] += categorical_lpmf(resp_choice[i, t] | Qeff);
      
      y_pred[i, t] = categorical_rng(Qeff);
    }
  }
}
