// Code is based on R package hBayesDM
// Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing neuro-computational mechanisms of reinforcement learning and decision-making with the hBayesDM package. 
// Computational Psychiatry, 1, 24-57. https://doi.org/10.1162/CPSY_a_00002
// Edited by Alexander Steinke
// 11th July, 2019 
// steinke.alexander@mh-hannover.de
data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Number of trials
  array[N] int<lower=1, upper=T> Tsubj; // Trials per subject
  array[N] int<lower=0, upper=1> holdout; // Holdout subjects
  array[N] int<lower=1, upper=3> group; // Group identifier
  array[N, T] int rew; // real
  array[N, T] int los; // real
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}
transformed data {
  
}
parameters {
  // Group-specific hyperparameters
  array[3] vector[8] mu_p; // One for each group
  array[3] vector<lower=0>[8] sigma; // One for each group
  
  // Subject-level raw parameters
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
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N] MB_Arew;
  vector<lower=0, upper=1>[N] MB_Apun;
  vector<lower=0, upper=1>[N] MB_gamma;
  vector<lower=0, upper=1>[N] MF_Arew;
  vector<lower=0, upper=1>[N] MF_Apun;
  vector<lower=0, upper=1>[N] MF_gamma;
  vector<lower=0>[N] temp;
  vector<lower=0, upper=1>[N] w;
  
  for (i in 1 : N) {
    int grp = group[i]; // Identify the group of the subject
    // Apply group-specific transformations
    MB_Arew[i] = Phi_approx(mu_p[grp][1] + sigma[grp][1] * MB_Arew_pr[i]);
    MB_Apun[i] = Phi_approx(mu_p[grp][2] + sigma[grp][2] * MB_Apun_pr[i]);
    MB_gamma[i] = Phi_approx(mu_p[grp][3] + sigma[grp][3] * MB_gamma_pr[i]);
    MF_Arew[i] = Phi_approx(mu_p[grp][4] + sigma[grp][4] * MF_Arew_pr[i]);
    MF_Apun[i] = Phi_approx(mu_p[grp][5] + sigma[grp][5] * MF_Apun_pr[i]);
    MF_gamma[i] = Phi_approx(mu_p[grp][6] + sigma[grp][6] * MF_gamma_pr[i]);
    temp[i] = Phi_approx(mu_p[grp][7] + sigma[grp][7] * temp_pr[i]) * 5;
    w[i] = Phi_approx(mu_p[grp][8] + sigma[grp][8] * w_pr[i]);
  }
}
model {
  // Priors for group-specific hyperparameters
  for (g in 1 : 3) {
    mu_p[g] ~ normal(0, 1.0);
    sigma[g] ~ cauchy(0, 5.0);
  }
  
  // Priors for subject-level parameters
  MB_Arew_pr ~ normal(0, 1.0);
  MB_Apun_pr ~ normal(0, 1.0);
  MB_gamma_pr ~ normal(0, 1.0);
  MF_Arew_pr ~ normal(0, 1.0);
  MF_Apun_pr ~ normal(0, 1.0);
  MF_gamma_pr ~ normal(0, 1.0);
  temp_pr ~ normal(0, 1.0);
  w_pr ~ normal(0, 1.0);
  
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      // Consider subject only if it is not in the holdout data set
      
      // Define values
      vector[4] Qeff; // Resulting Q values
      vector[3] QMB; // Model-based Q values
      vector[4] QMF; // Model-free Q values
      vector[3] PEMB; // Model-based prediction error 
      vector[4] PEMF; // Model-free prediction error 
      
      // Initialize values
      Qeff = rep_vector(0.0, 4);
      QMB = rep_vector(0.0, 3);
      QMF = rep_vector(0.0, 4);
      
      for (t in 1 : Tsubj[i]) {
        // Resulting Q values
        
        Qeff = rep_vector(-1.0, 4); // Set odd key cards to -1, all other key cards to the sum of matching Q-values
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        Qeff[resp_color[i, t]] = Qeff[resp_color[i, t]] + w[i] * QMB[1];
        Qeff[resp_shape[i, t]] = Qeff[resp_shape[i, t]] + w[i] * QMB[2];
        Qeff[resp_number[i, t]] = Qeff[resp_number[i, t]] + w[i] * QMB[3];
        Qeff = Qeff + (1 - w[i]) * QMF;
        
        // Temperature 
        Qeff = softmax(Qeff / temp[i]);
        
        // Softmax choice
        resp_choice[i, t] ~ categorical(Qeff);
        
        // Decay 
        QMB = MB_gamma[i] * QMB;
        QMF = MF_gamma[i] * QMF;
        
        // Prediction error signals
        PEMB = rep_vector(0.0, 3); // init PE vectors
        PEMF = rep_vector(0.0, 4);
        PEMB[rule_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMB[rule_choice[i, t]];
        PEMF[resp_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMF[resp_choice[i, t]];
        
        // Update Q values
        if (rew[i, t] == 1) {
          QMB = QMB + MB_Arew[i] * PEMB;
          QMF = QMF + MF_Arew[i] * PEMF;
        } else {
          QMB = QMB + MB_Apun[i] * PEMB;
          QMF = QMF + MF_Apun[i] * PEMF;
        }
      }
    }
  }
}
generated quantities {
  // Group-specific mean parameters
  array[3] real<lower=0, upper=1> mu_MB_Arew;
  array[3] real<lower=0, upper=1> mu_MB_Apun;
  array[3] real<lower=0, upper=1> mu_MB_gamma;
  
  array[3] real<lower=0, upper=1> mu_MF_Arew;
  array[3] real<lower=0, upper=1> mu_MF_Apun;
  array[3] real<lower=0, upper=1> mu_MF_gamma;
  
  array[3] real<lower=0> mu_temp;
  array[3] real<lower=0, upper=1> mu_w;
  
  // For log likelihood calculation
  array[N] real log_lik;
  
  // For posterior predictive check
  array[N, T] real y_pred;
  
  // Initialize posterior predictions
  for (i in 1 : N) {
    for (t in 1 : T) {
      y_pred[i, t] = -1;
    }
  }
  
  // Calculate group-specific means
  for (g in 1 : 3) {
    mu_MB_Arew[g] = Phi_approx(mu_p[g][1]);
    mu_MB_Apun[g] = Phi_approx(mu_p[g][2]);
    mu_MB_gamma[g] = Phi_approx(mu_p[g][3]);
    mu_MF_Arew[g] = Phi_approx(mu_p[g][4]);
    mu_MF_Apun[g] = Phi_approx(mu_p[g][5]);
    mu_MF_gamma[g] = Phi_approx(mu_p[g][6]);
    mu_temp[g] = Phi_approx(mu_p[g][7]) * 5;
    mu_w[g] = Phi_approx(mu_p[g][8]);
  }
  
  {
    // Local section for efficiency
    for (i in 1 : N) {
      // Define values
      vector[4] Qeff; // Resulting Q values
      vector[3] QMB; // Model-based Q values
      vector[4] QMF; // Model-free Q values
      vector[3] PEMB; // Model-based prediction error
      vector[4] PEMF; // Model-free prediction error
      
      // Initialize values
      Qeff = rep_vector(0.0, 4);
      QMB = rep_vector(0.0, 3);
      QMF = rep_vector(0.0, 4);
      log_lik[i] = 0.0;
      
      for (t in 1 : Tsubj[i]) {
        // Computation of Q values
        Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        Qeff[resp_color[i, t]] = Qeff[resp_color[i, t]] + w[i] * QMB[1];
        Qeff[resp_shape[i, t]] = Qeff[resp_shape[i, t]] + w[i] * QMB[2];
        Qeff[resp_number[i, t]] = Qeff[resp_number[i, t]] + w[i] * QMB[3];
        Qeff = Qeff + (1 - w[i]) * QMF;
        
        // Temperature adjustment
        Qeff = softmax(Qeff / temp[i]);
        
        // Compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_lpmf(resp_choice[i, t] | Qeff);
        
        // Draw posterior prediction for current trial
        y_pred[i, t] = categorical_rng(Qeff);
        
        // Decay and update Q values
        QMB = MB_gamma[i] * QMB;
        QMF = MF_gamma[i] * QMF;
        PEMB = rep_vector(0.0, 3);
        PEMF = rep_vector(0.0, 4);
        PEMB[rule_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMB[rule_choice[i, t]];
        PEMF[resp_choice[i, t]] = (rew[i, t] + los[i, t])
                                  - QMF[resp_choice[i, t]];
        if (rew[i, t] == 1) {
          QMB = QMB + MB_Arew[i] * PEMB;
          QMF = QMF + MF_Arew[i] * PEMF;
        } else {
          QMB = QMB + MB_Apun[i] * PEMB;
          QMF = QMF + MF_Apun[i] * PEMF;
        }
      }
    }
  }
}


