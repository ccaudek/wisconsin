// Code is based on R package hBayesDM
// Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing neuro-computational mechanisms of reinforcement learning and decision-making with the hBayesDM package. 
// Computational Psychiatry, 1, 24-57. https://doi.org/10.1162/CPSY_a_00002
// Edited by Alexander Steinke
// 11th July, 2019 
// steinke.alexander@mh-hannover.de
data {
  int<lower=1> N; // How many subjects? 
  int<lower=1> T; // How many trials?
  int<lower=1, upper=T>  Tsubj[N];   // How many trials per data set?
  int<lower=0, upper=1>  holdout[N]; // Should the participant be holdout for cross validation?

  real rew[N, T];
  real los[N, T];
  
  int rule_choice[N, T];
  int resp_choice[N, T];
  
  int resp_color[N, T];
  int resp_shape[N, T];
  int resp_number[N, T];
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters  
  vector[4] mu_p; 
  vector<lower=0>[4] sigma;

  // Individual-level raw parameters (for Matt trick)
  vector[N] MB_Arew_pr;
  vector[N] MB_Apun_pr;
  vector[N] MB_gamma_pr;
  vector[N] temp_pr;
}
transformed parameters {
  // Transform individual-level raw parameters
  vector<lower=0, upper=1>[N] MB_Arew;
  vector<lower=0, upper=1>[N] MB_Apun;
  vector<lower=0, upper=1>[N] MB_gamma;
  vector<lower=0> [N] temp;

  for (i in 1:N) {
    MB_Arew[i]  = Phi_approx( mu_p[1] + sigma[1] * MB_Arew_pr[i] );
    MB_Apun[i]  = Phi_approx( mu_p[2] + sigma[2] * MB_Apun_pr[i] );
    MB_gamma[i] = Phi_approx( mu_p[3] + sigma[3] * MB_gamma_pr[i] );
    temp[i]     = Phi_approx( mu_p[4] + sigma[4] * temp_pr[i] )*5;
  }
}
model {
  // Hyperparameters
  mu_p     ~ normal(0, 1.0);
  sigma    ~ cauchy(0, 5.0);
  
  // Individual-level parameters
  MB_Arew_pr  ~ normal(0, 1.0);
  MB_Apun_pr  ~ normal(0, 1.0);
  MB_gamma_pr ~ normal(0, 1.0);
  temp_pr     ~ normal(0, 1.0);
  
  for (i in 1:N) {
    if( holdout[i]==0 ) { // Consider subject only if it is not in the holdout data set
      // Define values
      vector[4] Qeff; // Resulting Q values
      vector[3] QMB;  // Model-based Q values
      vector[3] PEMB; // Model-based predicition error
      real sumQ; 
      
      // Initialize values
      Qeff  = rep_vector(0.0, 4); 
      QMB = rep_vector(0.0, 3); 
  
      for (t in 1:Tsubj[i]) {
      
        // Resulting Q values
        Qeff = rep_vector(-1.0, 4); // Set odd key cards to -1, all other key cards to the sum of matching Q-values
        Qeff[ resp_color[i,t] ]   = 0.0;
        Qeff[ resp_shape[i,t] ]   = 0.0;
        Qeff[ resp_number[i,t] ]  = 0.0;
        Qeff[ resp_color[i,t] ]   = Qeff[ resp_color[i,t] ]  + QMB[1];
        Qeff[ resp_shape[i,t] ]   = Qeff[ resp_shape[i,t] ]  + QMB[2];
        Qeff[ resp_number[i,t] ]  = Qeff[ resp_number[i,t] ] + QMB[3];
        
        // Temperature 
        Qeff = softmax(Qeff/temp[i]);
        
        // Softmax choice
        resp_choice[i, t] ~ categorical( Qeff );
        
        // Decay parameter
        QMB = MB_gamma[i] * QMB;
          
        // Prediction error signal
        PEMB = rep_vector(0.0, 3); // init PE vectors
        PEMB[rule_choice[i,t]] = (rew[i,t]+los[i,t]) - QMB[rule_choice[i,t]];
  
        // Update Q values
        if (rew[i,t]==1) {
          QMB = QMB + MB_Arew[i] * PEMB;
        } else {
          QMB = QMB + MB_Apun[i] * PEMB;
        }
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_MB_Arew;
  real<lower=0, upper=1> mu_MB_Apun;
  real<lower=0, upper=1> mu_MB_gamma;
  real<lower=0> mu_temp;

  // For log likelihood calculation
  real log_lik[N];
  
  // For posterior predictive check
  // real y_pred[N,T]; 
  
  // Set all posterior predictions to 0 (avoids NULL values)
  //for (i in 1:N) {
  //  for (t in 1:T) {
  //    y_pred[i,t] = -1;
  //  }
  //}

  mu_MB_Arew  = Phi_approx(mu_p[1]);
  mu_MB_Apun  = Phi_approx(mu_p[2]);
  mu_MB_gamma = Phi_approx(mu_p[3]);
  mu_temp     = Phi_approx(mu_p[4])*5;

  { // Local section, this saves time and space
    for (i in 1:N) {
      // Define values
      vector[4] Qeff; // Resulting Q values
      vector[3] QMB;  // Model-baesd Q values; 1=color; 2=form; 3=number, 4=odd
      vector[3] PEMB; // Model-based prediction error
      real sumQ; 

      // Initialize values
      Qeff = rep_vector(0.0, 4); 
      QMB  = rep_vector(0.0, 3); 
      log_lik[i] = 0.0;
  
      for (t in 1:Tsubj[i]) {
      
        // Resulting Q-values
        Qeff = rep_vector(-1.0, 4); // Set odd key cards to -1, all other key cards to the sum of matching Q-values
        Qeff[ resp_color[i,t] ]   = 0.0;
        Qeff[ resp_shape[i,t] ]   = 0.0;
        Qeff[ resp_number[i,t] ]  = 0.0;
        Qeff[ resp_color[i,t] ]   = Qeff[ resp_color[i,t] ]  + QMB[1];
        Qeff[ resp_shape[i,t] ]   = Qeff[ resp_shape[i,t] ]  + QMB[2];
        Qeff[ resp_number[i,t] ]  = Qeff[ resp_number[i,t] ] + QMB[3];
        
        // Temperature 
        Qeff = softmax(Qeff/temp[i]);
        
        // Compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_lpmf( resp_choice[i, t] | Qeff );
        
        // Draw posterior prediction for current trial
        // y_pred[i,t] = categorical_rng(Qeff);
      
        // Decay parameter
        QMB = MB_gamma[i] * QMB;
          
        // Prediction error signals
        PEMB = rep_vector(0.0, 3); // init PE vectors
        PEMB[rule_choice[i,t]] = (rew[i,t]+los[i,t]) - QMB[rule_choice[i,t]];
  
        // Update Q values
        if (rew[i,t]==1) {
            QMB = QMB + MB_Arew[i] * PEMB;
        } else {
            QMB = QMB + MB_Apun[i] * PEMB;
        }
      }
    }
  }
}
