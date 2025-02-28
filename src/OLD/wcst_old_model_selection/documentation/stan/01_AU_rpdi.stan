// Code is based on R package hBayesDM
// Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing neuro-computational mechanisms of reinforcement learning and decision-making with the hBayesDM package. 
// Computational Psychiatry, 1, 24-57. https://doi.org/10.1162/CPSY_a_00002
// Edited by Alexander Steinke
// 11th July, 2019 
// steinke.alexander@mh-hannover.de
data {
  int<lower=1> N; // how many subjects? 
  int<lower=1> T; // how many trials?
  int<lower=1, upper=T>  Tsubj[N];   //how many trials per data set?
  int<lower=0, upper=1>  holdout[N]; //should the subject be holdout?
  
  real rew[N, T];
  real los[N, T];
  
  int rule_choice[N, T];
  int resp_choice[N, T];
  
  int resp_color[N, T];
  int resp_shape[N, T];
  int resp_number[N, T];
  
  int odd_response[N, T];
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters  
  vector[4] mu_par; 
  vector<lower=0>[4] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[N] r_pr;
  vector[N] p_pr;
  vector[N] d_pr;
  vector[N] i_pr;

}
transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N] r;
  vector<lower=0, upper=1>[N] p;
  vector<lower=0, upper=5>[N] d;
  vector<lower=0, upper=5>[N] i;

  for (j in 1:N) {
    r[j] = Phi_approx( mu_par[1] + sigma[1] * r_pr[j] );
    p[j] = Phi_approx( mu_par[2] + sigma[2] * p_pr[j] );
    d[j] = Phi_approx( mu_par[3] + sigma[3] * d_pr[j] ) * 5;
    i[j] = Phi_approx( mu_par[4] + sigma[4] * i_pr[j] ) * 5;
  }
}
model {
  // Hyperparameters
  mu_par   ~ normal(0, 2.0);
  sigma    ~ cauchy(0, 5.0);
  
  // individual parameters
  r_pr  ~ normal(0, 2.0);
  p_pr  ~ normal(0, 2.0);
  d_pr  ~ normal(0, 2.0);
  i_pr  ~ normal(0, 2.0);

  for (j in 1:N) {
    if( holdout[j]==0 ) { // consider subject only if it is not in the holdout data set
      // Define values
      vector[4] Pr;  // Resulting response probailities
      vector[3] A; //rule level A values; 1=color; 2=form; 3=number
      vector[3] S; 
      
      real sumA;
      real sumS;
  
      // Initialize values
      Pr = rep_vector(0.0, 4); 
      A  = rep_vector(0.33, 3); 
      S  = rep_vector(0.0, 3); 
      
      // trial 1
      // get response probabilities
        sumA = pow(A[1],d[j]) + pow(A[2],d[j]) + pow(A[3],d[j]);
        Pr[ resp_color[j,1] ]   = Pr[ resp_color[j,1] ] + pow(A[1],d[j]);
        Pr[ resp_shape[j,1] ]   = Pr[ resp_shape[j,1] ] + pow(A[2],d[j]);
        Pr[ resp_number[j,1] ]  = Pr[ resp_number[j,1] ] + pow(A[3],d[j]);
        for (k in 1:4) Pr[k] = Pr[k]/sumA;
      
      resp_choice[j, 1] ~ categorical( Pr );
      
      for (t in 2:Tsubj[j]) {
          
        if (rew[j,t-1]==1) {
          S = rep_vector(0.0, 3);
          S[ rule_choice[j,t-1] ] = 1;
          
          sumS = S[1]*pow(A[1],i[j]) + S[2]*pow(A[2],i[j]) + S[3]*pow(A[3],i[j]);
          S[1] = S[1]*pow(A[1],i[j]) / sumS;
          S[2] = S[2]*pow(A[2],i[j]) / sumS;
          S[3] = S[3]*pow(A[3],i[j]) / sumS;
          
          A = (1-r[j]) * A + r[j] * S;
          
        } else {
          S = rep_vector(1.0, 3);
          S[ rule_choice[j,t-1] ] = 0;
          
          sumS = S[1]*pow(A[1],i[j]) + S[2]*pow(A[2],i[j]) + S[3]*pow(A[3],i[j]);
          S[1] = S[1]*pow(A[1],i[j]) / sumS;
          S[2] = S[2]*pow(A[2],i[j]) / sumS;
          S[3] = S[3]*pow(A[3],i[j]) / sumS;
            
          A = (1-p[j]) * A + p[j] * S;
        }
          
        // Response probabilities
        Pr = rep_vector(0.0, 4);
        sumA = pow(A[1],d[j]) + pow(A[2],d[j]) + pow(A[3],d[j]);
        Pr[ resp_color[j,t] ]   = Pr[ resp_color[j,t] ] + pow(A[1],d[j]);
        Pr[ resp_shape[j,t] ]   = Pr[ resp_shape[j,t] ] + pow(A[2],d[j]);
        Pr[ resp_number[j,t] ]  = Pr[ resp_number[j,t] ] + pow(A[3],d[j]);
        for (k in 1:4) Pr[k] = Pr[k]/sumA;
        
        resp_choice[j, t] ~ categorical( Pr );
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_r;
  real<lower=0, upper=1> mu_p;
  real<lower=0, upper=5> mu_d;
  real<lower=0, upper=5> mu_i;

  // For log likelihood calculation
  real log_lik[N];
  
  // For posterior predictive check
  real y_pred[N,T]; 

  // Set all posterior predictions to 0 (avoids NULL values)
  for (j in 1:N) {
    for (t in 1:T) {
      y_pred[j,t] = -1;
    }
  }

  mu_r  = Phi_approx(mu_par[1]);
  mu_p  = Phi_approx(mu_par[2]);
  mu_d  = Phi_approx(mu_par[3]) *5;
  mu_i  = Phi_approx(mu_par[4]) *5;

  { // local section, this saves time and space
    for (j in 1:N) {
      // Define values
      vector[4] Pr;  // Resulting response probailities
      vector[3] A; //rule level A values; 1=color; 2=form; 3=number
      vector[3] S; 
      
      real sumA;
      real sumS;
  
      // Initialize values
      Pr = rep_vector(0.0, 4); 
      A  = rep_vector(0.33, 3); 
      S  = rep_vector(0.0, 3);
      log_lik[j] = 0.0;
    
      // trial 1
      // get response probabilities
        sumA = pow(A[1],d[j]) + pow(A[2],d[j]) + pow(A[3],d[j]);
        Pr[ resp_color[j,1] ]   = Pr[ resp_color[j,1] ] + pow(A[1],d[j]);
        Pr[ resp_shape[j,1] ]   = Pr[ resp_shape[j,1] ] + pow(A[2],d[j]);
        Pr[ resp_number[j,1] ]  = Pr[ resp_number[j,1] ] + pow(A[3],d[j]);
        for (k in 1:4) Pr[k] = Pr[k]/sumA;
    
      // compute log likelihood of current trial // only if not odd choice
      log_lik[j] = log_lik[j] + categorical_lpmf( resp_choice[j, 1] | Pr );
      
      // generate posterior prediction for current trial
      y_pred[j,1] = categorical_rng(Pr);

      for (t in 2:Tsubj[j]) {
        
        if (rew[j,t-1]==1) {
          S = rep_vector(0.0, 3);
          S[ rule_choice[j,t-1] ] = 1;
        
          sumS = S[1]*pow(A[1],i[j]) + S[2]*pow(A[2],i[j]) + S[3]*pow(A[3],i[j]);
          S[1] = S[1]*pow(A[1],i[j]) / sumS;
          S[2] = S[2]*pow(A[2],i[j]) / sumS;
          S[3] = S[3]*pow(A[3],i[j]) / sumS;
        
          A = (1-r[j]) * A + r[j] * S;
        
        } else {
          S = rep_vector(1.0, 3);
          S[ rule_choice[j,t-1] ] = 0;
        
          sumS = S[1]*pow(A[1],i[j]) + S[2]*pow(A[2],i[j]) + S[3]*pow(A[3],i[j]);
          S[1] = S[1]*pow(A[1],i[j]) / sumS;
          S[2] = S[2]*pow(A[2],i[j]) / sumS;
          S[3] = S[3]*pow(A[3],i[j]) / sumS;
          
          A = (1-p[j]) * A + p[j] * S;
        }
        
        // response probabilities
        Pr = rep_vector(0.0, 4);
        sumA = pow(A[1],d[j]) + pow(A[2],d[j]) + pow(A[3],d[j]);
        Pr[ resp_color[j,t] ]   = Pr[ resp_color[j,t] ] + pow(A[1],d[j]);
        Pr[ resp_shape[j,t] ]   = Pr[ resp_shape[j,t] ] + pow(A[2],d[j]);
        Pr[ resp_number[j,t] ]  = Pr[ resp_number[j,t] ] + pow(A[3],d[j]);
        for (k in 1:4) Pr[k] = Pr[k]/sumA;

        // compute log likelihood of current trial
        log_lik[j] = log_lik[j] + categorical_lpmf( resp_choice[j, t] | Pr );
      
        // generate posterior prediction for current trial
        y_pred[j,t] = categorical_rng(Pr);
      }
    }
  }
}
