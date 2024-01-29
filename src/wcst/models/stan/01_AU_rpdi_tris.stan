data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Number of trials
  array[N] int<lower=1, upper=T> Tsubj; // Number of trials per subject
  array[N] int<lower=0, upper=1> holdout; // Holdout indicator for each subject
  array[N] int<lower=1, upper=2> group; // Group identifier: 1 for control, 2 for patient
  
  array[N, T] real rew;
  array[N, T] real los;
  
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}

parameters {
  // Group-level parameters for controls
  vector[4] mu_par_c;
  vector<lower=0>[4] sigma_c;
  
  // Group-level parameters for patients
  vector[4] mu_par_p;
  vector<lower=0>[4] sigma_p;
  
  // Subject-level raw parameters
  vector[N] r_pr;
  vector[N] p_pr;
  vector[N] d_pr;
  vector[N] i_pr;
}

transformed parameters {
  vector<lower=0, upper=1>[N] r;
  vector<lower=0, upper=1>[N] p;
  vector<lower=0, upper=5>[N] d;
  vector<lower=0, upper=5>[N] i;
  
  // Transform subject-level raw parameters using group-specific parameters
  for (j in 1:N) {
    if (group[j] == 1) { // Control
      r[j] = Phi_approx(mu_par_c[1] + sigma_c[1] * r_pr[j]);
      p[j] = Phi_approx(mu_par_c[2] + sigma_c[2] * p_pr[j]);
      d[j] = Phi_approx(mu_par_c[3] + sigma_c[3] * d_pr[j]) * 5;
      i[j] = Phi_approx(mu_par_c[4] + sigma_c[4] * i_pr[j]) * 5;
    } else { // Patient
      r[j] = Phi_approx(mu_par_p[1] + sigma_p[1] * r_pr[j]);
      p[j] = Phi_approx(mu_par_p[2] + sigma_p[2] * p_pr[j]);
      d[j] = Phi_approx(mu_par_p[3] + sigma_p[3] * d_pr[j]) * 5;
      i[j] = Phi_approx(mu_par_p[4] + sigma_p[4] * i_pr[j]) * 5;
    }
  }
}

model {
  // Priors for control group
  mu_par_c ~ normal(0, 2.0);
  sigma_c ~ cauchy(0, 5.0);
  
  // Priors for patient group
  mu_par_p ~ normal(0, 2.0);
  sigma_p ~ cauchy(0, 5.0);
  
  // Subject-level parameters
  r_pr ~ normal(0, 2.0);
  p_pr ~ normal(0, 2.0);
  d_pr ~ normal(0, 2.0);
  i_pr ~ normal(0, 2.0);
  
  // Likelihood
  for (j in 1 : N) {
    if (holdout[j] == 0) {
      // consider subject only if it is not in the holdout data set
      // Define values
      vector[4] Pr; // Resulting response probailities
      vector[3] A; //rule level A values; 1=color; 2=form; 3=number
      vector[3] S;
      
      real sumA;
      real sumS;
      
      // Initialize values
      Pr = rep_vector(0.0, 4);
      A = rep_vector(0.33, 3);
      S = rep_vector(0.0, 3);
      
      // trial 1
      // get response probabilities
      sumA = pow(A[1], d[j]) + pow(A[2], d[j]) + pow(A[3], d[j]);
      Pr[resp_color[j, 1]] = Pr[resp_color[j, 1]] + pow(A[1], d[j]);
      Pr[resp_shape[j, 1]] = Pr[resp_shape[j, 1]] + pow(A[2], d[j]);
      Pr[resp_number[j, 1]] = Pr[resp_number[j, 1]] + pow(A[3], d[j]);
      for (k in 1 : 4) {
        Pr[k] = Pr[k] / sumA;
      }
      
      resp_choice[j, 1] ~ categorical(Pr);
      
      for (t in 1:Tsubj[j]) {
        if (rew[j, t - 1] == 1) {
          S = rep_vector(0.0, 3);
          S[rule_choice[j, t - 1]] = 1;
          
          sumS = S[1] * pow(A[1], i[j]) + S[2] * pow(A[2], i[j])
                 + S[3] * pow(A[3], i[j]);
          S[1] = S[1] * pow(A[1], i[j]) / sumS;
          S[2] = S[2] * pow(A[2], i[j]) / sumS;
          S[3] = S[3] * pow(A[3], i[j]) / sumS;
          
          A = (1 - r[j]) * A + r[j] * S;
        } else {
          S = rep_vector(1.0, 3);
          S[rule_choice[j, t - 1]] = 0;
          
          sumS = S[1] * pow(A[1], i[j]) + S[2] * pow(A[2], i[j])
                 + S[3] * pow(A[3], i[j]);
          S[1] = S[1] * pow(A[1], i[j]) / sumS;
          S[2] = S[2] * pow(A[2], i[j]) / sumS;
          S[3] = S[3] * pow(A[3], i[j]) / sumS;
          
          A = (1 - p[j]) * A + p[j] * S;
        }
        
        // Response probabilities
        Pr = rep_vector(0.0, 4);
        sumA = pow(A[1], d[j]) + pow(A[2], d[j]) + pow(A[3], d[j]);
        Pr[resp_color[j, t]] = Pr[resp_color[j, t]] + pow(A[1], d[j]);
        Pr[resp_shape[j, t]] = Pr[resp_shape[j, t]] + pow(A[2], d[j]);
        Pr[resp_number[j, t]] = Pr[resp_number[j, t]] + pow(A[3], d[j]);
        for (k in 1 : 4) {
          Pr[k] = Pr[k] / sumA;
        }
        
        resp_choice[j, t] ~ categorical(Pr);
      }
    }
  }
}

generated quantities {
  real<lower=0, upper=1> mu_r_c = Phi_approx(mu_par_c[1]);
  real<lower=0, upper=1> mu_p_c = Phi_approx(mu_par_c[2]);
  real<lower=0, upper=5> mu_d_c = Phi_approx(mu_par_c[3]) * 5;
  real<lower=0, upper=5> mu_i_c = Phi_approx(mu_par_c[4]) * 5;
  
  real<lower=0, upper=1> mu_r_p = Phi_approx(mu_par_p[1]);
  real<lower=0, upper=1> mu_p_p = Phi_approx(mu_par_p[2]);
  real<lower=0, upper=5> mu_d_p = Phi_approx(mu_par_p[3]) * 5;
  real<lower=0, upper=5> mu_i_p = Phi_approx(mu_par_p[4]) * 5;
  
  // For log likelihood calculation
  array[N] real log_lik;
  
  // For posterior predictive check
  array[N, T] real y_pred;
  
  for (j in 1:N) {
    log_lik[j] = 0;
    for (t in 1:Tsubj[j]) {
      y_pred[j, t] = -1; // Initialize with a default value or use actual data
    }
  }
  
  mu_r_c = Phi_approx(mu_par_c[1]);
  mu_p_c = Phi_approx(mu_par_c[2]);
  mu_d_c = Phi_approx(mu_par_c[3]) * 5;
  mu_i_c = Phi_approx(mu_par_c[4]) * 5;
  
  mu_r_p = Phi_approx(mu_par_p[1]);
  mu_p_p = Phi_approx(mu_par_p[2]);
  mu_d_p = Phi_approx(mu_par_p[3]) * 5;
  mu_i_p = Phi_approx(mu_par_p[4]) * 5;
  
  {
    // local section, this saves time and space
    for (j in 1 : N) {
      // Define values
      vector[4] Pr; // Resulting response probailities
      vector[3] A; //rule level A values; 1=color; 2=form; 3=number
      vector[3] S;
      
      real sumA;
      real sumS;
      
      // Initialize values
      Pr = rep_vector(0.0, 4);
      A = rep_vector(0.33, 3);
      S = rep_vector(0.0, 3);
      log_lik[j] = 0.0;
      
      // trial 1
      // get response probabilities
      sumA = pow(A[1], d[j]) + pow(A[2], d[j]) + pow(A[3], d[j]);
      Pr[resp_color[j, 1]] = Pr[resp_color[j, 1]] + pow(A[1], d[j]);
      Pr[resp_shape[j, 1]] = Pr[resp_shape[j, 1]] + pow(A[2], d[j]);
      Pr[resp_number[j, 1]] = Pr[resp_number[j, 1]] + pow(A[3], d[j]);
      for (k in 1 : 4) {
        Pr[k] = Pr[k] / sumA;
      }
      
      // compute log likelihood of current trial // only if not odd choice
      log_lik[j] = log_lik[j] + categorical_lpmf(resp_choice[j, 1] | Pr);
      
      // generate posterior prediction for current trial
      y_pred[j, 1] = categorical_rng(Pr);
      
      for (t in 2 : Tsubj[j]) {
        if (rew[j, t - 1] == 1) {
          S = rep_vector(0.0, 3);
          S[rule_choice[j, t - 1]] = 1;
          
          sumS = S[1] * pow(A[1], i[j]) + S[2] * pow(A[2], i[j])
                 + S[3] * pow(A[3], i[j]);
          S[1] = S[1] * pow(A[1], i[j]) / sumS;
          S[2] = S[2] * pow(A[2], i[j]) / sumS;
          S[3] = S[3] * pow(A[3], i[j]) / sumS;
          
          A = (1 - r[j]) * A + r[j] * S;
        } else {
          S = rep_vector(1.0, 3);
          S[rule_choice[j, t - 1]] = 0;
          
          sumS = S[1] * pow(A[1], i[j]) + S[2] * pow(A[2], i[j])
                 + S[3] * pow(A[3], i[j]);
          S[1] = S[1] * pow(A[1], i[j]) / sumS;
          S[2] = S[2] * pow(A[2], i[j]) / sumS;
          S[3] = S[3] * pow(A[3], i[j]) / sumS;
          
          A = (1 - p[j]) * A + p[j] * S;
        }
        
        // response probabilities
        Pr = rep_vector(0.0, 4);
        sumA = pow(A[1], d[j]) + pow(A[2], d[j]) + pow(A[3], d[j]);
        Pr[resp_color[j, t]] = Pr[resp_color[j, t]] + pow(A[1], d[j]);
        Pr[resp_shape[j, t]] = Pr[resp_shape[j, t]] + pow(A[2], d[j]);
        Pr[resp_number[j, t]] = Pr[resp_number[j, t]] + pow(A[3], d[j]);
        for (k in 1 : 4) {
          Pr[k] = Pr[k] / sumA;
        }
        
        // compute log likelihood of current trial
        log_lik[j] = log_lik[j] + categorical_lpmf(resp_choice[j, t] | Pr);
        
        // generate posterior prediction for current trial
        y_pred[j, t] = categorical_rng(Pr);
      }
    }
  }
}
