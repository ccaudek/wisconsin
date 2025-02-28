data {
  int<lower=1> N;
  int<lower=1> T;
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
transformed data {
  int K = 6;
  int M = 10; // number of mixture components
}
parameters {
  real<lower=0> alpha_dp;
  vector<lower=0, upper=1>[M] v;
  
  array[M] vector[K] mu_cluster;
  array[M] vector<lower=0>[K] sigma_cluster;
  
  array[N] vector[K] theta; // subject-level parameters (no non-centered now)
}
transformed parameters {
  vector[M] w_m;
  {
    real prod = 1.0;
    for (m in 1 : (M - 1)) {
      w_m[m] = v[m] * prod;
      prod *= (1 - v[m]);
    }
    w_m[M] = prod;
  }
}
model {
  // Priors DP
  alpha_dp ~ gamma(2, 1);
  
  for (m in 1 : (M - 1)) 
    v[m] ~ beta(1, alpha_dp);
  v[M] ~ beta(1, alpha_dp); // though redundant
  
  // Priors cluster parameters
  for (m in 1 : M) {
    for (k in 1 : K) {
      mu_cluster[m][k] ~ normal(0, 2);
      sigma_cluster[m][k] ~ cauchy(0, 2);
    }
  }
  
  // Prior on subjects via mixture
  // p(theta[i]) = sum_m w_m * prod_k normal(theta[i,k]| mu_cluster[m,k], sigma_cluster[m,k])
  // We must write the likelihood in a marginalized form:
  // log p(theta[i]) = log( sum_m w_m * prod_k normal(...))
  // Use log_sum_exp trick or log_mix in Stan. 
  // We'll loop over subjects:
  for (i in 1 : N) {
    vector[M] lps; // log probabilities for each mixture component
    for (m in 1 : M) {
      real lp_component = 0;
      for (k in 1 : K) {
        lp_component += normal_lpdf(theta[i, k] | mu_cluster[m][k], sigma_cluster[m][k]);
      }
      lps[m] = log(w_m[m]) + lp_component;
    }
    target += log_sum_exp(lps);
  }
  
  // Now the Q-learning likelihood given the parameters theta:
  // Transform parameters
  // alpha_rew = inv_logit(theta[i,1])
  // alpha_pun = inv_logit(theta[i,2])
  // gamma_ = inv_logit(theta[i,3])
  // tau = exp(theta[i,4])
  // w = inv_logit(theta[i,5])
  // lapse = inv_logit(theta[i,6])
  
  // Add the behavioral likelihood
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      real alpha_rew_i = inv_logit(theta[i, 1]);
      real alpha_pun_i = inv_logit(theta[i, 2]);
      real gamma_i = inv_logit(theta[i, 3]);
      real tau_i = exp(theta[i, 4]);
      real w_i = inv_logit(theta[i, 5]);
      real lapse_i = inv_logit(theta[i, 6]);
      
      vector[4] Q = rep_vector(0.0, 4);
      
      for (t in 1 : Tsubj[i]) {
        vector[4] Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        
        Qeff[resp_color[i, t]] += w_i * Q[resp_color[i, t]];
        Qeff[resp_shape[i, t]] += w_i * Q[resp_shape[i, t]];
        Qeff[resp_number[i, t]] += w_i * Q[resp_number[i, t]];
        
        vector[4] p_choice = softmax(Qeff / tau_i);
        p_choice = (1 - lapse_i) * p_choice
                   + lapse_i * (rep_vector(1.0 / 4.0, 4));
        
        resp_choice[i, t] ~ categorical(p_choice);
        
        Q = gamma_i * Q;
        real outcome = rew[i, t] - los[i, t];
        int c = resp_choice[i, t];
        real PE = outcome - Q[c];
        
        if (rew[i, t] == 1) {
          Q[c] = Q[c] + alpha_rew_i * PE;
        } else if (los[i, t] == 1) {
          Q[c] = Q[c] + alpha_pun_i * PE;
        } else {
          Q[c] = Q[c] + 0.5 * (alpha_rew_i + alpha_pun_i) * PE;
        }
      }
    }
  }
}
generated quantities {
  array[N] real alpha_rew_i;
  array[N] real alpha_pun_i;
  array[N] real gamma_i;
  array[N] real tau_i;
  array[N] real w_i;
  array[N] real lapse_i;
  
  array[N] real log_lik;
  array[N, T] int y_pred;
  
  for (i in 1 : N) {
    // Calcola i parametri soggetto da theta
    alpha_rew_i[i] = inv_logit(theta[i, 1]);
    alpha_pun_i[i] = inv_logit(theta[i, 2]);
    gamma_i[i] = inv_logit(theta[i, 3]);
    tau_i[i] = exp(theta[i, 4]);
    w_i[i] = inv_logit(theta[i, 5]);
    lapse_i[i] = inv_logit(theta[i, 6]);
    
    vector[4] Q = rep_vector(0.0, 4);
    log_lik[i] = 0;
    
    for (t in 1 : Tsubj[i]) {
      vector[4] Qeff = rep_vector(-1.0, 4);
      Qeff[resp_color[i, t]] = 0.0;
      Qeff[resp_shape[i, t]] = 0.0;
      Qeff[resp_number[i, t]] = 0.0;
      
      Qeff[resp_color[i, t]] += w_i[i] * Q[resp_color[i, t]];
      Qeff[resp_shape[i, t]] += w_i[i] * Q[resp_shape[i, t]];
      Qeff[resp_number[i, t]] += w_i[i] * Q[resp_number[i, t]];
      
      vector[4] p_choice = softmax(Qeff / tau_i[i]);
      p_choice = (1 - lapse_i[i]) * p_choice
                 + lapse_i[i] * (rep_vector(1.0 / 4.0, 4));
      
      log_lik[i] += categorical_lpmf(resp_choice[i, t] | p_choice);
      y_pred[i, t] = categorical_rng(p_choice);
      
      Q = gamma_i[i] * Q;
      real outcome = rew[i, t] - los[i, t];
      int c = resp_choice[i, t];
      real PE = outcome - Q[c];
      
      if (rew[i, t] == 1) {
        Q[c] = Q[c] + alpha_rew_i[i] * PE;
      } else if (los[i, t] == 1) {
        Q[c] = Q[c] + alpha_pun_i[i] * PE;
      } else {
        Q[c] = Q[c] + 0.5 * (alpha_rew_i[i] + alpha_pun_i[i]) * PE;
      }
    }
    for (t in (Tsubj[i] + 1) : T) {
      y_pred[i, t] = -1;
    }
  }
}
