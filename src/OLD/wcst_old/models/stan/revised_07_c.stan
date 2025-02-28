data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Number of trials
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
  int K = 6; // alpha_rew, alpha_pun, gamma, tau, w, lapse
}
parameters {
  vector[K] mu_p;
  vector<lower=0>[K] sigma_param;
  
  array[N] vector[K] theta_raw;
}
transformed parameters {
  array[N] real alpha_rew;
  array[N] real alpha_pun;
  array[N] real gamma_;
  array[N] real tau;
  array[N] real w;
  array[N] real lapse;
  
  for (i in 1 : N) {
    vector[K] subj_params = mu_p + sigma_param .* theta_raw[i];
    
    // Nessun clamping, lasciamo libertà ai parametri
    alpha_rew[i] = inv_logit(subj_params[1]);
    alpha_pun[i] = inv_logit(subj_params[2]);
    gamma_[i] = inv_logit(subj_params[3]);
    tau[i] = exp(subj_params[4]);
    w[i] = inv_logit(subj_params[5]);
    lapse[i] = inv_logit(subj_params[6]);
  }
}
model {
  // Priors più larghi per mu_p
  mu_p ~ normal(0, 2);
  
  // Priors half-Cauchy per sigma_param
  // In Stan, una mezza-Cauchy si ottiene con cauchy(0,2) e constraint <lower=0>
  sigma_param ~ cauchy(0, 2);
  
  for (i in 1 : N) {
    theta_raw[i] ~ normal(0, 1);
  }
  
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      vector[4] Q = rep_vector(0.0, 4);
      for (t in 1 : Tsubj[i]) {
        vector[4] Qeff = rep_vector(-1.0, 4);
        Qeff[resp_color[i, t]] = 0.0;
        Qeff[resp_shape[i, t]] = 0.0;
        Qeff[resp_number[i, t]] = 0.0;
        
        Qeff[resp_color[i, t]] += w[i] * Q[resp_color[i, t]];
        Qeff[resp_shape[i, t]] += w[i] * Q[resp_shape[i, t]];
        Qeff[resp_number[i, t]] += w[i] * Q[resp_number[i, t]];
        
        vector[4] p_choice = softmax(Qeff / tau[i]);
        p_choice = (1 - lapse[i]) * p_choice
                   + lapse[i] * (rep_vector(1.0 / 4.0, 4));
        
        resp_choice[i, t] ~ categorical(p_choice);
        
        Q = gamma_[i] * Q;
        
        real outcome = rew[i, t] - los[i, t];
        int c = resp_choice[i, t];
        real PE = outcome - Q[c];
        
        if (rew[i, t] == 1) {
          Q[c] = Q[c] + alpha_rew[i] * PE;
        } else if (los[i, t] == 1) {
          Q[c] = Q[c] + alpha_pun[i] * PE;
        } else {
          Q[c] = Q[c] + 0.5 * (alpha_rew[i] + alpha_pun[i]) * PE;
        }
      }
    }
  }
}
generated quantities {
  array[N] real log_lik;
  array[N, T] int y_pred;
  
  for (i in 1 : N) {
    vector[4] Q = rep_vector(0.0, 4);
    log_lik[i] = 0;
    for (t in 1 : Tsubj[i]) {
      vector[4] Qeff = rep_vector(-1.0, 4);
      Qeff[resp_color[i, t]] = 0.0;
      Qeff[resp_shape[i, t]] = 0.0;
      Qeff[resp_number[i, t]] = 0.0;
      
      Qeff[resp_color[i, t]] += w[i] * Q[resp_color[i, t]];
      Qeff[resp_shape[i, t]] += w[i] * Q[resp_shape[i, t]];
      Qeff[resp_number[i, t]] += w[i] * Q[resp_number[i, t]];
      
      vector[4] p_choice = softmax(Qeff / tau[i]);
      p_choice = (1 - lapse[i]) * p_choice
                 + lapse[i] * (rep_vector(1.0 / 4.0, 4));
      log_lik[i] += categorical_lpmf(resp_choice[i, t] | p_choice);
      y_pred[i, t] = categorical_rng(p_choice);
      
      Q = gamma_[i] * Q;
      real outcome = rew[i, t] - los[i, t];
      int c = resp_choice[i, t];
      real PE = outcome - Q[c];
      
      if (rew[i, t] == 1) {
        Q[c] = Q[c] + alpha_rew[i] * PE;
      } else if (los[i, t] == 1) {
        Q[c] = Q[c] + alpha_pun[i] * PE;
      } else {
        Q[c] = Q[c] + 0.5 * (alpha_rew[i] + alpha_pun[i]) * PE;
      }
    }
    for (t in (Tsubj[i] + 1) : T) {
      y_pred[i, t] = -1;
    }
  }
}
