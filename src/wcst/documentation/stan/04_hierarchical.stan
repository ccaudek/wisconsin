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

parameters {
  // GROUP-LEVEL: prior molto vaghi
  real<lower=0, upper=1> mu_alpha;
  real<lower=0> mu_beta;
  real<lower=0> mu_inertia;
  
  // CRUCIALE: sigma più grandi = più variabilità tra soggetti
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_inertia;
  
  // SUBJECT-LEVEL: non-centered parameterization
  vector[N] alpha_raw;
  vector[N] beta_raw;
  vector[N] inertia_raw;
}

transformed parameters {
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0>[N] beta;
  vector<lower=0>[N] inertia;
  
  for (i in 1:N) {
    alpha[i] = inv_logit(logit(mu_alpha) + sigma_alpha * alpha_raw[i]);
    beta[i] = exp(log(mu_beta) + sigma_beta * beta_raw[i]);
    inertia[i] = exp(log(mu_inertia) + sigma_inertia * inertia_raw[i]);
  }
}

model {
  // PRIOR MOLTO VAGHI per le medie di gruppo
  mu_alpha ~ beta(1, 1);  // uniform su (0,1)
  mu_beta ~ gamma(1, 0.5);  // mean=2, sd=2.83
  mu_inertia ~ gamma(1, 0.5);  // mean=2, sd=2.83
  
  // PRIOR VAGHI per le SD: permettono molta variabilità
  sigma_alpha ~ normal(0, 2);  // molto più largo!
  sigma_beta ~ normal(0, 3);   // molto più largo!
  sigma_inertia ~ normal(0, 2);  // molto più largo!
  
  // Standard normal per raw parameters
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();
  inertia_raw ~ std_normal();
  
  // LIKELIHOOD
  for (i in 1:N) {
    if (holdout[i] == 0) {
      vector[3] Q = rep_vector(0.0, 3);
      int last_dim = 0;
      
      for (t in 1:Tsubj[i]) {
        if (resp_choice[i, t] > 0) {
          vector[3] dim_logits = beta[i] * Q;
          if (last_dim > 0) {
            dim_logits[last_dim] += inertia[i];
          }
          vector[3] dim_probs = softmax(dim_logits);
          
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) {
            chosen_dim = 1;
          } else if (resp_choice[i, t] == resp_shape[i, t]) {
            chosen_dim = 2;
          } else if (resp_choice[i, t] == resp_number[i, t]) {
            chosen_dim = 3;
          }
          
          if (chosen_dim > 0) {
            target += log(dim_probs[chosen_dim]);
            real outcome = rew[i, t] + los[i, t];
            Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
            last_dim = chosen_dim;
          }
        }
      }
    }
  }
}

generated quantities {
  array[N] real log_lik;
  array[N, T] int y_pred;
  
  for (i in 1:N) {
    log_lik[i] = 0;
    vector[3] Q = rep_vector(0.0, 3);
    int last_dim = 0;
    
    for (t in 1:Tsubj[i]) {
      y_pred[i, t] = -1;
      
      if (resp_choice[i, t] > 0) {
        vector[3] dim_logits = beta[i] * Q;
        if (last_dim > 0) {
          dim_logits[last_dim] += inertia[i];
        }
        vector[3] dim_probs = softmax(dim_logits);
        
        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) {
          chosen_dim = 1;
        } else if (resp_choice[i, t] == resp_shape[i, t]) {
          chosen_dim = 2;
        } else if (resp_choice[i, t] == resp_number[i, t]) {
          chosen_dim = 3;
        }
        
        if (chosen_dim > 0) {
          log_lik[i] += log(dim_probs[chosen_dim]);
          
          int sim_dim = categorical_rng(dim_probs);
          if (sim_dim == 1) {
            y_pred[i, t] = resp_color[i, t];
          } else if (sim_dim == 2) {
            y_pred[i, t] = resp_shape[i, t];
          } else {
            y_pred[i, t] = resp_number[i, t];
          }
          
          real outcome = rew[i, t] + los[i, t];
          Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
          last_dim = chosen_dim;
        }
      }
    }
  }
}
