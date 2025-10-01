data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Max number of trials
  array[N] int<lower=1, upper=T> Tsubj; // Trials per subject
  array[N] int<lower=0, upper=1> holdout; // Holdout indicator
  
  array[N, T] real rew; // Reward (e.g. +1 for correct)
  array[N, T] real los; // Loss   (e.g. -1 for incorrect)
  
  // Subject's rule choice and response:
  //   rule_choice is the correct rule dimension (1..3),
  //   resp_choice is the chosen card (1..4, or 0 if missing)
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  
  // Each dimension’s “matching card” on that trial
  array[N, T] int resp_color; // The card index that would match color
  array[N, T] int resp_shape; // The card index that would match shape
  array[N, T] int resp_number; // The card index that would match number
}
parameters {
  // ---------------------- GROUP-LEVEL MEANS & SDs ----------------------
  // alpha parameters
  real<lower=0, upper=1> mu_alpha; // learning rate mean
  real<lower=0> mu_beta; // inverse temperature mean
  real<lower=0> mu_inertia; // inertia mean (can be wide if you want)
  
  real<lower=0, upper=1> sigma_alpha; // SD of alpha across subjects
  real<lower=0> sigma_beta; // SD of beta across subjects
  real<lower=0> sigma_inertia; // SD of inertia across subjects
  
  // ---------------------- SUBJECT-LEVEL "RAW" PARAMETERS ----------------------
  vector[N] alpha_raw; // transforms to alpha in (0,1)
  vector[N] beta_raw; // transforms to beta > 0
  vector[N] inertia_raw; // transforms to inertia > 0
}
transformed parameters {
  // ---------------------- SUBJECT-LEVEL TRANSFORMED PARAMETERS ----------------------
  vector<lower=0, upper=1>[N] alpha; // learning rate in (0,1)
  vector<lower=0>[N] beta; // inverse temperature > 0
  vector<lower=0>[N] inertia; // “stickiness” or “inertia” > 0
  
  for (i in 1 : N) {
    // alpha in (0,1): logistic transform around mu_alpha
    alpha[i] = inv_logit(logit(mu_alpha) + sigma_alpha * alpha_raw[i]);
    
    // beta > 0: exponential transform around mu_beta
    beta[i] = exp(log(mu_beta) + sigma_beta * beta_raw[i]);
    
    // inertia > 0: exponential transform around mu_inertia
    inertia[i] = exp(log(mu_inertia) + sigma_inertia * inertia_raw[i]);
  }
}
model {
  // ---------------------- PRIORS ----------------------
  // Group-level priors
  mu_alpha ~ beta(2, 2);
  mu_beta ~ gamma(2, 1);        // CAMBIATO
  mu_inertia ~ gamma(1.5, 1);   // CAMBIATO
  
  sigma_alpha ~ exponential(2);  // CAMBIATO
  sigma_beta ~ exponential(1);   // CAMBIATO
  sigma_inertia ~ exponential(1); // CAMBIATO
  
  // Subject-level priors
  alpha_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  inertia_raw ~ normal(0, 1);
  
  // ---------------------- LIKELIHOOD ----------------------
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      // only use non-holdout data for log-likelihood
      // Initialize Q-values for the 3 dimensions
      vector[3] Q = rep_vector(0.0, 3);
      
      // We'll keep track of which dimension was chosen last
      int last_dim = 0; // 0 means "none" for the first trial
      
      for (t in 1 : Tsubj[i]) {
        if (resp_choice[i, t] > 0) {
          // Construct logits for each dimension = beta * Q
          vector[3] dim_logits = beta[i] * Q;
          
          // Add inertia to the dimension chosen on the previous trial
          if (last_dim > 0) {
            dim_logits[last_dim] += inertia[i];
          }
          
          // Convert logits to probabilities
          vector[3] dim_probs = softmax(dim_logits);
          
          // Figure out which dimension was chosen
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) {
            chosen_dim = 1;
          } else if (resp_choice[i, t] == resp_shape[i, t]) {
            chosen_dim = 2;
          } else if (resp_choice[i, t] == resp_number[i, t]) {
            chosen_dim = 3;
          }
          
          // If we can map the card to a dimension, update the log-likelihood
          if (chosen_dim > 0) {
            target += log(dim_probs[chosen_dim]);
            
            // RL update
            real outcome = rew[i, t] + los[i, t]; // e.g., +1 correct, -1 incorrect
            Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
            
            // Update last_dim
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
  
  for (i in 1 : N) {
    // Subject-level log_lik
    log_lik[i] = 0;
    
    // Re-run the same logic for posterior predictive
    vector[3] Q = rep_vector(0.0, 3);
    int last_dim = 0;
    
    for (t in 1 : Tsubj[i]) {
      y_pred[i, t] = -1; // placeholder
      
      if (resp_choice[i, t] > 0) {
        // Construct logits
        vector[3] dim_logits = beta[i] * Q;
        if (last_dim > 0) {
          dim_logits[last_dim] += inertia[i];
        }
        vector[3] dim_probs = softmax(dim_logits);
        
        // Which dimension was actually chosen by the subject?
        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) {
          chosen_dim = 1;
        } else if (resp_choice[i, t] == resp_shape[i, t]) {
          chosen_dim = 2;
        } else if (resp_choice[i, t] == resp_number[i, t]) {
          chosen_dim = 3;
        }
        
        // Log-likelihood
        if (chosen_dim > 0) {
          log_lik[i] += log(dim_probs[chosen_dim]);
          
          // Posterior predictive: draw dimension from dim_probs
          int sim_dim = categorical_rng(dim_probs);
          
          // Convert dimension to a card index
          if (sim_dim == 1) {
            y_pred[i, t] = resp_color[i, t];
          } else if (sim_dim == 2) {
            y_pred[i, t] = resp_shape[i, t];
          } else {
            y_pred[i, t] = resp_number[i, t];
          }
          
          // Update Q with the outcome
          real outcome = rew[i, t] + los[i, t];
          Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
          
          last_dim = chosen_dim;
        }
      }
    }
  }
}
