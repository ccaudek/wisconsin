data {
  int<lower=1> N; // Number of subjects
  int<lower=1> T; // Max number of trials
  array[N] int<lower=1, upper=T> Tsubj; // Trials per subject
  array[N] int<lower=0, upper=1> holdout; // Holdout indicator
  
  array[N, T] real rew; // Reward (1 for correct)
  array[N, T] real los; // Loss  (-1 for incorrect)
  
  array[N, T] int rule_choice; // Correct rule: 1=color, 2=shape, 3=number
  array[N, T] int resp_choice; // Subject's card choice: 1..4 (0 if missing)
  
  // Each dimension's "matching card" on that trial
  array[N, T] int resp_color; // The card index that would match color
  array[N, T] int resp_shape; // The card index that would match shape
  array[N, T] int resp_number; // The card index that would match number
}
parameters {
  // Group-level means and SDs for alpha, beta
  real<lower=0, upper=1> mu_alpha;
  real<lower=0> mu_beta;
  
  real<lower=0, upper=1> sigma_alpha; // SD for alpha
  real<lower=0> sigma_beta; // SD for beta
  
  // Subject-level parameters (raw scale, will transform)
  vector[N] alpha_raw;
  vector[N] beta_raw;
}
transformed parameters {
  // Transform subject-level parameters
  vector<lower=0, upper=1>[N] alpha; // learning rate
  vector<lower=0>[N] beta; // inverse temperature (softmax)
  
  for (i in 1 : N) {
    // alpha in (0,1), via a logit transform around mu_alpha
    // alpha_raw[i] ~ Normal(0,1) then scaled by sigma_alpha
    alpha[i] = inv_logit(logit(mu_alpha) + sigma_alpha * alpha_raw[i]);
    
    // beta > 0, via a log transform around mu_beta
    // beta_raw[i] ~ Normal(0,1) then scaled by sigma_beta
    beta[i] = exp(log(mu_beta) + sigma_beta * beta_raw[i]);
  }
}
model {
  // ------------------- PRIORS -------------------
  // Group-level priors
  mu_alpha ~ beta(2, 2); // learning rate typically around 0.5
  mu_beta ~ normal(0, 2); // broad prior for inverse temperature
  sigma_alpha ~ normal(0, 0.5); // half-normal
  sigma_beta ~ normal(0, 0.5); // half-normal
  
  // Subject-level priors
  alpha_raw ~ normal(0, 1);
  beta_raw ~ normal(0, 1);
  
  // ------------------- LIKELIHOOD -------------------
  for (i in 1 : N) {
    // Skip holdout data for log-likelihood
    if (holdout[i] == 0) {
      // Initialize Q-values for the 3 dimensions: color/shape/number
      vector[3] Q = rep_vector(0.0, 3);
      
      for (t in 1 : Tsubj[i]) {
        // Exclude missing choices
        if (resp_choice[i, t] > 0) {
          // 1) Compute dimension probabilities via softmax
          vector[3] dim_probs = softmax(beta[i] * Q);
          
          // 2) Identify the dimension chosen by the subject
          //    (based on which "matching card" is actually chosen)
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) {
            chosen_dim = 1; // color
          } else if (resp_choice[i, t] == resp_shape[i, t]) {
            chosen_dim = 2; // shape
          } else if (resp_choice[i, t] == resp_number[i, t]) {
            chosen_dim = 3; // number
          }
          
          // Only update likelihood if we can map the card to a dimension
          if (chosen_dim > 0) {
            // Probability of selecting that dimension
            target += log(dim_probs[chosen_dim]);
            
            // 3) RL Update for the chosen dimension
            real outcome = rew[i, t] + los[i, t]; // e.g., +1 correct, -1 incorrect
            Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
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
    // We'll store a subject-level log_lik for LOO or WAIC
    log_lik[i] = 0;
    
    // Re-run the same logic, but keep track of the probability
    // and generate predictions
    vector[3] Q = rep_vector(0.0, 3);
    
    for (t in 1 : Tsubj[i]) {
      y_pred[i, t] = -1; // default placeholder
      
      if (resp_choice[i, t] > 0) {
        vector[3] dim_probs = softmax(beta[i] * Q);
        
        // Which dimension was chosen by the subject?
        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) {
          chosen_dim = 1;
        } else if (resp_choice[i, t] == resp_shape[i, t]) {
          chosen_dim = 2;
        } else if (resp_choice[i, t] == resp_number[i, t]) {
          chosen_dim = 3;
        }
        
        // Log-likelihood contribution
        if (chosen_dim > 0) {
          log_lik[i] += log(dim_probs[chosen_dim]);
          
          // Posterior predictive: pick dimension from dim_probs
          int sim_dim = categorical_rng(dim_probs);
          
          // Convert dimension -> a card index
          if (sim_dim == 1) {
            y_pred[i, t] = resp_color[i, t];
          } else if (sim_dim == 2) {
            y_pred[i, t] = resp_shape[i, t];
          } else if (sim_dim == 3) {
            y_pred[i, t] = resp_number[i, t];
          }
          
          // Update Q
          real outcome = rew[i, t] + los[i, t];
          Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
        }
      }
    }
  }
}
