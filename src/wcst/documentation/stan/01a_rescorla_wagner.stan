data {
  int<lower=1> N; // number of subjects
  int<lower=1> T; // max number of trials per subject
  array[N] int<lower=1, upper=T> Tsubj; // actual trials per subject
  array[N] int<lower=0, upper=1> holdout; // holdout indicator: 0 => use for fitting
  
  // trial-level outcomes
  array[N, T] real rew; // e.g. +1 for correct
  array[N, T] real los; // e.g. -1 for incorrect
  
  // participant's rule / response
  array[N, T] int rule_choice; // correct rule (1..3)
  array[N, T] int resp_choice; // chosen card index (1..4, or 0 if missing)
  
  // card indices matching each dimension
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}
parameters {
  // ------------------------------ GROUP-LEVEL ------------------------------
  // Population-level means with wider ranges
  real<lower=0, upper=1> mu_alpha_pos; // learning rate for positive outcomes
  real<lower=0, upper=1> mu_alpha_neg; // learning rate for negative outcomes
  real<lower=0> mu_beta; // inverse temperature
  real<lower=0, upper=5> mu_persev; // perseveration strength
  
  // Population-level SDs - much wider to allow individual differences
  real<lower=0> sigma_alpha_pos;
  real<lower=0> sigma_alpha_neg;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_persev;
  
  // Correlation structure
  cholesky_factor_corr[4] L_Rho;
  
  // ------------------------------ SUBJECT-LEVEL ----------------------------
  // Individual-level parameters (non-centered)
  vector[N] z_alpha_pos;
  vector[N] z_alpha_neg;
  vector[N] z_beta;
  vector[N] z_persev;
}
transformed parameters {
  // Subject-level parameters
  vector<lower=0, upper=1>[N] alpha_pos; // learning rate for positive feedback
  vector<lower=0, upper=1>[N] alpha_neg; // learning rate for negative feedback
  vector<lower=0>[N] beta; // inverse temperature
  vector<lower=0>[N] persev; // perseveration strength
  
  // Transform parameters to appropriate scales 
  for (i in 1 : N) {
    alpha_pos[i] = inv_logit(logit(mu_alpha_pos)
                             + sigma_alpha_pos * z_alpha_pos[i]);
    alpha_neg[i] = inv_logit(logit(mu_alpha_neg)
                             + sigma_alpha_neg * z_alpha_neg[i]);
    beta[i] = exp(log(mu_beta) + sigma_beta * z_beta[i]);
    persev[i] = exp(log(mu_persev) + sigma_persev * z_persev[i]);
  }
}
model {
  // ---------------------- PRIORS ----------------------
  // Population means with minimally informative priors
  mu_alpha_pos ~ beta(1, 1); // uniform prior
  mu_alpha_neg ~ beta(1, 1); // uniform prior
  mu_beta ~ gamma(2, 0.2); // wider prior allowing higher inverse temps
  mu_persev ~ gamma(1, 0.5); // allow moderate perseveration
  
  // Standard deviations with wider priors
  sigma_alpha_pos ~ normal(0, 0.5);
  sigma_alpha_neg ~ normal(0, 0.5);
  sigma_beta ~ normal(0, 5); // much wider prior
  sigma_persev ~ normal(0, 2); // wider prior
  
  // Correlation structure
  L_Rho ~ lkj_corr_cholesky(2.0); // slight regularization toward smaller correlations
  
  // Non-centered parameterization with correlation structure
  {
    // Joint distribution for all z parameters
    matrix[N, 4] z = append_col(append_col(append_col(z_alpha_pos,
                                                      z_alpha_neg),
                                           z_beta),
                                z_persev);
    
    for (i in 1 : N) {
      z[i] ~ multi_normal_cholesky([0, 0, 0, 0]', L_Rho);
    }
  }
  
  // ---------------------- LIKELIHOOD --------------------
  for (i in 1 : N) {
    if (holdout[i] == 0) {
      // Initialize Q-values for 3 dimensions
      vector[3] Q = rep_vector(0.0, 3);
      
      // Initialize previous choice variables
      array[3] int prev_dim_chosen = {0, 0, 0};
      
      for (t in 1 : Tsubj[i]) {
        if (resp_choice[i, t] > 0) {
          // Compute choice values including Q-values and perseveration
          vector[3] choice_values = beta[i] * Q;
          
          // Add perseveration bonus for previously chosen dimensions
          for (d in 1 : 3) {
            if (prev_dim_chosen[d] == 1) {
              choice_values[d] += persev[i];
            }
          }
          
          // Calculate dimension probabilities via softmax
          vector[3] dim_probs = softmax(choice_values);
          
          // Find which dimension was chosen
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) 
            chosen_dim = 1;
          else if (resp_choice[i, t] == resp_shape[i, t]) 
            chosen_dim = 2;
          else if (resp_choice[i, t] == resp_number[i, t]) 
            chosen_dim = 3;
          
          // Update likelihood & Q-values if a valid dimension was chosen
          if (chosen_dim > 0) {
            target += log(dim_probs[chosen_dim]);
            
            // Get outcome
            real outcome = rew[i, t] + los[i, t];
            
            // Update Q-value differently based on positive or negative outcome
            if (outcome > 0) {
              Q[chosen_dim] += alpha_pos[i] * (outcome - Q[chosen_dim]);
            } else {
              Q[chosen_dim] += alpha_neg[i] * (outcome - Q[chosen_dim]);
            }
            
            // Reset perseveration markers and set for chosen dimension
            prev_dim_chosen = {0, 0, 0};
            prev_dim_chosen[chosen_dim] = 1;
          }
        }
      }
    }
  }
}
generated quantities {
  array[N] real log_lik;
  array[N, T] int y_pred;
  matrix[4, 4] Rho = multiply_lower_tri_self_transpose(L_Rho);
  
  // Additional parameter combinations for analysis
  vector[N] alpha_diff = alpha_pos - alpha_neg; // difference between positive/negative learning
  vector[N] alpha_sum = alpha_pos + alpha_neg; // overall learning sensitivity
  
  for (i in 1 : N) {
    log_lik[i] = 0;
    vector[3] Q = rep_vector(0.0, 3);
    array[3] int prev_dim_chosen = {0, 0, 0};
    
    for (t in 1 : Tsubj[i]) {
      y_pred[i, t] = -1; // default placeholder
      
      if (resp_choice[i, t] > 0) {
        // Compute choice values including Q-values and perseveration
        vector[3] choice_values = beta[i] * Q;
        
        // Add perseveration bonus for previously chosen dimensions
        for (d in 1 : 3) {
          if (prev_dim_chosen[d] == 1) {
            choice_values[d] += persev[i];
          }
        }
        
        // Calculate dimension probabilities
        vector[3] dim_probs = softmax(choice_values);
        
        // Find which dimension the subject chose
        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) 
          chosen_dim = 1;
        else if (resp_choice[i, t] == resp_shape[i, t]) 
          chosen_dim = 2;
        else if (resp_choice[i, t] == resp_number[i, t]) 
          chosen_dim = 3;
        
        // Calculate log-likelihood and predicted choices
        if (chosen_dim > 0) {
          log_lik[i] += log(dim_probs[chosen_dim]);
          
          // Generate posterior prediction
          int sim_dim = categorical_rng(dim_probs);
          if (sim_dim == 1) 
            y_pred[i, t] = resp_color[i, t];
          else if (sim_dim == 2) 
            y_pred[i, t] = resp_shape[i, t];
          else if (sim_dim == 3) 
            y_pred[i, t] = resp_number[i, t];
          
          // Update Q-value
          real outcome = rew[i, t] + los[i, t];
          if (outcome > 0) {
            Q[chosen_dim] += alpha_pos[i] * (outcome - Q[chosen_dim]);
          } else {
            Q[chosen_dim] += alpha_neg[i] * (outcome - Q[chosen_dim]);
          }
          
          // Reset perseveration markers and set for chosen dimension
          prev_dim_chosen = {0, 0, 0};
          prev_dim_chosen[chosen_dim] = 1;
        }
      }
    }
  }
}
