data {
  int<lower=1> T;  // number of trials for this subject
  array[T] int<lower=-1, upper=1> rew;
  array[T] int<lower=-1, upper=1> los;
  array[T] int resp_choice;
  array[T] int resp_color;
  array[T] int resp_shape;
  array[T] int resp_number;
}

parameters {
  real<lower=1e-4, upper=1-1e-4> alpha_pos;  // learning rate for positive feedback
  real<lower=1e-4, upper=1-1e-4> alpha_neg;  // learning rate for negative feedback
  real<lower=0> beta;                // inverse temperature (softmax)
}

model {
  // Priors
  alpha_pos ~ beta(2, 2);
  alpha_neg ~ beta(2, 2);
  beta ~ gamma(2, 1);
  
  // Initialize Q-values for 3 dimensions
  vector[3] Q = rep_vector(0.0, 3);
  
  for (t in 1:T) {
    if (resp_choice[t] > 0) {
      
      // Compute dimension probabilities via softmax
      vector[3] dim_probs = softmax(beta * Q);
      
      // Identify chosen dimension
      int chosen_dim = 0;
      if (resp_choice[t] == resp_color[t]) {
        chosen_dim = 1;
      } else if (resp_choice[t] == resp_shape[t]) {
        chosen_dim = 2;
      } else if (resp_choice[t] == resp_number[t]) {
        chosen_dim = 3;
      }
      
      // Update likelihood and Q-values
      if (chosen_dim > 0) {
        target += log(dim_probs[chosen_dim]);
        
        real outcome = rew[t] + los[t];
        
        // Use different learning rates for positive vs negative outcomes
        if (outcome > 0) {
          Q[chosen_dim] += alpha_pos * (outcome - Q[chosen_dim]);
        } else {
          Q[chosen_dim] += alpha_neg * (outcome - Q[chosen_dim]);
        }
      }
    }
  }
}

generated quantities {
  real log_lik = 0;
  array[T] int y_pred;
  
  vector[3] Q = rep_vector(0.0, 3);
  
  for (t in 1:T) {
    y_pred[t] = -1;
    
    if (resp_choice[t] > 0) {
      vector[3] dim_probs = softmax(beta * Q);
      
      int chosen_dim = 0;
      if (resp_choice[t] == resp_color[t]) {
        chosen_dim = 1;
      } else if (resp_choice[t] == resp_shape[t]) {
        chosen_dim = 2;
      } else if (resp_choice[t] == resp_number[t]) {
        chosen_dim = 3;
      }
      
      if (chosen_dim > 0) {
        log_lik += log(dim_probs[chosen_dim]);
        
        int sim_dim = categorical_rng(dim_probs);
        if (sim_dim == 1) {
          y_pred[t] = resp_color[t];
        } else if (sim_dim == 2) {
          y_pred[t] = resp_shape[t];
        } else {
          y_pred[t] = resp_number[t];
        }
        
        real outcome = rew[t] + los[t];
        if (outcome > 0) {
          Q[chosen_dim] += alpha_pos * (outcome - Q[chosen_dim]);
        } else {
          Q[chosen_dim] += alpha_neg * (outcome - Q[chosen_dim]);
        }
      }
    }
  }
}
