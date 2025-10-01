data {
  int<lower=1> T;  // number of trials for this subject
  array[T] int<lower=-1, upper=1> rew;  // +1 for correct
  array[T] int<lower=-1, upper=1> los;  // -1 for incorrect
  array[T] int resp_choice;  // chosen pile (1-4)
  array[T] int resp_color;   // pile matching color
  array[T] int resp_shape;   // pile matching shape
  array[T] int resp_number;  // pile matching number
}

parameters {
  real<lower=0, upper=1> r;  // reward learning rate
  real<lower=0, upper=1> p;  // punishment learning rate
  real<lower=0.01, upper=5> d;  // decision consistency
}

model {
  // Priors
  r ~ beta(1.1, 1.1);
  p ~ beta(1.1, 1.1);
  d ~ gamma(1, 0.3);
  
  // Initialize attention weights: [1/3, 1/3, 1/3]
  vector[3] a = rep_vector(1.0/3.0, 3);
  
  for (t in 1:T) {
    if (resp_choice[t] > 0) {
      
      // Create match vectors for each pile
      array[4] vector[3] m;
      for (k in 1:4) {
        m[k] = rep_vector(0.0, 3);
      }
      
      // Set match indicators
      if (resp_color[t] >= 1 && resp_color[t] <= 4) {
        m[resp_color[t]][1] = 1;  // dimension 1 = color
      }
      if (resp_shape[t] >= 1 && resp_shape[t] <= 4) {
        m[resp_shape[t]][2] = 1;  // dimension 2 = shape
      }
      if (resp_number[t] >= 1 && resp_number[t] <= 4) {
        m[resp_number[t]][3] = 1;  // dimension 3 = number
      }
      
      // Calculate choice probabilities (Equation 9 from Bishara et al.)
      // P_t,k = (sum_i (a_i^d * m_t,k,i)) / sum_j (sum_i (a_i^d * m_t,j,i))
      vector[4] pile_values;
      vector[3] a_powered;
      
      // Raise attention weights to power d
      for (dim in 1:3) {
        a_powered[dim] = pow(a[dim], d);
      }
      
      // Calculate value for each pile
      for (k in 1:4) {
        pile_values[k] = dot_product(a_powered, m[k]);
      }
      
      // Add to log-likelihood
      target += log(pile_values[resp_choice[t]]) - log_sum_exp(pile_values);
      
      // Update attention weights based on feedback
      
      // First, determine which dimension was chosen
      int chosen_dim = 0;
      if (resp_choice[t] == resp_color[t]) {
        chosen_dim = 1;  // color
      } else if (resp_choice[t] == resp_shape[t]) {
        chosen_dim = 2;  // shape
      } else if (resp_choice[t] == resp_number[t]) {
        chosen_dim = 3;  // number
      }
      
      if (chosen_dim > 0) {
        // Create feedback signal vector
        vector[3] s;
        int feedback_correct = (rew[t] == 1) ? 1 : 0;
        
        // Determine which dimensions are consistent with feedback
        int n_consistent = 0;
        array[3] int consistent_dims = {0, 0, 0};
        
        for (dim in 1:3) {
          if (feedback_correct == 1) {
            // For correct feedback, dimension is consistent if it matches
            if (m[resp_choice[t]][dim] == 1) {
              consistent_dims[dim] = 1;
              n_consistent += 1;
            }
          } else {
            // For incorrect feedback, dimension is consistent if it doesn't match
            if (m[resp_choice[t]][dim] == 0) {
              consistent_dims[dim] = 1;
              n_consistent += 1;
            }
          }
        }
        
        // Calculate feedback signal with f=1 (Equations 13-14 from paper)
        // Signal is proportional to current attention weights
        real denom = 0;
        for (dim in 1:3) {
          if (consistent_dims[dim] == 1) {
            denom += pow(a[dim], 1.0);  // f=1 means power is 1
          }
        }
        
        for (dim in 1:3) {
          if (consistent_dims[dim] == 1) {
            s[dim] = pow(a[dim], 1.0) / denom;
          } else {
            s[dim] = 0;
          }
        }
        
        // Update attention weights (Equations 6 and 8)
        if (feedback_correct == 1) {
          a = (1 - r) * a + r * s;  // Equation 6: reward learning
        } else {
          a = (1 - p) * a + p * s;  // Equation 8: punishment learning
        }
      }
    }
  }
}

generated quantities {
  real log_lik = 0;
  array[T] int y_pred;
  
  // Re-run same logic for posterior predictive
  vector[3] a = rep_vector(1.0/3.0, 3);
  
  for (t in 1:T) {
    y_pred[t] = -1;
    
    if (resp_choice[t] > 0) {
      
      array[4] vector[3] m;
      for (k in 1:4) {
        m[k] = rep_vector(0.0, 3);
      }
      
      if (resp_color[t] >= 1 && resp_color[t] <= 4) {
        m[resp_color[t]][1] = 1;
      }
      if (resp_shape[t] >= 1 && resp_shape[t] <= 4) {
        m[resp_shape[t]][2] = 1;
      }
      if (resp_number[t] >= 1 && resp_number[t] <= 4) {
        m[resp_number[t]][3] = 1;
      }
      
      vector[4] pile_values;
      vector[3] a_powered;
      
      for (dim in 1:3) {
        a_powered[dim] = pow(a[dim], d);
      }
      
      for (k in 1:4) {
        pile_values[k] = dot_product(a_powered, m[k]);
      }
      
      vector[4] pile_probs = softmax(log(pile_values));
      
      log_lik += log(pile_probs[resp_choice[t]]);
      y_pred[t] = categorical_rng(pile_probs);
      
      // Update attention
      int chosen_dim = 0;
      if (resp_choice[t] == resp_color[t]) {
        chosen_dim = 1;
      } else if (resp_choice[t] == resp_shape[t]) {
        chosen_dim = 2;
      } else if (resp_choice[t] == resp_number[t]) {
        chosen_dim = 3;
      }
      
      if (chosen_dim > 0) {
        vector[3] s;
        int feedback_correct = (rew[t] == 1) ? 1 : 0;
        
        int n_consistent = 0;
        array[3] int consistent_dims = {0, 0, 0};
        
        for (dim in 1:3) {
          if (feedback_correct == 1) {
            if (m[resp_choice[t]][dim] == 1) {
              consistent_dims[dim] = 1;
              n_consistent += 1;
            }
          } else {
            if (m[resp_choice[t]][dim] == 0) {
              consistent_dims[dim] = 1;
              n_consistent += 1;
            }
          }
        }
        
        real denom = 0;
        for (dim in 1:3) {
          if (consistent_dims[dim] == 1) {
            denom += pow(a[dim], 1.0);
          }
        }
        
        for (dim in 1:3) {
          if (consistent_dims[dim] == 1) {
            s[dim] = pow(a[dim], 1.0) / denom;
          } else {
            s[dim] = 0;
          }
        }
        
        if (feedback_correct == 1) {
          a = (1 - r) * a + r * s;
        } else {
          a = (1 - p) * a + p * s;
        }
      }
    }
  }
}
