data {
  int<lower=1> N;  // number of subjects
  int<lower=1> T;  // max trials
  array[N] int<lower=1, upper=T> Tsubj;
  array[N] int<lower=0, upper=1> holdout;
  
  array[N, T] int<lower=-1, upper=1> rew;  // +1 for correct
  array[N, T] int<lower=-1, upper=1> los;  // -1 for incorrect
  array[N, T] int rule_choice;
  array[N, T] int resp_choice;
  array[N, T] int resp_color;
  array[N, T] int resp_shape;
  array[N, T] int resp_number;
}

parameters {
  // Subject-level parameters (no hierarchy for maximum differentiation)
  vector<lower=0, upper=1>[N] r;  // reward learning rate
  vector<lower=0, upper=1>[N] p;  // punishment learning rate
  vector<lower=0.01, upper=5>[N] d;  // decision consistency
}

transformed parameters {
  // d_power[i] will store d[i]^d[i] for efficiency
  vector[N] log_d = log(d);
}

model {
  // Weakly informative priors
  r ~ beta(1.1, 1.1);
  p ~ beta(1.1, 1.1);
  d ~ gamma(1.0, 0.3);
  
  // Likelihood
  for (i in 1:N) {
    if (holdout[i] == 0) {
      // Attention weights start equal: [1/3, 1/3, 1/3]
      vector[3] a = rep_vector(1.0/3.0, 3);
      
      for (t in 1:Tsubj[i]) {
        if (resp_choice[i, t] > 0) {
          
          // Match vectors for each pile
          array[4] vector[3] m;
          for (k in 1:4) {
            m[k] = rep_vector(0.0, 3);
          }
          
          // Pile that matches color
          if (resp_color[i, t] >= 1 && resp_color[i, t] <= 4) {
            m[resp_color[i, t]][1] = 1;
          }
          // Pile that matches shape
          if (resp_shape[i, t] >= 1 && resp_shape[i, t] <= 4) {
            m[resp_shape[i, t]][2] = 1;
          }
          // Pile that matches number
          if (resp_number[i, t] >= 1 && resp_number[i, t] <= 4) {
            m[resp_number[i, t]][3] = 1;
          }
          
          // Calculate choice probabilities using Eq. 9 from paper
          // P_t,k = (sum_i (a_i^d * m_t,k,i)) / sum_j (sum_i (a_i^d * m_t,j,i))
          vector[4] pile_values;
          vector[3] a_powered;
          
          for (dim in 1:3) {
            a_powered[dim] = pow(a[dim], d[i]);
          }
          
          for (k in 1:4) {
            pile_values[k] = dot_product(a_powered, m[k]);
          }
          
          // Add to log-likelihood
          target += log(pile_values[resp_choice[i, t]]) - log_sum_exp(pile_values);
          
          // Update attention weights
          // First, determine which dimension was chosen
          int chosen_dim = 0;
          if (resp_choice[i, t] == resp_color[i, t]) {
            chosen_dim = 1;
          } else if (resp_choice[i, t] == resp_shape[i, t]) {
            chosen_dim = 2;
          } else if (resp_choice[i, t] == resp_number[i, t]) {
            chosen_dim = 3;
          }
          
          if (chosen_dim > 0) {
            // Feedback signal vector
            vector[3] s;
            int feedback_correct = (rew[i, t] == 1) ? 1 : 0;
            
            // Count how many dimensions are consistent with feedback
            int n_consistent = 0;
            array[3] int consistent_dims = {0, 0, 0};
            for (dim in 1:3) {
              if (feedback_correct == 1) {
                if (m[resp_choice[i, t]][dim] == 1) {
                  consistent_dims[dim] = 1;
                  n_consistent += 1;
                }
              } else {
                if (m[resp_choice[i, t]][dim] == 0) {
                  consistent_dims[dim] = 1;
                  n_consistent += 1;
                }
              }
            }
            
            // Calculate feedback signal with f=1 (proportional to attention)
            // Using Equations 13-14 from paper
            real denom = 0;
            for (dim in 1:3) {
              if (consistent_dims[dim] == 1) {
                denom += pow(a[dim], 1.0);  // f=1
              }
            }
            
            for (dim in 1:3) {
              if (consistent_dims[dim] == 1) {
                s[dim] = pow(a[dim], 1.0) / denom;
              } else {
                s[dim] = 0;
              }
            }
            
            // Update attention with appropriate learning rate
            if (feedback_correct == 1) {
              a = (1 - r[i]) * a + r[i] * s;  // Equation 6
            } else {
              a = (1 - p[i]) * a + p[i] * s;  // Equation 8
            }
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
    vector[3] a = rep_vector(1.0/3.0, 3);
    
    for (t in 1:Tsubj[i]) {
      y_pred[i, t] = -1;
      
      if (resp_choice[i, t] > 0) {
        
        // Match vectors
        array[4] vector[3] m;
        for (k in 1:4) {
          m[k] = rep_vector(0.0, 3);
        }
        if (resp_color[i, t] >= 1 && resp_color[i, t] <= 4) {
          m[resp_color[i, t]][1] = 1;
        }
        if (resp_shape[i, t] >= 1 && resp_shape[i, t] <= 4) {
          m[resp_shape[i, t]][2] = 1;
        }
        if (resp_number[i, t] >= 1 && resp_number[i, t] <= 4) {
          m[resp_number[i, t]][3] = 1;
        }
        
        // Choice probabilities
        vector[4] pile_values;
        vector[3] a_powered;
        for (dim in 1:3) {
          a_powered[dim] = pow(a[dim], d[i]);
        }
        for (k in 1:4) {
          pile_values[k] = dot_product(a_powered, m[k]);
        }
        
        vector[4] pile_probs = softmax(log(pile_values));
        
        // Log-likelihood
        log_lik[i] += log(pile_probs[resp_choice[i, t]]);
        
        // Posterior prediction
        y_pred[i, t] = categorical_rng(pile_probs);
        
        // Update attention
        int chosen_dim = 0;
        if (resp_choice[i, t] == resp_color[i, t]) {
          chosen_dim = 1;
        } else if (resp_choice[i, t] == resp_shape[i, t]) {
          chosen_dim = 2;
        } else if (resp_choice[i, t] == resp_number[i, t]) {
          chosen_dim = 3;
        }
        
        if (chosen_dim > 0) {
          vector[3] s;
          int feedback_correct = (rew[i, t] == 1) ? 1 : 0;
          
          int n_consistent = 0;
          array[3] int consistent_dims = {0, 0, 0};
          for (dim in 1:3) {
            if (feedback_correct == 1) {
              if (m[resp_choice[i, t]][dim] == 1) {
                consistent_dims[dim] = 1;
                n_consistent += 1;
              }
            } else {
              if (m[resp_choice[i, t]][dim] == 0) {
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
            a = (1 - r[i]) * a + r[i] * s;
          } else {
            a = (1 - p[i]) * a + p[i] * s;
          }
        }
      }
    }
  }
}

