data {
  int<lower=1> T;  // trials for THIS subject
  array[T] real rew;
  array[T] real los;
  array[T] int resp_choice;
  array[T] int resp_color;
  array[T] int resp_shape;
  array[T] int resp_number;
}

parameters {
  real<lower=0, upper=1> alpha;  // learning rate
  real<lower=0> beta;            // inverse temperature
}

model {
  // Weak priors (quasi-uniform)
  alpha ~ beta(1.1, 1.1);
  beta ~ gamma(1.5, 0.5);
  
  // Likelihood
  vector[3] Q = rep_vector(0.0, 3);
  
  for (t in 1:T) {
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
        target += log(dim_probs[chosen_dim]);
        real outcome = rew[t] + los[t];
        Q[chosen_dim] += alpha * (outcome - Q[chosen_dim]);
      }
    }
  }
}
