data {
  int<lower=1> T;
  array[T] int<lower=-1, upper=1> rew;   // -1,0,1
  array[T] int<lower=-1, upper=1> los;   // -1,0,1 (spesso 0 o -1)
  array[T] int resp_choice;
  array[T] int resp_color;
  array[T] int resp_shape;
  array[T] int resp_number;
}

parameters {
  real<lower=1e-4, upper=1-1e-4> alpha_pos; // LR esiti positivi
  real<lower=1e-4, upper=1-1e-4> alpha_neg; // LR esiti negativi
  real<lower=1e-6> beta;                     // inverse temperature
  real kappa;                                // sticky-choice (può essere +/-)
  real<lower=0, upper=0.5> lapse;            // probabilità di risposta casuale
}

model {
  // Priors regolarizzanti
  alpha_pos ~ beta(2, 2);
  alpha_neg ~ beta(2, 2);
  beta      ~ lognormal(log(2.5), 0.5);      // evita β≈0
  kappa     ~ normal(0, 1);                  // sticky moderato
  lapse     ~ beta(2, 20);                   // medio ~0.09

  vector[3] Q = rep_vector(0.0, 3);
  int last_dim = 0;

  for (t in 1:T) {
    if (resp_choice[t] > 0) {
      vector[3] logits = beta * Q;
      if (last_dim > 0) logits[last_dim] += kappa;

      vector[3] soft = softmax(logits);
      vector[3] dim_probs = (1 - lapse) * soft + lapse / 3;

      int chosen_dim = 0;
      if      (resp_choice[t] == resp_color[t])  chosen_dim = 1;
      else if (resp_choice[t] == resp_shape[t])  chosen_dim = 2;
      else if (resp_choice[t] == resp_number[t]) chosen_dim = 3;

      if (chosen_dim > 0) {
        target += log(dim_probs[chosen_dim]);

        real outcome = rew[t] + los[t]; // in {-1,0,1}
        if (outcome > 0) {
          Q[chosen_dim] += alpha_pos * (outcome - Q[chosen_dim]);
        } else {
          Q[chosen_dim] += alpha_neg * (outcome - Q[chosen_dim]);
        }
        last_dim = chosen_dim;
      }
    }
  }
}

generated quantities {
  real log_lik = 0;
  array[T] int y_pred;
  vector[3] Q = rep_vector(0.0, 3);
  int last_dim = 0;

  for (t in 1:T) {
    y_pred[t] = -1;
    if (resp_choice[t] > 0) {
      vector[3] logits = beta * Q;
      if (last_dim > 0) logits[last_dim] += kappa;
      vector[3] soft = softmax(logits);
      vector[3] dim_probs = (1 - lapse) * soft + lapse / 3;

      int chosen_dim = 0;
      if      (resp_choice[t] == resp_color[t])  chosen_dim = 1;
      else if (resp_choice[t] == resp_shape[t])  chosen_dim = 2;
      else if (resp_choice[t] == resp_number[t]) chosen_dim = 3;

      if (chosen_dim > 0) {
        log_lik += log(dim_probs[chosen_dim]);
        int sim_dim = categorical_rng(dim_probs);
        if      (sim_dim == 1) y_pred[t] = resp_color[t];
        else if (sim_dim == 2) y_pred[t] = resp_shape[t];
        else                   y_pred[t] = resp_number[t];

        real outcome = rew[t] + los[t];
        if (outcome > 0) {
          Q[chosen_dim] += alpha_pos * (outcome - Q[chosen_dim]);
        } else {
          Q[chosen_dim] += alpha_neg * (outcome - Q[chosen_dim]);
        }
        last_dim = chosen_dim;
      }
    }
  }
}
