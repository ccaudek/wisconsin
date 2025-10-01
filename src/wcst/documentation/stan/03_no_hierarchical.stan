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
  // NESSUNA GERARCHIA - stima indipendente per ogni soggetto
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0>[N] beta;
  vector<lower=0>[N] inertia;
}

model {
  // Prior indipendenti (non informativi)
  alpha ~ beta(1, 1);  // uniform
  beta ~ gamma(1, 0.5);  // mean=2, molto vago
  inertia ~ gamma(1, 0.5);
  
  // Likelihood identico
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
  
  for (i in 1:N) {
    log_lik[i] = 0;
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
          log_lik[i] += log(dim_probs[chosen_dim]);
          real outcome = rew[i, t] + los[i, t];
          Q[chosen_dim] += alpha[i] * (outcome - Q[chosen_dim]);
          last_dim = chosen_dim;
        }
      }
    }
  }
}
