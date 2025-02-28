data {
  int<lower=1> N; // Total number of trials
  int<lower=1> S; // Number of subjects
  int<lower=1> K; // Number of sessions
  array[N] int<lower=0, upper=1> choice; // Choices (0 or 1)
  array[N] int<lower=0, upper=1> feedback; // Feedback (0 or 1)
  array[N] int<lower=1, upper=S> subject; // Subject IDs
  array[N] int<lower=1, upper=K> session; // Session IDs
  array[S] int<lower=1, upper=2> condition; // Condition IDs (1 or 2)
  array[N] int<lower=1, upper=5> mood; // Discrete trial-by-trial mood ratings (1 to 5)
}
parameters {
  real<lower=0, upper=1> grand_mu_alpha;
  real<lower=0> grand_mu_beta;
  vector<lower=0, upper=1>[2] condition_mu_alpha; // Condition-specific effects on learning rate
  vector<lower=0>[2] condition_mu_beta; // Condition-specific effects on decision-making
  vector<lower=0>[2] sigma_condition_alpha; // SD of condition-specific alpha deviations
  vector<lower=0>[2] sigma_condition_beta; // SD of condition-specific beta deviations
  
  vector[S] z_alpha_raw; // Raw deviations in alpha for subjects
  vector[S] z_beta_raw; // Raw deviations in beta for subjects
  
  vector[K] session_alpha_raw; // Raw session-level deviations in alpha
  vector[K] session_beta_raw; // Raw session-level deviations in beta
  
  // Thresholds for the ordinal mood effect
  ordered[4] mood_thresholds_alpha; // Ordered thresholds for alpha (4 thresholds separate 5 mood levels)
  ordered[4] mood_thresholds_beta; // Ordered thresholds for beta (4 thresholds separate 5 mood levels)
  
  // Standard deviations for mood effects
  real<lower=0> sigma_mood_alpha;
  real<lower=0> sigma_mood_beta;
  
  // Subject-specific mood effects
  vector[S] mood_effect_alpha; // Mood effects on alpha for each subject
  vector[S] mood_effect_beta; // Mood effects on beta for each subject
  
  // Interaction terms for mood × condition
  real mood_condition_alpha_interaction; // Interaction term for mood × condition on alpha
  real mood_condition_beta_interaction; // Interaction term for mood × condition on beta
  
  // Resistance terms for both conditions
  array[2] real resistance_autobiographical; // Resistance for both groups
}
transformed parameters {
  vector<lower=0, upper=1>[S] alpha;
  vector<lower=0>[S] beta;
  
  for (s in 1 : S) {
    // Non-centered parameterization for subject-level alpha and beta
    real subject_alpha = condition_mu_alpha[condition[s]]
                         + sigma_condition_alpha[condition[s]]
                           * z_alpha_raw[s];
    real subject_beta = condition_mu_beta[condition[s]]
                        + sigma_condition_beta[condition[s]] * z_beta_raw[s];
    
    // Non-centered parameterization for session-level alpha and beta
    real session_alpha_s = session_alpha_raw[session[s]]; // Assumes session_alpha_raw is scaled properly
    real session_beta_s = session_beta_raw[session[s]]; // Assumes session_beta_raw is scaled properly
    
    // Condition-specific resistance
    real resistance_term = (condition[s] == 1)
                           ? resistance_autobiographical[1]
                           : resistance_autobiographical[2];
    
    // Ordinal mood effect based on the thresholds
    real mood_effect_alpha_s = mood_effect_alpha[s]
                               + mood_condition_alpha_interaction
                                 * (condition[s] == 1);
    real mood_effect_beta_s = mood_effect_beta[s]
                              + mood_condition_beta_interaction
                                * (condition[s] == 1);
    
    // Learning rate (alpha) calculation with bounding to prevent NaNs
    real temp_alpha = logit(subject_alpha) + session_alpha_s
                      + mood_effect_alpha_s - resistance_term;
    
    // Bound to avoid extreme values that cause NaNs
    temp_alpha = fmin(fmax(temp_alpha, -15), 15);
    
    // Apply inverse logit transformation to get alpha in [0, 1]
    alpha[s] = inv_logit(temp_alpha);
    
    // Print debug information for alpha
    if (is_nan(alpha[s])) {
      print("alpha for subject ", s, " is NaN with inputs: ",
            "condition_mu_alpha=", condition_mu_alpha[condition[s]],
            ", sigma_condition_alpha=", sigma_condition_alpha[condition[s]],
            ", z_alpha_raw=", z_alpha_raw[s], ", session_alpha_raw=",
            session_alpha_raw[session[s]], ", mood_effect_alpha_s=",
            mood_effect_alpha_s, ", resistance_term=", resistance_term);
    }
    
    // Inverse temperature (beta) calculation
    real temp_beta = subject_beta + session_beta_s + mood_effect_beta_s;
    temp_beta = fmax(temp_beta, -15); // Lower bound to avoid NaNs
    beta[s] = exp(temp_beta);
    
    // Print debug information for beta
    if (is_nan(beta[s])) {
      print("beta for subject ", s, " is NaN with inputs: ",
            "condition_mu_beta=", condition_mu_beta[condition[s]],
            ", sigma_condition_beta=", sigma_condition_beta[condition[s]],
            ", z_beta_raw=", z_beta_raw[s], ", session_beta_raw=",
            session_beta_raw[session[s]], ", mood_effect_beta_s=",
            mood_effect_beta_s);
    }
  }
}
model {
  // Priors
  grand_mu_alpha ~ beta(2, 2);
  grand_mu_beta ~ normal(1, 1);
  
  condition_mu_alpha ~ beta(2, 2);
  condition_mu_beta ~ normal(1, 1);
  
  sigma_condition_alpha ~ normal(0, 1);
  sigma_condition_beta ~ normal(0, 1);
  
  // Non-centered parameterization: raw values from normal(0, 1)
  z_alpha_raw ~ std_normal();
  z_beta_raw ~ std_normal();
  
  session_alpha_raw ~ std_normal(); // Raw values for session-level deviations
  session_beta_raw ~ std_normal(); // Raw values for session-level deviations
  
  // Ordinal mood effect thresholds (must be ordered)
  mood_thresholds_alpha ~ normal(0, 1);
  mood_thresholds_beta ~ normal(0, 1);
  
  // Priors for mood effects
  sigma_mood_alpha ~ normal(0, 1);
  sigma_mood_beta ~ normal(0, 1);
  
  mood_effect_alpha ~ normal(0, sigma_mood_alpha);
  mood_effect_beta ~ normal(0, sigma_mood_beta);
  
  // Mood × condition interaction terms
  mood_condition_alpha_interaction ~ normal(0, 1);
  mood_condition_beta_interaction ~ normal(0, 1);
  
  // Resistance term prior
  resistance_autobiographical ~ normal(0, 1);
  
  // Likelihood
  for (n in 1 : N) {
    vector[2] v = rep_vector(0.5, 2); // Initialize values
    
    real feedback_mood_effect = feedback[n] * mood[n]; // Feedback-mood interaction
    
    if (condition[subject[n]] == 1) {
      // Autobiographical group
      v[choice[n] + 1] += (alpha[subject[n]] + feedback_mood_effect)
                          * (feedback[n] - v[choice[n] + 1]);
    } else {
      // Fractal group
      v[choice[n] + 1] += alpha[subject[n]]
                          * (feedback[n] - v[choice[n] + 1]);
    }
    
    real p = inv_logit(beta[subject[n]] * (v[2] - v[1]));
    choice[n] ~ bernoulli(p); // Likelihood of choice
  }
}
generated quantities {
  real mood_alpha_diff; // Difference in mood effect on learning (alpha) between the two groups
  real mood_beta_diff; // Difference in mood effect on decision-making (beta) between the two groups
  
  mood_alpha_diff = mood_condition_alpha_interaction; // Mood × condition effect on alpha
  mood_beta_diff = mood_condition_beta_interaction; // Mood × condition effect on beta
}
