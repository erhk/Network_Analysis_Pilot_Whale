// Discrete Phase Influence Model - Adding pairwise influence

data {
  int<lower=1> N_events;
  int<lower=1> N_whales;
  array[N_events] int<lower=1, upper=N_whales> caller;
  vector[N_events] time;
  int<lower=1> N_bouts;
  array[N_events] int<lower=1, upper=N_bouts> bout_id;

  int<lower=1> N_pairs;
  array[N_pairs] int<lower=1, upper=N_whales> pair_i;
  array[N_pairs] int<lower=1, upper=N_whales> pair_j;
}

//parameters {
//  vector<lower=0.1>[N_whales] omega;      // intrinsic calling rates
//  vector[N_pairs] A_raw;                  // pairwise influence (unconstrained)
//  real<lower=0.01> sigma;                 // observation noise
//}

//transformed parameters {
//  matrix[N_whales, N_whales] A = rep_matrix(0, N_whales, N_whales);
//  for (k in 1:N_pairs)
//    A[pair_i[k], pair_j[k]] = A_raw[k];
//}

//model {
  // Priors
  //omega ~ normal(1, 1);
  //A_raw ~ normal(0, 1);
  //sigma ~ normal(0, 1) T[0.01,];

// ---- Try to model omega a bit differently. Above chunks do work though
// Changes: Whales’ baseline call tendencies (omega) are drawn from a group-level 
// distribution with mean mu_omega and spread sigma_omega.
parameters {
  real mu_omega;
  real<lower=0> sigma_omega;
  vector<lower=0>[N_whales] omega_raw;
}
transformed parameters {
  vector[N_whales] omega = mu_omega + sigma_omega * omega_raw;
}
model {
  mu_omega ~ normal(0.01, 0.01);
  sigma_omega ~ exponential(1);
  omega_raw ~ normal(0, 1);
}
// ------
  // Likelihood: event timing depends on omega and summed influence from others
  for (n in 2:N_events) {
    if (bout_id[n] == bout_id[n - 1]) {
      int i = caller[n];
      int j_prev = caller[n - 1];

      // Influence from previous caller to current caller
      real influence = A[i, j_prev];
      real predicted_delay = 1.0 / (omega[i] + influence);

      target += normal_lpdf(time[n] - time[n - 1] | predicted_delay, sigma);
    }
  }
}
