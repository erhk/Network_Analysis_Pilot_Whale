data {
  int<lower=1> N_whales;
  int<lower=1> N_events;
  array[N_events] int<lower=1, upper=N_whales> caller;
  vector[N_events] time;
  int<lower=1> N_bouts;
  array[N_events] int<lower=1, upper=N_bouts> bout_id;

  int<lower=1> N_pairs;
  array[N_pairs] int<lower=1, upper=N_whales> pair_i; // source
  array[N_pairs] int<lower=1, upper=N_whales> pair_j; // target
}

parameters {
  vector[N_whales] omega;            // Intrinsic calling rate
  real<lower=0> K;                   // Global coupling strength
  vector[N_pairs] A_raw;             // Pairwise influences
  real<lower=0.01> sigma;            // Observation noise
}

transformed parameters {
  matrix[N_whales, N_whales] A;
  A = rep_matrix(0, N_whales, N_whales);
  for (k in 1:N_pairs) {
    A[pair_i[k], pair_j[k]] = A_raw[k];
  }
}

model {
  // Regularized priors
  omega ~ normal(0, 1);
  K ~ normal(0, 0.5);
  A_raw ~ normal(0, 0.25);
  sigma ~ normal(0, 1) T[0.01, ];

  vector[N_whales] theta_prev = rep_vector(0, N_whales);

  for (n in 2:N_events) {
    real dt = time[n] - time[n - 1];

    if (bout_id[n] != bout_id[n - 1]) {
      theta_prev = rep_vector(0, N_whales);
    }

    vector[N_whales] dtheta;
    for (i in 1:N_whales) {
      real influence = 0;
      for (j in 1:N_whales) {
        influence += A[i, j] * sin(theta_prev[j] - theta_prev[i]);
      }
      dtheta[i] = omega[i] + K * influence;
    }

    theta_prev += dt * dtheta;

    if (bout_id[n] == bout_id[n - 1]) {
      target += normal_lpdf(theta_prev[caller[n]] | theta_prev[caller[n - 1]], sigma);
    }
  }
}