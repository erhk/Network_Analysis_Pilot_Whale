data {
  int<lower=1> N_whales;         // Number of whales
  int<lower=1> N_events;         // Number of call events
  int<lower=1> N_pairs;          // Number of pairwise connections

  array[N_events] int<lower=1, upper=N_whales> caller;  // Whale making the call
  vector[N_events] time;                                // Call time (since bout start)
  array[N_pairs] int<lower=1, upper=N_whales> pair_i;   // From whale (in A[i,j])
  array[N_pairs] int<lower=1, upper=N_whales> pair_j;   // To whale (in A[i,j])

  int<lower=1> N_bouts;
  array[N_events] int<lower=1, upper=N_bouts> bout_id;  // Bout index
}

parameters {
  vector[N_whales] omega;            // Intrinsic call frequency
  real<lower=0> K;                   // Global coupling strength
  vector[N_pairs] A_raw;             // Pairwise connection strengths
  real<lower=0.01> sigma;            // Observation noise, bounded to prevent NaNs
}

transformed parameters {
  matrix[N_whales, N_whales] A;
  A = rep_matrix(0, N_whales, N_whales);
  for (k in 1:N_pairs) {
    A[pair_i[k], pair_j[k]] = A_raw[k];
  }
}

model {
  // Regularized priors for better geometry
  omega ~ normal(0, 1);
  A_raw ~ normal(0, 0.5);  // tighter prior to avoid overfitting
  K ~ normal(0, 0.5);
  sigma ~ normal(0, 1) T[0.01, ];  // reinforce lower bound

  // Phase tracker
  vector[N_whales] theta_prev = rep_vector(0, N_whales);

  for (n in 2:N_events) {
    real dt = time[n] - time[n - 1];

    // Reset phase if bout changes
    if (bout_id[n] != bout_id[n - 1]) {
      theta_prev = rep_vector(0, N_whales);
    }

    // Phase update
    vector[N_whales] dtheta;
    for (i in 1:N_whales) {
      real influence = 0;
      for (j in 1:N_whales) {
        influence += A[i, j] * sin(theta_prev[j] - theta_prev[i]);
      }
      dtheta[i] = omega[i] + K * influence;
    }

    theta_prev += dt * dtheta;

    // Likelihood: only if caller[n] and caller[n - 1] are in same bout
    if (bout_id[n] == bout_id[n - 1]) {
      target += normal_lpdf(theta_prev[caller[n]] | theta_prev[caller[n - 1]], sigma);
    }
  }
}
