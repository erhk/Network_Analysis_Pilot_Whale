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
  vector[N_whales]  omega;          // intrinsic call rates
  // real<lower=0>     K;              // global coupling scale

  // --- hierarchical shrinkage for the pairwise influences ---
  vector[N_pairs]   z_A;            // standard-normal
  // real<lower=0>     tau;            // global SD  (to be learned)
  
  
  real<lower=0.1, upper=2> tau; // Try to keep from 0
  real<lower=0.1, upper=2> K;   // Same as tau

  real<lower=0.05>  sigma;          // observation noise
}

transformed parameters {
  // non-centred parameterisation:  A_raw = tau * z_A
  vector[N_pairs]   A_raw = tau * z_A;

  // rebuild A matrix
  matrix[N_whales,N_whales] A = rep_matrix(0, N_whales, N_whales);
  for (k in 1:N_pairs)
    A[pair_i[k], pair_j[k]] = A_raw[k];
}

model {
  // priors
  omega ~ normal(0, 1);
  K     ~ normal(0.5, 0.3); // soften prior here

  z_A   ~ normal(0, 1);              // standard-normal
  //tau   ~ exponential(0.1);       // strong pooling; adjust if too tight. Leaving out for now, tryong normal
  // tau ~ normal(0.3, 0.2) T[0.05,];  
  tau ~ normal(0.5, 0.3); // soften prior
  sigma ~ normal(0, 1) T[0.01,];

  // -------- likelihood (unchanged) --------
  vector[N_whales] theta_prev = rep_vector(0, N_whales);

  for (n in 2:N_events) {
    real dt = time[n] - time[n-1];

    if (bout_id[n] != bout_id[n-1])
      theta_prev = rep_vector(0, N_whales);

    vector[N_whales] dtheta;
    for (i in 1:N_whales) {
      real influence = 0;
      for (j in 1:N_whales)
        influence += A[i,j] * sin(theta_prev[j] - theta_prev[i]);
        // DEBUG PRINT - will crash R when run
  // print("event ", n, ", whale ", i, ", influence = ", influence, ", dtheta = ", omega[i] + K * influence);

      dtheta[i] = omega[i] + K * influence;
    }
    theta_prev += dt * dtheta;

    if (bout_id[n] == bout_id[n-1])
      target += normal_lpdf(theta_prev[caller[n]] |
                            theta_prev[caller[n-1]], sigma);
  }
}
