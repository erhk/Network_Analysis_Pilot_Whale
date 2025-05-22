data {
  int<lower=1> N_whales;         // Number of whales
  int<lower=1> N_events;         // Number of call events
  array[N_events] int<lower=1, upper=N_whales> caller;  // Whale making the call
  vector[N_events] time;                                // Call time (relative to bout start)
  int<lower=1> N_bouts;
  array[N_events] int<lower=1, upper=N_bouts> bout_id;  // Bout index for each call
}

parameters {
  vector[N_whales] omega;           // Intrinsic call frequencies
  real<lower=0> K;                  // Global coupling strength
  real<lower=0> A_shared;           // Shared pairwise coupling strength
  real<lower=0.01> sigma;           // Observation noise (bounded to avoid instability)
}

model {
  // Priors
  omega ~ normal(0, 1);
  K ~ normal(0, 0.5);
  A_shared ~ normal(0, 0.5);
  sigma ~ normal(0, 1) T[0.01, ];

  // Phase state for each whale
  vector[N_whales] theta_prev = rep_vector(0, N_whales);

  for (n in 2:N_events) {
    real dt = time[n] - time[n - 1];

    // Reset phase if new bout starts
    //if (bout_id[n] != bout_id[n - 1]) {
      //theta_prev = rep_vector(0, N_whales);
    //}
    if (bout_id[n] != bout_id[n - 1]) {
      theta_prev = rep_vector(0, N_whales);
      continue;  // Skip to next step — don't compute anything across bouts
    }

    // Phase dynamics for each whale
    vector[N_whales] dtheta;
    for (i in 1:N_whales) {
      real influence = 0;
      for (j in 1:N_whales) {
        influence += sin(theta_prev[j] - theta_prev[i]);
      }
      dtheta[i] = omega[i] + K * A_shared * influence;
    }

    // Update phase
    theta_prev += dt * dtheta;

    // Likelihood for within-bout transitions
    if (bout_id[n] == bout_id[n - 1]) {
      target += normal_lpdf(theta_prev[caller[n]] | theta_prev[caller[n - 1]], sigma);
    }
  }
}
