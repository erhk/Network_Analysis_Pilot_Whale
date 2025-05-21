// Discrete Phase Influence Model - Stepwise version

data {
  int<lower=1> N_events;
  int<lower=1> N_whales;
  array[N_events] int<lower=1, upper=N_whales> caller;
  vector[N_events] time;
  int<lower=1> N_bouts;
  array[N_events] int<lower=1, upper=N_bouts> bout_id;
}

parameters {
  vector[N_whales] omega;               // intrinsic calling rate
  real<lower=0.01> sigma;               // observation noise
}

model {
  // Priors
  omega ~ normal(0, 1);
  sigma ~ normal(0, 1) T[0.01,];

  // Likelihood (only temporal smoothness for now)
  for (n in 2:N_events) {
    if (bout_id[n] == bout_id[n-1]) {
      target += normal_lpdf(time[n] - time[n-1] |
                            1.0 / omega[caller[n]], sigma);
    }
  }
}
