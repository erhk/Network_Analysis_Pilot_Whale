// Discrete Bout-Wide Influence model
// Models how each whale's calling time is influenced by all prior calls in the same bout,
// with decaying influence over time between calls pr whale.

data {
  int<lower=1> N_events;// Total number of call events
  int<lower=1> N_whales;// Number of individual whales
  array[N_events] int<lower=1, upper=N_whales> caller; // WhaleID for each call event (indexed so stan can have fun)
  vector[N_events] time;// Time of each call
  int<lower=1> N_bouts; // Total number of bouts
  array[N_events] int<lower=1, upper=N_bouts> bout_id; // Bout ID for each event

  int<lower=1> N_pairs; // Number of directed influence pairs
  array[N_pairs] int<lower=1, upper=N_whales> pair_i; // Index for influenced whale
  array[N_pairs] int<lower=1, upper=N_whales> pair_j; // Index for influencial whale
}
                                                  
parameters {
  
  vector<lower=0.2>[N_whales] omega;   // Baseline call rate per whale (must be positive)
  vector[N_pairs] A_raw; // Raw innfluence strength from j to i
  real<lower=0.01> sigma; // Observation noise in delay durations
  //real<lower=0> lambda; // Time decay rate for influence. this might be too much, not sure if it'll work
  vector<lower=0>[N_whales] lambda;
}

transformed parameters {
  // Influence matrix A[i, j] is still from whale j to i. Start at 0, no whale starts with any influence on any other.
  matrix[N_whales, N_whales] A = rep_matrix(0, N_whales, N_whales);

  // Fill in the nonzero directed influence values from A_raw. 
  for (k in 1:N_pairs) {
    A[pair_i[k], pair_j[k]] = A_raw[k];
  }
}

model {
  // Priors --------------------------------
  // Baseline calling tendency: omega ~ N(1, 1)
  omega ~ normal(1, 1);
  
  // Pairwise influence strengths: A_raw ~ N(0, 0.5)
  // (More conservative prior to avoid model exploding in my face when fitting)
  A_raw ~ normal(0, 0.5);

  // Noise in observed call timing: sigma ~ truncated normal, for same above reason
  sigma ~ normal(0, 1) T[0.01,];

  // Decay rate for time-based influence (positive exponential prior) - no clue what will happen! Might need adjusted!
  //lambda ~ exponential(1);
  lambda ~ normal(1, 0.5) T[0, ];
  //lambda ~ lognormal(-1, 1); would like to try this in the future. No time to run this 

  

  // Likelihood --------------------------------
  for (n in 2:N_events) { // Starts at call 2, because it's alwats a comparison to previous ones, can't compare to nothing!
    int i = caller[n]; // The whale calling at time[n]
    real influence_sum = 0; // Accumulated influence on whale [i]

    // You look at all earlier call events (m < n). Loop through all previous events (before current event n). See if any earlier calls (influence of these) affect currnt ones call
    for (m in 1:(n - 1)) { 
      if (bout_id[m] == bout_id[n]) { // Seperate into bouts. Only include calls from the same bout
        int j = caller[m]; // Whale[j] who called at time[m]
        real delta_t = time[n] - time[m];  // Time since previous call

        if (delta_t > 0) {
        // Use whale-specific decay rate lambda[i]
          real decay = exp(-lambda[i] * delta_t);

        // Accumulate influence from j to i, decayed by i's memory
          influence_sum += A[i, j] * decay;
        }
        // Previous model, models global decay, not whale level 
        //if (delta_t > 0) { // Guard against zero or negative lag (same timestamp)
          
          // Close to riccardo coursenotes, but with help from chatgpt
          // Exponential decay: recent calls have stronger effec. 
          // Apply exponential decay:
          //More recent calls → higher decay value → stronger influence.
          //Older calls → lower decay → diminished effect. Lambda controls how fast this decay happens.
          //real decay = exp(-lambda * delta_t);

          // Accumulate influence from j to i, decayed by time lag. Multiply by decay to downweight older influence.
          // Add to the running total influence_sum, that represents all the cumulative influence on whale i before event n in this bout.
          //influence_sum += A[i, j] * decay;
        //}
      }
    }

    // Predict expected delay between calls based on current whale's base rate and influence
    // (Shorter expected delay if influence is strong)
    real predicted_delay = 1.0 / (omega[i] + influence_sum);

    // Observed delay between current and previous call
    real observed_delay = time[n] - time[n - 1];

    // Add to likelihood: observed delay ~ Normal(predicted_delay, sigma)
    target += normal_lpdf(observed_delay | predicted_delay, sigma);
  }
}
