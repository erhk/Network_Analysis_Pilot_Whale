// Stan model: Response-call latency influenced by initiator whale
// Goal: Estimate directional influence A_ij based on latency of response from i to j

// The response-focused Stan model is ready. It estimates how much each initiator whale 
// influences the latency of responses from other whales.

// This model fits:

// - omega[i]: baseline latency of whale i as a responder

// - alpha[j]: how much whale j reduces others' response time when initiating

// -  sigma: observation noise


data {
  int<lower=1> N;                       // number of response events
  vector[N] latency;                   // latency between initiation and response
  int<lower=1> N_whales;               // number of whales
  array[N] int<lower=1, upper=N_whales> responder;  // whale who responded
  array[N] int<lower=1, upper=N_whales> initiator;  // whale who initiated the bout
}

parameters {
  vector[N_whales] omega;              // baseline latency per responder
  vector[N_whales] alpha;              // influence level per initiator
  real<lower=0> sigma;                 // observation noise
}

model {
  // Priors
  omega ~ normal(0, 1);                // each responder has their own baseline latency
  alpha ~ normal(0, 1);                // each initiator has their influence strength
  sigma ~ exponential(1);

  // Likelihood: latency reduced by initiator's influence
  for (n in 1:N)
    latency[n] ~ normal(omega[responder[n]] - alpha[initiator[n]], sigma);
}
