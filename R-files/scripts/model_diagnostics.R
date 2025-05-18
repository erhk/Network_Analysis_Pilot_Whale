# Model Diagnostics for Kuramoto Model Fit

library(cmdstanr)
library(dplyr)
library(bayesplot)
library(ggplot2)
library(tidyr)
pacman::p_load(posterior)

# Load the fitted model object
fit <- readRDS("fit_kuramoto_G2_hierarchical.rds")
fit <- readRDS("fit_response_model_compx.rds")
# ─── SUMMARY DIAGNOSTICS ──────────────────────
print("\n---- Convergence Summary ----")
summary_stats <- fit$summary()
print(summary_stats %>% filter(rhat > 1.01 | ess_bulk < 100))

# ─── TRACEPLOTS ───────────────────────────────
posterior_array <- fit$draws(format = "draws_array")
mcmc_trace(posterior_array, pars = c("K", "sigma"), n_warmup = 0)

# Optionally trace a few A_raw (alpha for reponse model) and omega parameters
mcmc_trace(posterior_array, pars = c("lp__", "alpha[1]", "omega[1]", "omega[2]"), n_warmup = 0)

# ─── POSTERIOR DISTRIBUTIONS of HIERARCHICAL COUPLING MODEL: K AND SIGMA ──────────────────
posterior_df <- fit$draws(format = "draws_df")

# Plot K and sigma
mcmc_areas(posterior_array, pars = c("K", "sigma"))

# Histogram of global coupling strength K
ggplot(posterior_df, aes(x = K)) +
  geom_histogram(binwidth = 0.005, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Posterior Distribution of K (Coupling Strength)", x = "K")

# Histogram of noise parameter
ggplot(posterior_df, aes(x = sigma)) +
  geom_histogram(binwidth = 0.01, fill = "grey", color = "white") +
  theme_minimal() +
  labs(title = "Posterior Distribution of Sigma (Noise)", x = "sigma")

# ─── POSTERIOR DISTRIBUTIONS RESPONSE MODEL: OMEGA AND ALPHA ──────────────────

posterior_df %>%
  select(starts_with("omega[")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "tomato", color = "white") +
  facet_wrap(~param, scales = "free", ncol = 3) +
  labs(title = "Posterior Distributions of Omega (Responder Baselines)")

posterior_df %>%
  select(starts_with("alpha[")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~param, scales = "free", ncol = 3) +
  labs(title = "Posterior Distributions of Alpha (Initiator Influence)")


# ─── DIVERGENCE AND TREEDPTH CHECKS ───────────
print("\n---- CmdStanR Diagnostics ----")
fit$cmdstan_diagnose()

# Show treedepth hits (from metadata)
out <- fit$output_files()
cmdstanr::read_cmdstan_csv(out)$metadata
