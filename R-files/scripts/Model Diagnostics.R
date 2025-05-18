# Model Diagnostics for Kuramoto Model Fit
# ---------------------- Libraries
pacman::p_load(cmdstanr, dplyr, posterior, bayesplot, ggplot2)

# Load the fitted model object
fit <- readRDS("fit_kuramoto_G1_small.rds")

# ─── SUMMARY DIAGNOSTICS ──────────────────────
print("\n---- Convergence Summary ----")
summary_stats <- fit$summary()
print(summary_stats %>% filter(rhat > 1.01 | ess_bulk < 100))

# ─── TRACEPLOTS ───────────────────────────────
posterior_array <- fit$draws(format = "draws_array")
mcmc_trace(posterior_array, pars = c("K", "sigma"), n_warmup = 0)

# Optionally trace a few A_raw and omega parameters
mcmc_trace(posterior_array, pars = c("A_raw[1]", "A_raw[2]", "omega[1]", "omega[2]"), n_warmup = 0)

# ─── POSTERIOR DISTRIBUTIONS ──────────────────
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

# ─── DIVERGENCE AND TREEDPTH CHECKS ───────────
print("\n---- CmdStanR Diagnostics ----")
fit$cmdstan_diagnose()

# Show treedepth hits (from metadata)
out <- fit$output_files()
cmdstanr::read_cmdstan_csv(out)$metadata
