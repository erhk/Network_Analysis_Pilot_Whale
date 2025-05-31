# Discrete bout-wide Influence Model

pacman::p_load(tidyverse, cmdstanr, tidybayes, dplyr, ggplot2, igraph, ggraph, bayesplot
)

# Load data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")


#### --------------- STAN data prep ---------------- ####
# Filter for Group G1
g1_data <- prep %>%
  filter(Group == "G1") %>%
  arrange(BoutStartTime, time = StartTime)  # Ensure temporal order

g1_data <- g1_data %>%
  filter(Group == "G1") %>%
  arrange(bout, StartTime) %>%
  mutate(
    WhaleIndex = as.integer(factor(WhaleID)), # tan can't handle characters
    BoutIndex = bout  # Not really needed to change, but it's just to keep the "index" seperate
  )


# Define number of whales and events -  G3
whales <- sort(unique(g1_data$WhaleIndex))
N_whales <- length(whales)
N_events <- nrow(g1_data)

# Generate influence pairs 
pairs <- expand.grid(i = whales, j = whales) %>%
  filter(i != j)  # skip self-influence (might matter, but i'm looking at social influence. So only A[i,j], skipping A[i,i]

# Stan data - G3
stan_data_g1 <- list(
  N_events = N_events,
  N_whales = N_whales,
  caller = g1_data$WhaleIndex,
  time = g1_data$StartTime,
  N_bouts = max(g1_data$BoutIndex),
  bout_id = g1_data$BoutIndex,
  N_pairs = nrow(pairs),
  pair_i = pairs$i,
  pair_j = pairs$j
)

# Save stan data list
#saveRDS(stan_data_g3, "g3_stan_data.rds")
#saveRDS(stan_data_g2, "g2_stan_data.rds")
#saveRDS(stan_data_g1, "g1_stan_data.rds")
#

# Save whale Id mapping
whale_id_map <- g1_data %>%
  select(WhaleID, WhaleIndex) %>%
  distinct() %>%
  arrange(WhaleIndex)

#saveRDS(whale_id_map, "g1_whale_id_map.rds")

# Save OG subsetted data
#saveRDS(g1_data, "g1_data.rds")




#### --------------- Fit Models  ------------------- ####

mod <- cmdstan_model("../models/discrete bout-wide Influence model.stan", cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = stan_data_g2,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.999, # went as close to 1 as possible
  max_treedepth = 20,
  refresh = 5,
  seed = 1990
)

#fit$save_object("fit_G1_lambda_normprior_all.rds")



#### --- Check model fits, especially lambda, because it's being silly :(

# Plot Posterior for lambda 
fit <- readRDS("fit_G1_lambda_normprior_all.rds")


# Extract posterior draws for lambda
lambda_draws <- fit_g1$draws("lambda", format = "draws_df") %>%
  spread_draws(lambda[i])

# Plot posterior density
ggplot(lambda_draws, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 1) +
  labs(
    title = "Posterior Distribution of λ (Decay Rate)",
    x = expression(lambda),
    y = "Posterior Density"
  ) +
  theme_minimal(base_size = 14)



# ---- Model Checks, filter for bad values 
# Diagnostics for lambda because it looks very strange, bimodal, with two peaks near 0.2 and 1.7, valley close to 1

fit$summary("lambda")
fit$cmdstan_diagnose()
posterior::rhat(fit$draws("lambda"))

# Check chain mixing
mcmc_trace(fit$draws("lambda"))

# check rhat - ideally all above 1
fit$summary() |>
  dplyr::filter(rhat > 1.01)

# check ess - ideally below 400
fit$summary() |>
  dplyr::filter(ess_bulk < 400 | ess_tail < 400)



