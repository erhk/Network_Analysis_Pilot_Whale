# Discrete Phase Influence Model - Adding pairwise influence
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr)
library(ggplot2)

# Load data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")

# Preprocess function
prepare_discrete_phase_with_influence <- function(group_name, prep_data, time_scale = 1) {
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime) %>%
    group_by(bout) %>%
    arrange(StartTime, .by_group = TRUE) %>%
    mutate(
      CallTime   = StartTime * time_scale,
      WhaleIndex = as.integer(factor(WhaleID)),
      bout_id    = as.integer(factor(bout))
    ) %>%
    ungroup()
  
  unique_whales <- sort(unique(group_data$WhaleIndex))
  
  # Create all unique i ≠ j pairs (directional influence)
  pair_data <- expand.grid(i = unique_whales, j = unique_whales) %>%
    filter(i != j)
  
  stan_data <- list(
    N_events = nrow(group_data),
    N_whales = length(unique_whales),
    caller   = group_data$WhaleIndex,
    time     = group_data$CallTime,
    N_bouts  = length(unique(group_data$bout_id)),
    bout_id  = group_data$bout_id,
    N_pairs  = nrow(pair_data),
    pair_i   = pair_data$i,
    pair_j   = pair_data$j
  )
  
  return(list(data = stan_data, group_data = group_data, pair_data = pair_data))
}

# Subset o 2 bouts
prep_G1_small <- prep %>%
  filter(Group == "G1") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:2]) %>%
  ungroup()

data_G1_step2 <- prepare_discrete_phase_with_influence("G1", prep_G1_small, time_scale = 5)

# Full dataset on each group G1, G2, G3
# G1
data_G1_full <- prepare_discrete_phase_with_influence("G1", prep, time_scale = 5)

# G2
data_G2_full <- prepare_discrete_phase_with_influence("G2", prep, time_scale = 5)

# G3
data_G3_full <- prepare_discrete_phase_with_influence("G3", prep, time_scale = 5)

# ------------ Made changes to the stan model! -------------------------
# Check model compiling - Changed parameter chunk in model for omega, because i had omega[3],
# which was active with 34 calls, but seemingly couldn't be picked up by the model
# Whales’ baseline call tendencies (omega) are drawn from a group-level distribution with 
# mean mu_omega and spread sigma_omega.

# Gives robust estimates even for whales with unusual or sparse data,
# Pulls outliers like omega[3] closer to the center (but still allows flexibility)
# Dramatically improved convergence and ESS.

mod_step2 <- cmdstan_model("../models/discrete_phase_model.stan", cpp_options = list(stan_threads = TRUE))


fit_G1 <- mod_step2$sample(
  data = data_G1_full$data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 20,
  refresh = 5,
  init = function() list(
    omega = rep(1, data_G1_full$data$N_whales),
    A_raw = rep(0, data_G1_full$data$N_pairs),
    sigma = 0.1
  )
)
# === Save fit ===
# Full data on G1 with 4 chains 
fit_G1$save_object("fit_G1_pairwise_full.rds")


#And dramatically improves convergence and ESS.
fit_G1$summary() |>
  dplyr::filter(rhat > 1.01)

# 
fit_G1$summary() |>
  dplyr::filter(ess_bulk < 400 | ess_tail < 400)

# Great
posterior::as_draws_rvars(fit_G1$draws()) |>
  bayesplot::mcmc_trace(pars = c("omega[1]", "omega[2]", "sigma"))

# fixed omega[3] - posterior is better captured now. Previouls omega[3]
# caused issues due to: Whale 3’s calling being highly inconsistent or 
# reactive to others, resulting in a highly skewed posterior on its baseline rate.
posterior::as_draws_df(fit_G1$draws("omega[3]")) |>
  ggplot(aes(x = `omega[3]`)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "white") +
  theme_minimal() +
  labs(title = "Posterior of omega[3]")

data_G1_full$group_data |>
  dplyr::filter(WhaleIndex == 3) |>
  dplyr::summarise(n_calls = dplyr::n())


