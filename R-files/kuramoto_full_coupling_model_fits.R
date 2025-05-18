# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr) # for ucloud: install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))) 
library(ggplot2)
library(reshape2)

# --------- Load the data
prep <- read.delim("../Data/PrepData.csv", sep = ",")
prep$CallType <- ifelse(is.na(prep$CallType), "InitiationCall", prep$CallType)

# --------- Models
mod <- cmdstan_model("stan model/kuramoto_model.stan")#, cpp_options = list(stan_threads = TRUE))
mod_hier <- cmdstan_model("stan model/kuramoto_hierarchical.stan")

# --------- Prepare data function
prepare_group_data <- function(group_name, prep_data) {
  group_data <- prep %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime) %>%
    group_by(bout) %>%
    arrange(StartTime, .by_group = TRUE) %>%
    mutate(
      CallIndex = row_number(),
      CallTime = StartTime
    ) %>%
    ungroup() %>%
    mutate(
      WhaleIndex = as.integer(factor(WhaleID)),
      bout_id = as.integer(factor(bout))  # Re-index bouts 1:N
    )
  
  # Unique whale IDs in this group
  unique_whales <- sort(unique(group_data$WhaleIndex))
  
  # Build pairwise interactions only for bouts with 2+ whales
  pair_data <- group_data %>%
    group_by(bout_id) %>%
    filter(n_distinct(WhaleIndex) >= 2) %>%
    summarise(Pairs = list(as.data.frame(t(combn(unique(WhaleIndex), 2)))), .groups = "drop") %>%
    unnest(Pairs)
  
  # Rename V1 and V2 to i and j
  pair_data <- pair_data %>% rename(i = V1, j = V2)
  
  # Create reversed pairs
  pair_data_rev <- pair_data %>% rename(i = j, j = i)
  
  # Combine original, reversed, and remove duplicates
  pair_data <- bind_rows(pair_data, pair_data_rev) %>%
    distinct() %>%
    arrange(i, j)
  
  # Add self-pairs for all whales
  self_pairs <- group_data %>%
    distinct(WhaleIndex) %>%
    mutate(i = WhaleIndex, j = WhaleIndex) %>%
    select(i, j)
  
  pair_data <- bind_rows(pair_data, self_pairs) %>%
    distinct() %>%
    arrange(i, j)
  
  # Build Stan data list
  stan_data <- list(
    N_whales = length(unique_whales),
    N_events = nrow(group_data),
    caller = group_data$WhaleIndex,
    time = group_data$CallTime,
    N_bouts = length(unique(group_data$bout_id)),
    bout_id = group_data$bout_id,
    N_pairs = nrow(pair_data),
    pair_i = pair_data$i,
    pair_j = pair_data$j
  )
  
  return(list(
    data = stan_data,
    group_data = group_data,
    pair_data = pair_data
  ))
}

# Create subsets --------- Create a small subset (2 bouts from G1)
prep_small <- prep %>%
  filter(Group == "G1") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:2]) %>%
  ungroup() %>%
  arrange(bout, StartTime)

# Use prep_group_data on subset
data_G1_small <- prepare_group_data("G1", prep_small)
save(data_G1_small, file = "data_G1_small.RData")

# Check time makes sense! min = 0 
summary(diff(data_G1_small$data$time))


# kuramoto_pairwise model ------------------- Fit model on G1

# Fit full kuramoto model. It's heavy and struggles at higher settings
fit_kuramoto_full <- mod$sample(
  data = data_G1_small$data,
  chains = 1,
  iter_warmup = 250,
  iter_sampling = 250,
  refresh = 10,
  max_treedepth = 12,
  adapt_delta = 0.9,
  init = function() list(
    omega = rep(0, data_G1_small$data$N_whales),
    K = 0.1,
    A_raw = rep(0, data_G1_small$data$N_pairs),
    sigma = 0.05
  )
)

# Save model fit
saveRDS(fit_kuramoto_full, file = "fit_kuramoto_G1_small.rds")

# Check output - Currently it isn't amazing. 
fit_kuramoto_full$summary(variables = c("K", "sigma", "lp__"))  # Key indicators
fit_kuramoto_full$summary() %>%
  dplyr::filter(rhat > 1.01 | ess_bulk < 100)  # Look for trouble


# Kuramoto_hierarchical-------------------- Fit model on G1 - took 5.5 hour to run :(

fit_G2_hier <- mod_hier$sample(
  data            = data_G2_small$data,
  chains          = 1,
  iter_warmup     = 200,
  iter_sampling   = 300,
  adapt_delta     = 0.95,
  max_treedepth   = 15,
  refresh         = 20,
  init = function() list(
    omega = rep(0, data_G2_small$data$N_whales),
    K     = 0.05,
    z_A   = rep(0, data_G2_small$data$N_pairs),
    tau   = 0.05,
    sigma = 0.05
  )
)
saveRDS(fit_G2_hier, "fit_kuramoto_G2_hierarchical.rds")
load("data_G2_small.RData")
fit <- fit_G2_hier
source("model_diagnostics.R")
fit_G2_hier$save_object("fit_G2_hier_full.rds")




# Group 2 --------- Create a small subset (2 bouts from G2)
prep_small <- prep %>%
  filter(Group == "G2") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:2]) %>%
  ungroup() %>%
  arrange(bout, StartTime)

# Use prep_group_data on subset
data_G2_small <- prepare_group_data("G2", prep_small)
save(data_G2_small, file = "data_G2_small.RData")

# Check time makes sense! min = 0 
summary(diff(data_G2_small$data$time))

# Fit full kuramoto model. It's heavy and struggles at higher settings
fit_kuramoto_full_G2 <- mod$sample(
  data = data_G2_small$data,
  chains = 1,
  iter_warmup = 250,
  iter_sampling = 250,
  refresh = 10,
  max_treedepth = 12,
  adapt_delta = 0.9,
  init = function() list(
    omega = rep(0, data_G2_small$data$N_whales),
    K = 0.1,
    A_raw = rep(0, data_G2_small$data$N_pairs),
    sigma = 0.05
  )
)

# Save model fit
saveRDS(fit_kuramoto_full_G2, file = "fit_kuramoto_G2_small.rds")

# Check output - Currently it isn't amazing. 
fit_kuramoto_full_G2$summary(variables = c("K", "sigma", "lp__"))  # Key indicators
fit_kuramoto_full_G2$summary() %>%
  dplyr::filter(rhat > 1.01 | ess_bulk < 100)  # Look for trouble


