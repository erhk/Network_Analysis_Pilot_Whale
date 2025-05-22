# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr)
library(ggplot2)

# Load and clean data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")
prep$CallType <- ifelse(is.na(prep$CallType), "InitiationCall", prep$CallType)

# Compile model
mod_hier <- cmdstan_model("../models/kuramoto_hierarchical.stan")

# Helper: prepare Stan data for a given group
prepare_group_data <- function(group_name, prep_data, time_scale = 1) {
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime) %>%
    group_by(bout) %>%
    arrange(StartTime, .by_group = TRUE) %>%
    mutate(
      CallIndex  = row_number(),
      WhaleIndex = as.integer(factor(WhaleID)),
      bout_id    = as.integer(factor(bout)),
      CallTime   = StartTime * time_scale  # ere we safely apply the scale
    ) %>%
    ungroup()
  
  unique_whales <- sort(unique(group_data$WhaleIndex))
  
  pair_data <- group_data %>%
    group_by(bout_id) %>%
    filter(n_distinct(WhaleIndex) >= 2) %>%
    summarise(Pairs = list(as.data.frame(t(combn(unique(WhaleIndex), 2)))), .groups = "drop") %>%
    unnest(Pairs) %>%
    rename(i = V1, j = V2)
  
  pair_data_rev <- pair_data %>% rename(i = j, j = i)
  
  pair_data <- bind_rows(pair_data, pair_data_rev) %>%
    filter(i != j) %>%         # Removes self-pairs
    distinct() %>%
    arrange(i, j)
  
  stan_data <- list(
    N_whales = length(unique_whales),
    N_events = nrow(group_data),
    caller   = group_data$WhaleIndex,
    time     = group_data$CallTime,        # Time now scaled
    N_bouts  = length(unique(group_data$bout_id)),
    bout_id  = group_data$bout_id,
    N_pairs  = nrow(pair_data),
    pair_i   = pair_data$i,
    pair_j   = pair_data$j
  )
  
  return(list(
    data = stan_data,
    group_data = group_data,
    pair_data = pair_data
  ))
}

# ----- Use this block to fit the hierarchical model on G2 small subset -----

# Subset G1 to 2 bouts
prep_G1_small <- prep %>%
  filter(Group == "G1") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:5]) %>%
  ungroup() %>%
  arrange(bout, StartTime)

# Scale time by 10x to boost theta evolution
data_G1_small <- prepare_group_data("G1", prep_G1_small, time_scale = 0.5)

save(data_G1_small, file = "data_G1_small.RData")

# sanity check that data is a-okay
range(data_G1_small$data$caller)      # Must be between 1 and N_whales
range(data_G1_small$data$pair_i)      # Also 1 to N_whales
range(data_G1_small$data$pair_j)
any(is.na(data_G1_small$data$time)) 


# Fit model
fit_G1_hier <- mod_hier$sample(
  data            = data_G1_small$data,
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
saveRDS(fit_G1_hier, "fit_kuramoto_G1_hierarchical.rds")
fit_G1_hier$save_object("fit_G1_hier_full.rds")




# Faster sampling version (plumbing check)
fit_light <- mod_hier$sample(
  data = data_G1_small$data,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 200,
  max_treedepth = 20,
  adapt_delta = 0.99,
  refresh = 10,
 
)
# Plumping check
saveRDS(fit_light, "fit_light_plumb.rds")
fit_light$save_object("fit_light_plumb.rds")
