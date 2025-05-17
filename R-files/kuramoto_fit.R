# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr) # for ucloud: install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))) 
library(ggplot2)
library(reshape2)

# Load the data
prep <- read.delim("../Data/PrepData.csv", sep = ",")
prep$CallType <- ifelse(is.na(prep$CallType), "InitiationCall", prep$CallType)

# Model
mod <- cmdstan_model("stan model/kuramoto_model.stan")#, cpp_options = list(stan_threads = TRUE))

# Prepare data function
prepare_group_data <- function(group_name, prep_data) {
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, Latency) %>%
    group_by(bout) %>%
    mutate(
      CallIndex = row_number(),
      CallTime = cumsum(replace_na(Latency, 0))
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




# --------- Prep Small

# Try only the first 2 bouts
prep_small <- prep %>%
  filter(Group == "G1") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:2]) %>%
  ungroup()

data_G1_small <- prepare_group_data("G1", prep_small)


system.time({
  fit_G1_small <- mod$sample(
    data = data_G1_small$data,
    chains = 1,
    iter_warmup = 250,
    iter_sampling = 250,
    max_treedepth = 10,
    adapt_delta = 0.9,
    init = function() list(
      omega = rep(0, data_G1_small$data$N_whales),
      K = 0.1,
      A_raw = rep(0, data_G1_small$data$N_pairs),
      sigma = 0.5
    )
  )
})

fit_G1_small$summary(variables = c("K", "sigma", "lp__"))

fit_G1_small$summary() %>%
  dplyr::filter(rhat > 1.01 | ess_bulk < 100)



# --------- Prep Big - too big to run at the moment ----------------- #


# data_G1 <- prepare_group_data("G1", prep)
# data_G2 <- prepare_group_data("G2", prep)
# data_G3 <- prepare_group_data("G3", prep)


# Running G1 on computer
# fit_G1 <- mod$sample(
#   data = data_G2$data,
#   chains = 1,
#   iter_warmup = 1000,
#   iter_sampling = 1000,
#   max_treedepth = 20,
#   adapt_delta = 0.95,
#   init = function() list(
#     omega = rep(0, data_G1$data$N_whales),
#     K = 0.1,
#     A_raw = rep(0, data_G1$data$N_pairs),
#     sigma = 0.5
#   )
# )



#fit_G1$cmdstan_diagnose()
#fit_G1$summary() %>% dplyr::filter(rhat > 1.01 | ess_bulk < 400)
