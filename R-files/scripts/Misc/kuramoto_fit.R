# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr) # for ucloud: install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))) 
library(ggplot2)
library(reshape2)
library(purrr)


# Load the data
prep_data <- read.delim("../../data/PrepData.csv", sep = ",")
prep_data$CallType <- ifelse(is.na(prep_data$CallType), "InitiationCall", prep_data$CallType)



prepare_group_data_batch <- function(group_name, prep_data) {
  num_bouts <- 5
  
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime)
  
  first_n_bouts <- group_data %>%
    distinct(bout) %>%
    slice_head(n = num_bouts) %>%
    pull(bout)
  
  group_data <- group_data %>%
    filter(bout %in% first_n_bouts)
  
  
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime) 
    # drop small bouts
    min_whales <- 2
    min_events <- 3
  
    group_data <- prep_data %>%
      filter(Group == group_name) %>%
      arrange(bout, StartTime)
    
    # Drop small bouts
    valid_bouts <- group_data %>%
      group_by(bout) %>%
      summarise(
        n_whales = n_distinct(WhaleID),
        n_events = n(),
        .groups = "drop"
      ) %>%
      filter(n_whales >= min_whales, n_events >= min_events) %>%
      pull(bout)
    
    group_data <- group_data %>%
      filter(bout %in% valid_bouts) %>%
      group_by(bout, StartTime) %>%
      mutate(time_group = cur_group_id()) %>%
      ungroup() %>%
      mutate(
        CallTime = StartTime,
        WhaleIndex = as.integer(factor(WhaleID)),
        bout_id = as.integer(factor(bout))
      )
    
  # Map to time[g] and bout_id[g] in Stan.
  time_by_group <- group_data %>%
    group_by(time_group) %>%
    summarise(
      time = first(CallTime),
      bout_id = first(bout_id),
      .groups = "drop"
    )
  
  event_group_map <- group_data %>%
    mutate(event_id = row_number()) %>%
    group_by(time_group) %>%
    summarise(events = list(event_id), .groups = "drop")
  
  # Collects all events (calls) in each time_group.
  # These populate event_matrix[g, k] and event_counts[g], which allow the model 
  # to process all call events at the same time step, needed for Kuramoto batch updates.
  max_events <- max(map_int(event_group_map$events, length))
  event_matrix <- matrix(0, nrow = nrow(event_group_map), ncol = max_events)
  event_counts <- integer(nrow(event_group_map))
  
  for (i in seq_len(nrow(event_group_map))) {
    idxs <- event_group_map$events[[i]]
    event_matrix[i, seq_along(idxs)] <- idxs
    event_counts[i] <- length(idxs)
  }
  # Pairwise influence structure
  # For each bout, generate all possible pairs of whales that co-occur, which represent possible directed influence.
  pair_data <- group_data %>%
    group_by(bout_id) %>%
    filter(n_distinct(WhaleIndex) >= 2) %>%
    summarise(Pairs = list(as.data.frame(t(combn(unique(WhaleIndex), 2)))), .groups = "drop") %>%
    unnest(Pairs) %>%
    rename(i = V1, j = V2)
  
  # Reverse direction in pairs. Might need to remove for computational time
#  pair_data_rev <- pair_data %>% rename(i = j, j = i)
  
  # Pairs final, forward + reverse. No self-pair
 # pair_data <- bind_rows(pair_data, pair_data_rev) %>%
  #  distinct() %>%
   # arrange(i, j)
  
  # Stan list 
  stan_data <- list(
    N_whales = length(unique(group_data$WhaleIndex)),
    N_events = nrow(group_data),
    caller = group_data$WhaleIndex,
    time = time_by_group$time,
    bout_id = time_by_group$bout_id,
    time_group = group_data$time_group,
    N_time_groups = nrow(time_by_group),
    N_pairs = nrow(pair_data),
    pair_i = pair_data$i,
    pair_j = pair_data$j,
    event_matrix = event_matrix,
    event_counts = event_counts
  )
  # Add max_events for Stan array dimension declaration
  stan_data$max_events <- max(event_counts)
  
  return(list(
    data = stan_data,
    group_data = group_data,
    pair_data = pair_data
  ))
}



# === Run for Group G1 ===
batch_data <- prepare_group_data_batch("G1", prep_data)
# === Run for Group G2 ===
batch_data_2 <- prepare_group_data_batch("G2", prep_data)
# === Run for Group G3 ===
batch_data_3 <- prepare_group_data_batch("G3", prep_data)


# === Save prepped data ===
saveRDS(batch_data, file = "data_batch_G1.rds")
data_batch_G1 <- readRDS("data_batch_G1.rds")

saveRDS(batch_data_3, file = "data_batch_G2.rds")
data_batch_G3 <- readRDS("data_batch_G2.rds")

saveRDS(batch_data_3, file = "data_batch_G3.rds")
data_batch_G3 <- readRDS("data_batch_G3.rds")
# === Load model ===
mod <- cmdstan_model("../models/kuramoto_model.stan")#, cpp_options = list(stan_threads = TRUE))

# === Fit model ===
fit_G1 <- mod$sample(
  data = data_batch_G1$data,
  seed = 123,
  chains = 2,
  parallel_chains = 2,
  #threads_per_chain = 2,
  iter_warmup =500,
  iter_sampling = 500,
  refresh = 1,
  max_treedepth = 20,
  adapt_delta = 0.99,
  init = 0.5
  
  
)

# === Save fit ===
fit_G1$save_object("fit_batch_G1.rds")
fit_G2$save_object("fit_batch_G2.rds")
fit_G3$save_object("fit_batch_G3.rds")

