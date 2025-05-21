# Baby discrete omega model -  it runs
library(dplyr)
library(tidyr)
library(readr)
library(cmdstanr)
library(ggplot2)



# preprocess function
prepare_discrete_phase_data <- function(group_name, prep_data, time_scale = 1) {
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
  
  stan_data <- list(
    N_events = nrow(group_data),
    N_whales = length(unique(group_data$WhaleIndex)),
    caller   = group_data$WhaleIndex,
    time     = group_data$CallTime,
    N_bouts  = length(unique(group_data$bout_id)),
    bout_id  = group_data$bout_id
  )
  
  return(list(data = stan_data, group_data = group_data))
}


# Load and preprocess
prep <- read.delim("../../Data/PrepData.csv", sep = ",")
prep$CallType <- ifelse(is.na(prep$CallType), "InitiationCall", prep$CallType)

prep_small <- prep %>%
  filter(Group == "G1") %>%
  group_by(bout) %>%
  filter(bout %in% unique(bout)[1:2]) %>%
  ungroup()

data_G1_phase <- prepare_discrete_phase_data("G1", prep_small, time_scale = 5)

# Compile and fit
mod_omega_only <- cmdstan_model("../models/discrete_phase_model.stan")

fit_omega <- mod_omega_only$sample(
  data = data_G1_phase$data,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 300,
  refresh = 20,
  
)
