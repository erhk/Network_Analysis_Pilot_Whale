library(tidyverse)
library(cmdstanr)

# Load your dataset
df <- read.delim("../../Data/PrepData.csv", sep = ",")

# Function to preprocess data for one group
prep_stan_data <- function(df_group) {
  df_group <- df_group %>%
    mutate(
      WhaleNum = as.integer(factor(WhaleID)),
      BoutNum = as.integer(factor(bout))
    ) %>%
    group_by(bout) %>%
    mutate(TimeRel = StartTime - min(StartTime)) %>%
    ungroup() %>%
    arrange(bout, TimeRel)
  
  list(
    N_whales = n_distinct(df_group$WhaleNum),
    N_events = nrow(df_group),
    caller = df_group$WhaleNum,
    time = df_group$TimeRel,
    N_bouts = n_distinct(df_group$BoutNum),
    bout_id = df_group$BoutNum
  )
}

# Compile the model
mod <- cmdstan_model("../models/kuramoto_shared_coupling.stan")

# Fit the model to each group
df_g1 <- df %>% filter(Group == "G1")
stan_data_g1 <- prep_stan_data(df_g1)
saveRDS(stan_data_g1, file = "stan_data_g1.rds")

fit_g1 <- mod$sample(
  data = stan_data_g1,
  seed = 123,
  chains = 1,
  #parallel_chains = 2,
  #threads_per_chain = 2,
  iter_warmup =500,
  iter_sampling = 500,
  refresh = 1,
  max_treedepth = 20,
  adapt_delta = 0.99,
  init = 0.5
  
  
)


fit_g1$save_object("fit_sharedcoupl_G1.rds")



