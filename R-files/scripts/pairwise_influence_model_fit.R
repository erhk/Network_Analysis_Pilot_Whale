# Discrete Phase Influence Model - Adding pairwise influence

pacman::p_load(dplyr, tidyr, readr, cmdstanr, ggplot2, posterior)


# Load data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")


# # # ----------------- Preprocess function --------------------- # # #
 
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

# Subset to 2 bouts
# prep_G1_small <- prep %>%
#   filter(Group == "G1") %>%
#   group_by(bout) %>%
#   filter(bout %in% unique(bout)[1:2]) %>%
#   ungroup()
# 
# data_G1_step2 <- prepare_discrete_phase_with_influence("G1", prep_G1_small, time_scale = 5)

# Preporcess full dataset on each group G1, G2, G3
# G1
data_G1_full <- prepare_discrete_phase_with_influence("G1", prep, time_scale = 5)

# G2
data_G2_full <- prepare_discrete_phase_with_influence("G2", prep, time_scale = 5)

# G3
data_G3_full <- prepare_discrete_phase_with_influence("G3", prep, time_scale = 5)



# # # ----------------- Made changes to the stan model! --------------------- # # #

# --- model compiling issues

# Group 1 issues
# Changed parameter chunk in model for omega, because i had omega[3],
# which was active with 34 calls, but seemingly couldn't be picked up by the model
# Wanted to change the whales’ baseline call tendencies (omega), drawn from a group-level distribution with 
# mean mu_omega and spread sigma_omega. This however impaired model's ability to sample from A[i,j] matrix. 
# I've left omega as previous, and will accept this as an instability in the fitted model, likely due to many calling/overlapping
# with other members in bouts.


mod_step2 <- cmdstan_model("../models/discrete_phase_model_pairwise_influence.stan", cpp_options = list(stan_threads = TRUE))


fit_G2 <- mod_step2$sample(
  data = data_G2_full$data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 20,
  refresh = 5,
  init = function() list(
    omega = rep(1, data_G2_full$data$N_whales),
    A_raw = rep(0, data_G2_full$data$N_pairs),
    sigma = 0.1
  )
)
# Save fit 
fit_G1_fix$save_object("fit_G1_pairwise_fix.rds") 
fit_G2$save_object("fit_G2_pairwise.rds")
fit_G3$save_object("fit_G3_pairwise.rds")


# # # ----------------- Check Model convergence and fit ----------------- # # #
fit_G1_fix <- readRDS("fit_G1_pairwise_fix.rds")
fit_G2 <- readRDS("fit_G2_pairwise.rds")
fit_G3 <- readRDS("fit_G3_pairwise.rds")

# check rhat - ideally all above 1
fit_G3$summary() |>
  dplyr::filter(rhat > 1.01)

# check ess - ideally below 400
fit_G3$summary() |>
  dplyr::filter(ess_bulk < 400 | ess_tail < 400)

# diagnostics summary
fit_G3$diagnostic_summary()


# Check chains
posterior::as_draws_rvars(fit_G2$draws()) |>
  bayesplot::mcmc_trace(pars = c("omega[1]", "omega[2]", "sigma"))

posterior::as_draws_rvars((fit_G2$draws())) |>
  bayesplot::mcmc_trace(pars = c("A_raw[5]", "A_raw[10]", "A_raw[20]"))



# # # ----------------- Analysis and visualisation ----------------- # # #

# Create a lookup table from WhaleIndex to WhaleID
# Get WhaleIndex-to-WhaleID from the original preprocessing input
index_map <- data_G1_full$group_data |>
  distinct(WhaleIndex, WhaleID) |>
  mutate(WhaleID = as.character(WhaleID))

# Create omega summary
omega_summary <- fit_G1$summary(variables = "omega") |>
  mutate(whale = as.integer(sub("omega\\[", "", sub("\\]", "", variable))))

# Omega match whaleid to index. Omega is a single vector - mapped to index
# omega[i]	i (1D index)	Single join: whale = WhaleIndex
omega_summary_named <- omega_summary |> 
  filter(!is.na(whale)) |> 
  left_join(index_map, by = c("whale" = "WhaleIndex"))

# Confirm all matched
table(is.na(omega_summary_named$WhaleID))

library(posterior)
library(ggplot2)
library(reshape2)

variables <- fit_G1$variables()
print(variables)

# Extract A_raw posterior means
a_draws <- as_draws_df(fit_G1_fix$draws(variables = "A_raw"))
a_means <- colMeans(a_draws)

# Reconstruct A matrix
N <- data_G1_full$data$N_whales
A_matrix <- matrix(0, N, N)
for (k in seq_along(a_means)) {
  i <- data_G1_full$data$pair_i[k]
  j <- data_G1_full$data$pair_j[k]
  A_matrix[i, j] <- a_means[k]
}

# Map whale indices to WhaleIDs
id_lookup <- omega_summary_named |>
  select(whale, WhaleID) |>
  arrange(whale)

rownames(A_matrix) <- id_lookup$WhaleID
colnames(A_matrix) <- id_lookup$WhaleID

# Heatmap
A_df <- melt(A_matrix, varnames = c("Target", "Source"), value.name = "Influence")

ggplot(A_df, aes(x = Source, y = Target, fill = Influence)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(title = "Whale Influence Matrix (A[i,j])", x = "Calling Whale", y = "Influenced Whale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Extract posterior draws
omega_draws <- as_draws_df(fit_G1$draws(variables = "omega"))

ggplot(omega_summary_named, aes(x = reorder(WhaleID, -mean), y = mean)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.2) +
  labs(
    title = "Intrinsic Calling Rates by WhaleID (omega[i])",
    x = "Whale ID",
    y = "Estimated Rate (1 / Expected Delay)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
