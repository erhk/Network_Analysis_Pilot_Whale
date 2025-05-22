# Prepare data for response-focused model
# Assumes `prep` is your full call dataset

library(dplyr)
library(tidyr)

# Step 1: Identify Initiation → Response pairs within bouts
response_pairs <- prep %>%
  filter(CallType %in% c("InitiationCall", "ResponseCall")) %>%
  arrange(Group, bout, StartTime) %>%
  group_by(Group, bout) %>%
  mutate(
    InitiatorID = first(WhaleID[CallType == "InitiationCall"]),
    InitiatorTime = first(StartTime[CallType == "InitiationCall"])
  ) %>%
  filter(CallType == "ResponseCall") %>%
  mutate(
    Latency = StartTime - InitiatorTime,
    Initiator = InitiatorID,
    Responder = WhaleID
  ) %>%
  ungroup()

# Step 2: Convert WhaleIDs to indices
whales <- sort(unique(c(response_pairs$Responder, response_pairs$Initiator)))
whale_index <- setNames(seq_along(whales), whales)

response_pairs <- response_pairs %>%
  mutate(
    ResponderIndex = whale_index[Responder],
    InitiatorIndex = whale_index[Initiator]
  )

# Step 3: Build Stan data list
stan_response_data <- list(
  N = nrow(response_pairs),
  latency = response_pairs$Latency,
  N_whales = length(whales),
  responder = response_pairs$ResponderIndex,
  initiator = response_pairs$InitiatorIndex
)

# Save for fitting
saveRDS(stan_response_data, file = "response_model_data.rds")
saveRDS(response_pairs, file = "response_pairs_longform.rds")


# Fit model - full

# Load the prepared data
stan_response_data <- readRDS("response_model_data.rds")

# Compile the response model
mod_response <- cmdstan_model("stan model/response_model.stan", cpp_options = list(stan_threads = TRUE))

# Fit it - simple
fit_response <- mod_response$sample(
  data = stan_response_data,
  chains = 1,
  iter_warmup = 250,
  iter_sampling = 250,
  refresh = 1
)

# Fit it - complex
fit_response_complexer <- mod_response$sample(
  data = stan_response_data,
  chains = 2,
  parallel_chains = 2,
  threads_per_chain = 2,
  iter_warmup = 1000,
  iter_sampling = 2000,
  refresh = 10,
  max_treedepth = 20,
  adapt_delta = 0.99
)
# Save for diagnostics and plotting
saveRDS(fit_response_complexer, "fit_response_model_compx.rds")
fit_comx <- readRDS("fit_response_model_compx.rds")
fit_comx$summary(variables = c("omega", "alpha", "sigma"))

posterior_df <- fit_comx$draws(format = "draws_df")


# Save for diagnostics and plotting
saveRDS(fit_response, "fit_response_model.rds")

fit <- readRDS("fit_response_model.rds")
fit$summary(variables = c("omega", "alpha", "sigma"))

posterior_df <- fit$draws(format = "draws_df")

library(dplyr)
library(ggplot2)

alpha_df <- posterior_df %>%
  select(starts_with("alpha[")) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = "param", values_to = "mean") %>%
  mutate(whale_index = as.integer(gsub("alpha\\[|\\]", "", param)))

# Map index back to WhaleID
response_pairs <- readRDS("response_pairs_longform.rds")
id_lookup <- response_pairs %>%
  select(Initiator, InitiatorIndex) %>%
  distinct() %>%
  arrange(InitiatorIndex)

alpha_df <- alpha_df %>%
  left_join(id_lookup, by = c("whale_index" = "InitiatorIndex"))

# Plot
ggplot(alpha_df, aes(x = reorder(Initiator, mean), y = mean)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Initiator Influence (alpha)", x = "Initiator Whale", y = "Mean Influence on Response Speed")
