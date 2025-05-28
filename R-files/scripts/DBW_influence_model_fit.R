# Discrete bout-wide (cumulative) Influence Model

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

# Generate influence pairs (dense by default)
pairs <- expand.grid(i = whales, j = whales) %>%
  filter(i != j)  # skip self-influence, or keep if you want A[i,i]

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
  data = stan_data_g1,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.999,
  max_treedepth = 20,
  refresh = 5,
  seed = 1990
)
#fit$save_object("fit_G1_dpi_cum.rds") 
#fit$save_object("fit_G2_dpi_cum.rds")
#fit$save_object("fit_G3_dpi_cum.rds")
fit$save_object("fit_G2_lambda_normprior_all.rds")
fit$save_object("fit_G1_lambda_normprior_all.rds")


#### -------- fitting model subsets to fix lambda ---- ####
# 
# # 1. Select the first 5 bouts in order
# first_bouts_g2 <- g2_data %>%
#   arrange(StartTime) %>%
#   distinct(bout) %>%
#   slice(1:5) %>%
#   pull(bout)
# 
# # 2. Filter G2 data to just those bouts
# g2_subset_df <- g2_data %>%
#   filter(bout %in% first_bouts_g2) %>%
#   arrange(StartTime) %>%
#   mutate(
#     WhaleIndex = as.integer(factor(WhaleID)),
#     BoutIndex = as.integer(factor(bout))
#   )
# 
# # 3. Define pair_i and pair_j (all possible caller pairs present in data)
# whale_ids <- unique(g2_subset_df$WhaleIndex)
# pair_grid <- expand.grid(i = whale_ids, j = whale_ids) %>%
#   filter(i != j)  # exclude self-influence if not modeling A_ii
# 
# # 4. Build Stan data list
# stan_data_g2_subset <- list(
#   N_events = nrow(g2_subset_df),
#   N_whales = length(unique(g2_subset_df$WhaleIndex)),
#   caller = g2_subset_df$WhaleIndex,
#   time = g2_subset_df$StartTime,
#   N_bouts = length(unique(g2_subset_df$BoutIndex)),
#   bout_id = g2_subset_df$BoutIndex,
#   N_pairs = nrow(pair_grid),
#   pair_i = pair_grid$i,
#   pair_j = pair_grid$j
# )
# 
# 
# whale_id_map_g2_subset <- g2_subset_df %>%
#   select(WhaleID, WhaleIndex) %>%
#   distinct() %>%
#   arrange(WhaleIndex)
# 
# mod <- cmdstan_model("../models/discrete bout-wide Influence model.stan", cpp_options = list(stan_threads = TRUE))
# 
# fit <- mod$sample(
#   data = stan_data_g2_subset,
#   chains = 4,
#   parallel_chains = 4,
#   threads_per_chain = 4,
#   iter_warmup = 2000,
#   iter_sampling = 2000,
#   adapt_delta = 0.999,
#   max_treedepth = 20,
#   refresh = 100,
#   seed = 1990
# )
#### --- ####


# ---- Plot Posterior for lambda 
fit_g2 <- readRDS("fit_G2_lambda_normprior_all.rds")
fit <- readRDS("fit_G2_lambda_all.rds")


# Extract posterior draws for lambda
lambda_draws <- fit$draws("lambda", format = "draws_df") %>%
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

#ggsave("lambda_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# Diagnostics for lambda because it looks very strange, bimodal, with two peaks near 0.2 and 1.7, valley close to 1
fit$summary("lambda")
fit$cmdstan_diagnose()
posterior::rhat(fit$draws("lambda"))

mcmc_trace(fit$draws("lambda"))

#### --------------- Model Checks ------------------ ####

# check rhat - ideally all above 1
fit$summary() |>
  dplyr::filter(rhat > 1.01)

# check ess - ideally below 400
fit$summary() |>
  dplyr::filter(ess_bulk < 400 | ess_tail < 400)

# diagnostics summary
fit$diagnostic_summary()



#### ------ !!! Junk below. Just to get a quick overview - Visualising in different script !!! ---------------- ####



# # Permanently map whaleid back onto whaleindex, so we have their actual identities again
# #g1_data <- g1_data %>%
# #  mutate(WhaleIndex = WhaleID)
# 
# # Temporary remapping, but it's kinda annoying
# whale_id_map <- g1_data %>%
#   select(WhaleIndex, WhaleID) %>%
#   distinct()
# 
# # Draws
# draws <- as_draws_df(fit$draws())
# 
# # pairs
# pairs <- tibble(
#   i = stan_data_g1$pair_i,
#   j = stan_data_g1$pair_j
# )
# 
# # Get draws of A_raw with pair labels
# # Get posterior draws of A_raw with i,j labels
# draws_A <- fit$draws("A_raw") %>%
#   spread_draws(A_raw[k]) %>%
#   mutate(
#     i = pairs$i[k],
#     j = pairs$j[k],
#     pair = paste0("A[", i, ",", j, "]")
#   )
# 
# # Summarize the posterior
# summary_A <- draws_A %>%
#   group_by(i, j, pair) %>%
#   median_qi(A_raw)
# 
# summary_A_named <- summary_A %>%
#   left_join(whale_map, by = c("i" = "WhaleIndex")) %>%
#   rename(i_label = WhaleID) %>%
#   left_join(whale_map, by = c("j" = "WhaleIndex")) %>%
#   rename(j_label = WhaleID) %>%
#   mutate(pair = paste0("A[", i_label, ", ", j_label, "]"))
# 
# 
# 
# # Example: extract full draws of A_raw[k] with labels
# draws_A_full <- fit$draws("A_raw") %>%
#   spread_draws(A_raw[k]) %>%
#   mutate(
#     i = pairs$i[k],
#     j = pairs$j[k]
#   ) %>%
#   left_join(whale_map, by = c("i" = "WhaleIndex")) %>%
#   rename(i_label = WhaleID) %>%
#   left_join(whale_map, by = c("j" = "WhaleIndex")) %>%
#   rename(j_label = WhaleID) %>%
#   mutate(pair = paste0(j_label, " → ", i_label))
# 
# # Choose a few pairs to plot
# pairs_to_plot <- c("F → E", "E → F", "G → E")
# 
# draws_A_full %>%
#   filter(pair %in% pairs_to_plot) %>%
#   ggplot(aes(x = A_raw, fill = pair)) +
#   geom_density(alpha = 0.6) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   facet_wrap(~pair, scales = "free") +
#   labs(title = "Posterior Distributions of A[i,j] Influence",
#        x = "Influence Value", y = "Density") +
#   theme_minimal()
# 
# 
# 
# 
# 
# ### --------- Heatmap
# ggplot(summary_A_named, aes(x = reorder(pair, -A_raw), y = A_raw)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
#   coord_flip() +
#   labs(title = "Posterior estimates of A[i,j]", y = "Influence", x = "")
# 
# 
# 
# ggplot(summary_A_named, aes(x = j_label, y = i_label, fill = A_raw)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limits = c(min(summary_A_named$A_raw), max(summary_A_named$A_raw)),
#                        name = "Influence") +
#   labs(title = "Whale Influence Matrix (A[i,j])", x = "Influencer (j)", y = "Receiver (i)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ### --------- Network
# 
# # Step 1: Build the vertex (node) list correctly
# nodes <- tibble(name = unique(c(summary_A_named$i_label, summary_A_named$j_label)))
# 
# # Step 2: Build the graph
# g <- summary_A_named %>%
#   filter(!is.na(A_raw)) %>%
#   select(from = j_label, to = i_label, weight = A_raw) %>%
#   graph_from_data_frame(vertices = nodes, directed = TRUE)
# 
# # Add edge weights (median influence)
# E(g)$weight <- summary_A_named$A_raw
# E(g)$color <- ifelse(E(g)$weight > 0, "red", "blue")
# 
# 
# # Plot with ggraph
# ggraph(g, layout = "circle") +
#   geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight), edge_color = weight),
#                  arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
#   geom_node_point(size = 6) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 5) +
#   scale_edge_color_gradient2(low = "blue", mid = "gray80", high = "red", midpoint = 0) +
#   labs(title = "Directed Influence Graph (Whale Calls)") +
#   theme_void()
# 
