# Analysis and visualisation of discrete bout-wide (cumulative) influence model: g2, G2, G3

# Add all the things
pacman::p_load(tidyverse, cmdstanr, tidybayes, dplyr, ggplot2, igraph, ggraph, bayesplot, ggrepel)  
               


#### --------------- Load all the needed files ---------------- ####
# Load fitted models
fit_g2 <- readRDS("../fits/cumulative_influence/fit_g2_dpi_cum.rds")
fit_g2 <- readRDS("../fits/cumulative_influence/fit_G2_dpi_cum.rds")
fit_g3 <- readRDS("../fits/cumulative_influence/fit_G3_dpi_cum.rds")

# Load stan data list
stan_data_g2 <- readRDS("saved_rds/g2_stan_data.rds")
stan_data_g2 <- readRDS("saved_rds/g2_stan_data.rds")
stan_data_g3 <- readRDS("saved_rds/g3_stan_data.rds")

# Load subsetted data
g2_data <- readRDS("saved_rds/g2_data.rds")
g2_data <- readRDS("saved_rds/g2_data.rds")
g3_data<- readRDS("saved_rds/g3_data.rds")

# Id map
whale_id_map_g2 <- readRDS("saved_rds/g2_whale_id_map.rds")
whale_id_map_g2 <- readRDS("saved_rds/g2_whale_id_map.rds")
whale_id_map_g3 <- readRDS("saved_rds/g3_whale_id_map.rds")

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g2 <- setNames(whale_id_map_g2$WhaleID, whale_id_map_g2$WhaleIndex)

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g2 <- setNames(whale_id_map_g2$WhaleID, whale_id_map_g2$WhaleIndex)

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g3 <- setNames(whale_id_map_g3$WhaleID, whale_id_map_g3$WhaleIndex)






#### --------------- Group 1 ---------------- ####

# ---- Posterior distribution of baseline calling rates ωᵢ
# Extract posterior draws for omega[i]
omega_plot_df_g2 <- fit_g2$draws("omega") %>%
  spread_draws(omega[i]) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])  # map index -> ID

# Plot posterior densities. They are bunched up at close to zero, lower bound for omega was set at 0.2
# That constraint is being hit by the posterior, meaning your data favor low intrinsic rates
# Suggests most of these whales are not initiating calls independently, they're socially reactive

ggplot(omega_plot_df_g2, aes(x = omega, fill = Whale)) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "red")+  # shows Stan’s lower bound
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Group 1 - Posterior of Baseline Calling Rates (ωᵢ)",
    x = expression(omega[i]),
    y = "Density",
    fill = "Whale"
  ) +
  theme_minimal(base_size = 14)

ggsave("omega_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# ---- Overlay influence and omega
omega_summary_g2 <- fit_g2$draws("omega") %>%
  spread_draws(omega[i]) %>%
  group_by(i) %>%
  summarise(omega_median = median(omega)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

A_summary_g2 <- fit_g2$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g2$pair_i[k], j = stan_data_g2$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Sum influence received by each whale (i.e., incoming influence to i)
incoming_influence_g2 <- A_summary_g2 %>%
  group_by(i) %>%
  summarise(total_influence_g2 = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

plot_overlay_g2 <- omega_summary_g2 %>%
  left_join(incoming_influence_g2, by = c("i", "Whale"))

# Added colour by whale, and icon by subgroup
plot_overlay_g2_col <- plot_overlay_g2 %>%
  left_join(g2_df %>% select(WhaleID, Source_subgroup) %>% distinct(), 
            by = c("Whale" = "WhaleID"))

ggplot(plot_overlay_g2_col, aes(x = total_influence_g2, y = omega_median)) +
  geom_point(aes(color = as.factor(Source_subgroup))) +
  geom_text_repel(aes(label = Whale), max.overlaps = 10, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Whales' Baseline Calling Rate vs. Social Influence (Group g2)",
    legend = element_text("Subgroup"),
    x = "Social influence: Sum of all incoming Aᵢⱼ values",
    y = expression("Median "*omega[i]*" (Baseline Calling Rate)"),
    color = "Subgroup"
  ) +
  theme_minimal(base_size = 14)

ggsave("omega_influce_overlay_g2.png", width = 10, height = 10, units = "in", dpi = 300)


# ---- Plot Posterior for lambda 

# Extract posterior draws for lambda
lambda_draws <- fit_g2$draws("lambda", format = "draws_df") %>%
  spread_draws(lambda)

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

ggsave("lambda_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# ---- Plot influence accross specific bouts, longer than 5 calls

long_bouts_g2 <- g2_df %>%
  group_by(bout) %>%
  filter(n() > 5) %>%
  pull(bout) %>%
  unique()


# Step: summarize A[i,j]
A_matrix <- fit_g2$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g2$pair_i[k], j = stan_data_g2$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Convert to matrix for fast lookup
N <- max(stan_data_g2$caller)
A_mat <- matrix(0, N, N)
for (row in 1:nrow(A_matrix)) {
  A_mat[A_matrix$i[row], A_matrix$j[row]] <- A_matrix$mean_A[row]
}

lambda_mean <- summarise_draws(fit_g2$draws("lambda"))$mean

lambda <- lambda_mean

influence_traces <- list()

for (b in long_bouts_g2[1:10]) {
  sub_bout <- g2_df %>%
    filter(bout == b) %>%
    arrange(StartTime) %>%
    mutate(event = row_number())
  
  n_events <- nrow(sub_bout)
  
  traces <- tibble()
  
  for (n in 2:n_events) {
    i <- sub_bout$WhaleIndex[n]
    t_n <- sub_bout$StartTime[n]
    
    influence_sum <- 0
    
    for (m in 1:(n - 1)) {
      j <- sub_bout$WhaleIndex[m]
      t_m <- sub_bout$StartTime[m]
      delta_t <- t_n - t_m
      decay <- exp(-lambda * delta_t)
      
      influence_sum <- influence_sum + A_mat[i, j] * decay
    }
    
    traces <- bind_rows(traces, tibble(
      Bout = b,
      Event = n,
      Time = t_n,
      Whale = sub_bout$WhaleID[n],
      Influence = influence_sum
    ))
  }
  
  influence_traces[[as.character(b)]] <- traces
}

influence_df <- bind_rows(influence_traces)
call_lines_df <- g2_df %>%
  filter(bout %in% long_bouts_g2[1:10]) %>%
  select(Bout = bout, Time = StartTime)


ggplot(influence_df, aes(x = Time, y = Influence, color = Whale)) +
  #geom_step(direction = "hv") +
  geom_point(aes(color = Whale, size = abs(Influence))) +
  geom_vline(data = call_lines_df, aes(xintercept = Time), color = "black", linetype = "dotted", size = 0.4) +

  facet_wrap(~ Bout, scales = "free_x") +
  labs(
    title = "Cumulative Influence Across Longer Bouts",
    x = "Time within Bout",
    y = "Summed Influence from Previous Calls"
  ) +
  theme_minimal()

# ---- Network graph


# Create a tidy edge list from A[i,j]
edge_df_g2 <- A_summary_g2 %>%
  mutate(
    from_whale = whale_index_to_id_g2[as.character(j)],
    to_whale   = whale_index_to_id_g2[as.character(i)]
  ) %>%
  select(from = from_whale, to = to_whale, mean_A)
# Optional threshold

# Create graph object
g2_influence <- graph_from_data_frame(edge_df_g2, directed = TRUE)


# Compute strength (weighted degree)
V(g2_influence)$in_strength  <- strength(g2_influence, mode = "in", weights = E(g2_influence)$mean_A)
V(g2_influence)$out_strength <- strength(g2_influence, mode = "out", weights = E(g2_influence)$mean_A)

# Centrality
V(g2_influence)$eigen_centrality <- eigen_centrality(g2_influence, weights = E(g2_influence)$mean_A)$vector
V(g2_influence)$betweenness <- betweenness(g2_influence, directed = TRUE)


ggraph(g2_influence, layout = "nicely") +
  geom_edge_fan(aes(width = abs(mean_A), color = mean_A),
                 arrow = arrow(length = unit(3, "mm")), end_cap = circle(3, "mm")) +
  geom_node_point(aes(size = in_strength)) +
  geom_node_label(aes(label = name)) +
  scale_edge_color_gradient2(low = "red", mid = "gray80", high = "green", midpoint = 0) +
  labs(title = "Whale Influence Network (from Aᵢⱼ)",
       subtitle = "Edge color = influence direction, width = strength, node size = total received influence") +
  theme_void()

### -- Add node strength to graph

# Step 1: Format edge list for plotting
edges <- A_summary_g2 %>%
  mutate(
    Source = whale_index_to_id_g2[as.character(j)],  # j = influencer
    Target = whale_index_to_id_g2[as.character(i)],  # i = influenced
    Influence = mean_A
  ) %>%
  filter(abs(Influence) > 0.01) %>% # Filter close to zero influence
  select(Source, Target, Influence)

# Graph and layout
g <- graph_from_data_frame(edges, directed = TRUE)

layout <- create_layout(g, layout = 'linear', circular = TRUE)

edges <- edges %>%
  mutate(
    from_x = layout$x[match(Source, layout$name)],
    from_y = layout$y[match(Source, layout$name)],
    to_x   = layout$x[match(Target, layout$name)],
    to_y   = layout$y[match(Target, layout$name)],
    dx = to_x - from_x,
    dy = to_y - from_y,
    norm = sqrt(dx^2 + dy^2),
    offset_x = -dy / norm * 0.05,
    offset_y = dx / norm * 0.05,
    mid_x = (from_x + to_x) / 2 + offset_x,
    mid_y = (from_y + to_y) / 2 + offset_y
  )

node_strength <- edges %>%
  group_by(Source) %>%
  summarise(InfluenceStrength = sum(abs(Influence))) %>%
  ungroup()

layout <- layout %>%
  left_join(node_strength, by = c("name" = "Source")) %>%
  mutate(
    InfluenceStrength = replace_na(InfluenceStrength, 0),
    Highlight = InfluenceStrength == max(InfluenceStrength)
  )

ggraph(layout) +
  geom_edge_fan(
    aes(width = abs(Influence), color = Influence),
    arrow = arrow(length = unit(3, 'mm')),
    end_cap = circle(4, 'mm'),
    lineend = "round"
  ) +
  geom_node_point(aes(size = InfluenceStrength, fill = Highlight),
                  color = "black", shape = 21, stroke = 1.2) +
  scale_fill_manual(values = c(`TRUE` = "orange", `FALSE` = "white")) +
  geom_node_text(aes(label = name, x = x * 1.15, y = y * 1.15),
                 size = 4, fontface = "bold") +
  geom_text(data = edges,
            aes(x = mid_x, y = mid_y, label = sprintf("%.2f", Influence)),
            size = 3, color = "black", fontface = "italic") +
  scale_edge_color_gradient2(low = "red", mid = "gray95", high = "green", midpoint = 0) +
  theme_void() +
  labs(
    title = "Whale Influence Network (g2)",
    fill = "Most Influential",
    size = "Total Influence Received"
  )
# save as svg to fix in inkscape
ggsave("dbw_infl_whale_network_g2.svg", width = 10, height = 10, units = "in", dpi = 300) 


#### ----------- NOT USING THESE Compare posterior and actual delays g2 ---------- ####


# Ensure events are ordered
g2_data <- g2_data %>%
  filter(Group == "g2") %>%
  arrange(bout, StartTime) %>%
  mutate(
    WhaleIndex = as.integer(factor(WhaleID)),
    BoutIndex = as.integer(factor(bout)),
    PreviousStartTime = lag(StartTime),
    ObservedDelay = StartTime - PreviousStartTime
  ) %>%
  filter(!is.na(ObservedDelay))  # remove first call in each bout


draws <- fit_g2$draws()
omega_draws <- draws %>% spread_draws(omega[i])
lambda_draws <- draws %>% spread_draws(lambda)
a_raw_draws <- draws %>% spread_draws(A_raw[k])


pair_lookup <- tibble(
  k = 1:length(stan_data_g2$pair_i),
  i = stan_data_g2$pair_i,
  j = stan_data_g2$pair_j
)

a_matrix <- a_raw_draws %>%
  left_join(pair_lookup, by = "k")


# Get posterior means
omega_mean <- omega_draws %>% group_by(i) %>% summarise(omega = mean(omega))
lambda_mean <- mean(lambda_draws$lambda)
a_mean <- a_matrix %>% group_by(i, j) %>% summarise(Aij = mean(A_raw))

# Join means into g2_data
g2_joined <- g2_data %>%
  left_join(omega_mean, by = c("WhaleIndex" = "i"))

# We'll build this loop-style (you can vectorize later)
mu_list <- vector("numeric", length = nrow(g2_joined))

for (n in seq_len(nrow(g2_joined))) {
  row_n <- g2_joined[n, ]
  i <- row_n$WhaleIndex
  t_n <- row_n$StartTime
  bout_n <- row_n$BoutIndex
  
  # All prior events in the same bout
  prior_calls <- g2_joined[1:(n - 1), ] %>%
    filter(BoutIndex == bout_n)
  
  # Sum decayed influences
  influence <- 0
  
  for (m in seq_len(nrow(prior_calls))) {
    j <- prior_calls$WhaleIndex[m]
    t_m <- prior_calls$StartTime[m]
    delta_t <- t_n - t_m
    
    Aij <- a_mean %>%
      filter(i == i, j == j) %>%
      pull(Aij)
    
    # But the above still uses `i` and `j` as loop variables — shadowing
    # Fix this by renaming inside the loop
    
    i_call <- row_n$WhaleIndex
    t_n <- row_n$StartTime
    bout_n <- row_n$BoutIndex
    
    for (m in seq_len(nrow(prior_calls))) {
      j_call <- prior_calls$WhaleIndex[m]
      t_m <- prior_calls$StartTime[m]
      delta_t <- t_n - t_m
      
      Aij <- a_mean %>%
        filter(i == i_call, j == j_call) %>%
        pull(Aij)
      
      influence <- influence + (Aij * exp(-lambda_mean * delta_t))
    }
  }
}
g2_joined$PredictedDelay <- mu_list


ggplot(g2_joined, aes(x = ObservedDelay, y = PredictedDelay)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Observed vs Predicted Call Delays (Posterior Mean)",
    x = "Observed Delay (sec)",
    y = "Predicted Delay (μₙ)"
  ) +
  theme_minimal()


#### -------- Influence exerced and recieced test plot, most of these are extracted all over th place already-- ####


# 1. Extract A[i,j] summary
A_summary_g2 <- fit_g2$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g2$pair_i[k], j = stan_data_g2$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# 2. Influence received (row sums)
influence_received_g2 <- A_summary_g2 %>%
  group_by(i) %>%
  summarise(total_received = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

# 3. Influence exerted (column sums)
influence_exerted_g2 <- A_summary_g2 %>%
  group_by(j) %>%
  summarise(total_exerted = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(j)])

# 4. ωᵢ summary
omega_summary_g2 <- fit_g2$draws("omega") %>%
  spread_draws(omega[i]) %>%
  group_by(i) %>%
  summarise(omega_median = median(omega)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

# 5. Combine all metrics
influence_profile_g2 <- full_join(influence_received_g2, influence_exerted_g2, by = "Whale") %>%
  full_join(omega_summary_g2, by = "Whale") %>%
  replace_na(list(total_received = 0, total_exerted = 0))

# 6. Plot: ωᵢ vs. influence exerted and received
ggplot(influence_profile_g2, aes(x = total_exerted, y = total_received, label = Whale)) +
  geom_point(aes(size = omega_median), color = "steelblue", alpha = 0.7) +
  geom_text_repel(size = 4) +
  scale_size_continuous(name = expression(omega[i]), range = c(3, 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  labs(
    title = "Whale Social Role: Influence vs. Baseline Calling",
    x = "Total Influence Exerted (Outgoing Aᵢⱼ)",
    y = "Total Influence Received (Incoming Aᵢⱼ)"
  ) +
  theme_minimal(base_size = 14)

ggsave("Influence_give_receive_g2.png", width = 10, height = 10, units = "in", dpi = 300)

#### --------------- Group 2 ---------------- ####

# ---- Posterior distribution of baseline calling rates ωᵢ
# Extract posterior draws for omega[i]
omega_plot_df_g2 <- fit_g2$draws("omega") %>%
  spread_draws(omega[i]) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])  # map index -> ID

# Plot posterior densities. They are bunched up at close to zero, lower bound for omega was set at 0.2
# That constraint is being hit by the posterior, meaning your data favor low intrinsic rates
# Suggests most of these whales are not initiating calls independently, they're socially reactive

ggplot(omega_plot_df_g2, aes(x = omega, fill = Whale)) +
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "red")+  # shows Stan’s lower bound
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Group 2 - Posterior of Baseline Calling Rates (ωᵢ)",
    x = expression(omega[i]),
    y = "Density",
    fill = "Whale"
  ) +
  theme_minimal(base_size = 14)

ggsave("omega_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# ---- Overlay influence and omega
omega_summary_g2 <- fit_g2$draws("omega") %>%
  spread_draws(omega[i]) %>%
  group_by(i) %>%
  summarise(omega_median = median(omega)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

A_summary_g2 <- fit_g2$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g2$pair_i[k], j = stan_data_g2$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Sum influence received by each whale (i.e., incoming influence to i)
incoming_influence_g2 <- A_summary_g2 %>%
  group_by(i) %>%
  summarise(total_influence_g2 = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g2[as.character(i)])

plot_overlay_g2 <- omega_summary_g2 %>%
  left_join(incoming_influence_g2, by = c("i", "Whale"))

# Added colour by whale, and icon by subgroup
plot_overlay_g2_col <- plot_overlay_g2 %>%
  left_join(g2_df %>% select(WhaleID, Source_subgroup) %>% distinct(), 
            by = c("Whale" = "WhaleID"))

ggplot(plot_overlay_g2_col, aes(x = total_influence_g2, y = omega_median)) +
  geom_point(aes(color = as.factor(Source_subgroup))) +
  geom_text_repel(aes(label = Whale), max.overlaps = 10, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Whales' Baseline Calling Rate vs. Social Influence (Group 2)",
    legend = element_text("Subgroup"),
    x = "Social influence: Sum of all incoming Aᵢⱼ values",
    y = expression("Median "*omega[i]*" (Baseline Calling Rate)"),
    color = "Subgroup"
  ) +
  theme_minimal(base_size = 14)

ggsave("omega_influce_overlay_g2.png", width = 10, height = 10, units = "in", dpi = 300)


# ---- Plot Posterior for lambda 

# Extract posterior draws for lambda
lambda_draws <- fit_g2$draws("lambda", format = "draws_df") %>%
  spread_draws(lambda)

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

ggsave("lambda_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# Diagnostics for lambda because it looks very strange, bimodal, with two peaks near 0.2 and 1.7, valley close to 1
fit_g2$summary("lambda")
fit_g2$cmdstan_diagnose()
posterior::rhat(fit_g2$draws("lambda"))

mcmc_trace(fit_g2$draws("lambda"))

### Diagnose G3
# Extract posterior draws for lambda
lambda_draws_g3 <- fit_g3$draws("lambda", format = "draws_df") %>%
  spread_draws(lambda)

# Plot posterior density
ggplot(lambda_draws_g3, aes(x = lambda)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "gray40", linewidth = 1) +
  labs(
    title = "G3 Posterior Distribution of λ (Decay Rate)",
    x = expression(lambda),
    y = "Posterior Density"
  ) +
  theme_minimal(base_size = 14)

ggsave("lambda_post_g2.png", width = 10, height = 10, units = "in", dpi = 300)

# Diagnostics for lambda because it looks very strange, bimodal, with two peaks near 0.2 and 1.7, valley close to 1
fit_g3$summary("lambda")
fit_g3$cmdstan_diagnose()
posterior::rhat(fit_g3$draws("lambda"))

mcmc_trace(fit_g3$draws("lambda"))

# ---- Plot influence accross specific bouts, longer than 5 calls

long_bouts_g2 <- g2_df %>%
  group_by(bout) %>%
  filter(n() > 5) %>%
  pull(bout) %>%
  unique()


# Step: summarize A[i,j]
A_matrix <- fit_g2$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g2$pair_i[k], j = stan_data_g2$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Convert to matrix for fast lookup
N <- max(stan_data_g2$caller)
A_mat <- matrix(0, N, N)
for (row in 1:nrow(A_matrix)) {
  A_mat[A_matrix$i[row], A_matrix$j[row]] <- A_matrix$mean_A[row]
}

lambda_mean <- summarise_draws(fit_g2$draws("lambda"))$mean

lambda <- lambda_mean

influence_traces <- list()

for (b in long_bouts_g2[1:10]) {
  sub_bout <- g2_df %>%
    filter(bout == b) %>%
    arrange(StartTime) %>%
    mutate(event = row_number())
  
  n_events <- nrow(sub_bout)
  
  traces <- tibble()
  
  for (n in 2:n_events) {
    i <- sub_bout$WhaleIndex[n]
    t_n <- sub_bout$StartTime[n]
    
    influence_sum <- 0
    
    for (m in 1:(n - 1)) {
      j <- sub_bout$WhaleIndex[m]
      t_m <- sub_bout$StartTime[m]
      delta_t <- t_n - t_m
      decay <- exp(-lambda * delta_t)
      
      influence_sum <- influence_sum + A_mat[i, j] * decay
    }
    
    traces <- bind_rows(traces, tibble(
      Bout = b,
      Event = n,
      Time = t_n,
      Whale = sub_bout$WhaleID[n],
      Influence = influence_sum
    ))
  }
  
  influence_traces[[as.character(b)]] <- traces
}

influence_df <- bind_rows(influence_traces)
call_lines_df <- g2_df %>%
  filter(bout %in% long_bouts_g2[1:10]) %>%
  select(Bout = bout, Time = StartTime)


ggplot(influence_df, aes(x = Time, y = Influence, color = Whale)) +
  #geom_step(direction = "hv") +
  geom_point(aes(color = Whale, size = abs(Influence))) +
  geom_vline(data = call_lines_df, aes(xintercept = Time), color = "black", linetype = "dotted", size = 0.4) +
  
  facet_wrap(~ Bout, scales = "free_x") +
  labs(
    title = "Cumulative Influence Across Longer Bouts",
    x = "Time within Bout",
    y = "Summed Influence from Previous Calls"
  ) +
  theme_minimal()

# ---- Network graph


# Create a tidy edge list from A[i,j]
edge_df_g2 <- A_summary_g2 %>%
  mutate(
    from_whale = whale_index_to_id_g2[as.character(j)],
    to_whale   = whale_index_to_id_g2[as.character(i)]
  ) %>%
  select(from = from_whale, to = to_whale, mean_A)
# Optional threshold

# Create graph object
g2_influence <- graph_from_data_frame(edge_df_g2, directed = TRUE)


# Compute strength (weighted degree)
V(g2_influence)$in_strength  <- strength(g2_influence, mode = "in", weights = E(g2_influence)$mean_A)
V(g2_influence)$out_strength <- strength(g2_influence, mode = "out", weights = E(g2_influence)$mean_A)

# Centrality
V(g2_influence)$eigen_centrality <- eigen_centrality(g2_influence, weights = E(g2_influence)$mean_A)$vector
V(g2_influence)$betweenness <- betweenness(g2_influence, directed = TRUE)


ggraph(g2_influence, layout = "nicely") +
  geom_edge_fan(aes(width = abs(mean_A), color = mean_A),
                arrow = arrow(length = unit(3, "mm")), end_cap = circle(3, "mm")) +
  geom_node_point(aes(size = in_strength)) +
  geom_node_label(aes(label = name)) +
  scale_edge_color_gradient2(low = "red", mid = "gray80", high = "green", midpoint = 0) +
  labs(title = "Whale Influence Network (from Aᵢⱼ)",
       subtitle = "Edge color = influence direction, width = strength, node size = total received influence") +
  theme_void()

### -- Add node strength to graph

# Step 1: Format edge list for plotting
edges <- A_summary_g2 %>%
  mutate(
    Source = whale_index_to_id_g2[as.character(j)],  # j = influencer
    Target = whale_index_to_id_g2[as.character(i)],  # i = influenced
    Influence = mean_A
  ) %>%
  filter(abs(Influence) > 0.01) %>% # Filter close to zero influence
  select(Source, Target, Influence)

# Graph and layout
g <- graph_from_data_frame(edges, directed = TRUE)

layout <- create_layout(g, layout = 'linear', circular = TRUE)

edges <- edges %>%
  mutate(
    from_x = layout$x[match(Source, layout$name)],
    from_y = layout$y[match(Source, layout$name)],
    to_x   = layout$x[match(Target, layout$name)],
    to_y   = layout$y[match(Target, layout$name)],
    dx = to_x - from_x,
    dy = to_y - from_y,
    norm = sqrt(dx^2 + dy^2),
    offset_x = -dy / norm * 0.05,
    offset_y = dx / norm * 0.05,
    mid_x = (from_x + to_x) / 2 + offset_x,
    mid_y = (from_y + to_y) / 2 + offset_y
  )

node_strength <- edges %>%
  group_by(Source) %>%
  summarise(InfluenceStrength = sum(abs(Influence))) %>%
  ungroup()

layout <- layout %>%
  left_join(node_strength, by = c("name" = "Source")) %>%
  mutate(
    InfluenceStrength = replace_na(InfluenceStrength, 0),
    Highlight = InfluenceStrength == max(InfluenceStrength)
  )

ggraph(layout) +
  geom_edge_fan(
    aes(width = abs(Influence), color = Influence),
    arrow = arrow(length = unit(3, 'mm')),
    end_cap = circle(4, 'mm'),
    lineend = "round"
  ) +
  geom_node_point(aes(size = InfluenceStrength, fill = Highlight),
                  color = "black", shape = 21, stroke = 1.2) +
  scale_fill_manual(values = c(`TRUE` = "orange", `FALSE` = "white")) +
  geom_node_text(aes(label = name, x = x * 1.15, y = y * 1.15),
                 size = 4, fontface = "bold") +
  geom_text(data = edges,
            aes(x = mid_x, y = mid_y, label = sprintf("%.2f", Influence)),
            size = 3, color = "black", fontface = "italic") +
  scale_edge_color_gradient2(low = "red", mid = "gray95", high = "green", midpoint = 0) +
  theme_void() +
  labs(
    title = "Whale Influence Network (g2)",
    fill = "Most Influential",
    size = "Total Influence Received"
  )
# save as svg to fix in inkscape
ggsave("dbw_infl_whale_network_g2.svg", width = 10, height = 10, units = "in", dpi = 300) 















# Group 3
# Extract posterior draws for omega[i]
omega_plot_df_g3 <- fit_g3$draws("omega") %>%
  spread_draws(omega[i]) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])  # map index → ID

# Plot posterior densities
ggplot(omega_plot_df_g3, aes(x = omega, fill = Whale)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Group 3 - Posterior of Baseline Calling Rates (ωᵢ)",
    x = expression(omega[i]),
    y = "Density",
    fill = "Whale"
  ) +
  theme_minimal(base_size = 14)


#### --------------- Influence Matrix — Aij ----------------- ####

# --- Group 1
pair_i <- stan_data_g2$pair_i
pair_j <- stan_data_g2$pair_j


fit_g2$draws() %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = pair_i[k], j = pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop") %>%
  left_join(whale_id_map_g2, by = c("i" = "WhaleIndex")) %>%
  rename(From = WhaleID) %>%
  left_join(whale_id_map_g2, by = c("j" = "WhaleIndex")) %>%
  rename(To = WhaleID) %>%
  ggplot(aes(x = From, y = To, fill = mean_A)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(title = "Group 1: Influence Matrix A[i,j]: Who Drives Whom?",
       x = "Influenced Whale (i)", y = "Influencing Whale (j)",
       fill = "Mean Aᵢⱼ") +
  theme_minimal()

# --- Group 2 
pair_i <- stan_data_g2$pair_i
pair_j <- stan_data_g2$pair_j


fit_g2$draws() %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = pair_i[k], j = pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop") %>%
  left_join(whale_id_map_g2, by = c("i" = "WhaleIndex")) %>%
  rename(From = WhaleID) %>%
  left_join(whale_id_map_g2, by = c("j" = "WhaleIndex")) %>%
  rename(To = WhaleID) %>%
  ggplot(aes(x = From, y = To, fill = mean_A)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(title = "Group 2: Influence Matrix A[i,j]: Who Drives Whom?",
       x = "Influenced Whale (i)", y = "Influencing Whale (j)",
       fill = "Mean Aᵢⱼ") +
  theme_minimal()


# --- Group 3 
pair_i <- stan_data_g3$pair_i
pair_j <- stan_data_g3$pair_j


fit_g3$draws() %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = pair_i[k], j = pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop") %>%
  left_join(whale_id_map_g3, by = c("i" = "WhaleIndex")) %>%
  rename(From = WhaleID) %>%
  left_join(whale_id_map_g3, by = c("j" = "WhaleIndex")) %>%
  rename(To = WhaleID) %>%
  ggplot(aes(x = From, y = To, fill = mean_A)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(title = "Group 3: Influence Matrix A[i,j]: Who Drives Whom?",
       x = "Influenced Whale (i)", y = "Influencing Whale (j)",
       fill = "Mean Aᵢⱼ") +
  theme_minimal()


#### --------------- Influence Network Graph ----------------- ####
levels <- levels(factor(g2_df$WhaleID))
data.frame(WhaleIndex = seq_along(levels), WhaleID = levels)


# Assume: pair_i and pair_j are already defined from stan_data_g2
# Pull WhaleID labels into named vectors (not join twice)
from_labels <- whale_id_map_g2 %>% rename(from = WhaleID)
to_labels   <- whale_id_map_g2 %>% rename(to = WhaleID)

# Build edge_df with minimal joins to avoid name conflicts
edge_df <- A_df %>%
  left_join(from_labels, by = c("i" = "WhaleIndex")) %>%
  left_join(to_labels, by = c("j" = "WhaleIndex")) %>%
  select(from, to, mean_A)  # keep only what you need

# Create the igraph object
g <- graph_from_data_frame(edge_df, directed = TRUE)

# Now plot with ggraph
ggraph(g, layout = "circle") +
  geom_edge_link(aes(width = abs(mean_A), color = mean_A),
                 arrow = arrow(length = unit(3, "mm")), end_cap = circle(3, "mm")) +
  geom_node_label(aes(label = name)) +
  scale_edge_color_gradient2(low = "blue", mid = "gray80", high = "red", midpoint = 0) +
  theme_void() +
  labs(title = "Group 1: Directed Influence Network (Aᵢⱼ)")


#### --------------- Posterior Predictive Checks ----------------- ####