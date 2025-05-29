# Analysis and visualisation of discrete bout-wide (cumulative) influence model: g2, g2, G3

# Add all the things
pacman::p_load(tidyverse, cmdstanr, tidybayes, dplyr, ggplot2, igraph, ggraph, bayesplot, ggrepel)  
               


#### --------------- Load all the needed files ---------------- ####
# Load fitted models
fit_g1 <- readRDS("../fits/cumulative_influence/fit_g1_lambda_normprior_all.rds")
fit_g2 <- readRDS("../fits/cumulative_influence/fit_g2_lambda_normprior_all.rds")
fit_g3 <- readRDS("../fits/cumulative_influence/fit_g3_lambda_normprior_all.rds")

# Load stan data list
stan_data_g1 <- readRDS("saved_rds/g1_stan_data.rds")
stan_data_g2 <- readRDS("saved_rds/g2_stan_data.rds")
stan_data_g3 <- readRDS("saved_rds/g3_stan_data.rds")

# Load subsetted data
g1_data <- readRDS("saved_rds/g1_data.rds")
g2_data <- readRDS("saved_rds/g2_data.rds")
g3_data<- readRDS("saved_rds/g3_data.rds")

# Id map
whale_id_map_g1 <- readRDS("saved_rds/g1_whale_id_map.rds")
whale_id_map_g2 <- readRDS("saved_rds/g2_whale_id_map.rds")
whale_id_map_g3 <- readRDS("saved_rds/g3_whale_id_map.rds")

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g1 <- setNames(whale_id_map_g1$WhaleID, whale_id_map_g2$WhaleIndex)

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g2 <- setNames(whale_id_map_g2$WhaleID, whale_id_map_g2$WhaleIndex)

# Create vector: names are indices, values are WhaleIDs
whale_index_to_id_g3 <- setNames(whale_id_map_g3$WhaleID, whale_id_map_g3$WhaleIndex)






#### --------------- Plots, just find and replace g(x) for different groups ---------------- ####

# ---- Posterior distribution of baseline calling rates ωᵢ
# Extract posterior draws for omega[i]
omega_plot_df_g3 <- fit_g3$draws("omega") %>%
  spread_draws(omega[i]) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])  # map index -> ID

# Plot posterior densities. They are bunched up at close to zero, lower bound for omega was set at 0.2
# That constraint is being hit by the posterior, meaning your data favor low intrinsic rates
# Suggests most of these whales are not initiating calls independently, they're socially reactive

ggplot(omega_plot_df_g3, aes(x = omega, fill = Whale)) +
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

ggsave("omega_post_g3_lambdamodel.png", width = 10, height = 10, units = "in", dpi = 300)

# ---- Overlay influence and omega
omega_summary_g3 <- fit_g3$draws("omega") %>%
  spread_draws(omega[i]) %>%
  group_by(i) %>%
  summarise(omega_median = median(omega)) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])

A_summary_g3 <- fit_g3$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g3$pair_i[k], j = stan_data_g3$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Sum influence received by each whale (i.e. incoming influence to i)
incoming_influence_g3 <- A_summary_g3 %>%
  group_by(i) %>%
  summarise(total_influence_g3 = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])

plot_overlay_g3 <- omega_summary_g3 %>%
  left_join(incoming_influence_g3, by = c("i", "Whale"))

# Added colour by whale, and icon by subgroup
plot_overlay_g3_col <- plot_overlay_g3 %>%
  left_join(g3_data %>% select(WhaleID, Source_subgroup) %>% distinct(), 
            by = c("Whale" = "WhaleID"))

ggplot(plot_overlay_g3_col, aes(x = total_influence_g3, y = omega_median)) +
  geom_point(aes(color = as.factor(Source_subgroup))) +
  geom_text_repel(aes(label = Whale), max.overlaps = 10, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Whales' Baseline Calling Rate vs. Social Influence (Group g3)",
    legend = element_text("Subgroup"),
    x = "Social influence: Sum of all incoming Aᵢⱼ values",
    y = expression("Median "*omega[i]*" (Baseline Calling Rate)"),
    color = "Subgroup"
  ) +
  theme_minimal(base_size = 14)

#ggsave("omega_influce_overlay_g3_lambdamodel.png", width = 10, height = 10, units = "in", dpi = 300)


# ---- Plot Posterior for lambda. Model is wonky, lambda looks weird, we already knew this from fit diagnostics! proceed with care

# Extract posterior draws for lambda
lambda_draws <- fit_g3$draws("lambda", format = "draws_df") %>%
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

#ggsave("lambda[i]_post_g3.png", width = 10, height = 10, units = "in", dpi = 300)


#### --------------- Plot influence accross specific bouts, longer than 5 calls - So hard to read sadly ---------------- ####

# Identify longer bouts (more than 5 calls) 
# When not using bouts, all long_bouts_g2[1:10] are replaced by g2_data$bout  
# long_bouts_g2 <- g2_data %>%
#   group_by(bout) %>%
#   filter(n() > 5) %>%
#   pull(bout) %>%
#   unique()


# ---- Summarize A[i,j] = mean across posterior draws
A_matrix <- fit_g3$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g3$pair_i[k], j = stan_data_g3$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Convert A[i,j] to matrix for fast lookup
N <- max(stan_data_g3$caller)
A_mat <- matrix(0, N, N)
for (row in 1:nrow(A_matrix)) {
  A_mat[A_matrix$i[row], A_matrix$j[row]] <- A_matrix$mean_A[row]
}

# ---- Extract per-whale lambda[i] posterior means 
lambda_df <- fit_g3$draws("lambda") %>%
  spread_draws(lambda[i]) %>%
  group_by(i) %>%
  summarise(lambda_i = mean(lambda))

# Named vector: whale index -> lambda
lambda_vec <- setNames(lambda_df$lambda_i, as.character(lambda_df$i))

# ---- Trace influence over time in each bout: cumulative influence a whale receives during each call it makes

# Within a bout: Calculating the cumulative influence [i] recieves during each call: 
# Based on previous whales who called before them:
# The decay effect of earlier calls (time delay and pr-whale decay rate lambda[i])
# The strength of influence connection from the influence matrix A_mat[i,j]

# to store influence traces pr bout
influence_traces <- list()

# loop opver each bout unique bout number
for (b in unique(g3_data$bout)) {
  # filter data to current bout, sort chronologically by Starttime, add "event" coloumn to index call events (nth call in time order)
  sub_bout <- g3_data %>%
    filter(bout == b) %>%
    arrange(StartTime) %>%
    mutate(event = row_number())
  
  # Get number of calll events in the bout, and create a tibble to store influence values for each event in this bout
  n_events <- nrow(sub_bout)
  traces <- tibble()
  
  # Start at second call (can't be influenced by nothing)
  for (n in 2:n_events) {
    # index for whale at call n
    i <- sub_bout$WhaleIndex[n]
    # t_n the time of that call (event)
    t_n <- sub_bout$StartTime[n]
    # add sum or accumumulated influence on this whale from earlier calls
    influence_sum <- 0
    
    # loop through all earlier calls m (1 to n-1) in bout
    for (m in 1:(n - 1)) {
      # Ahle who made earlier call
      j <- sub_bout$WhaleIndex[m]
      # bout start at call of whale m, so previous call to n
      t_m <- sub_bout$StartTime[m]
      # time difference between earlier call and current one
      delta_t <- t_n - t_m
      # decay rate for receiver whale i
      lambda_i <- lambda_vec[as.character(i)]
      # how much the influence from the earlier call has decayed by the time of i's call
      decay <- exp(-lambda_i * delta_t)
      # A_mat[i,j the strength of influence from whale j to whale i. The product of that times decayed influence of that earlier call on this
      # whale, add to total influence sum
      influence_sum <- influence_sum + A_mat[i, j] * decay
    }
    
    # Store all results: bout nr, event nr, callt ime, ID, total influence from previous calls
    traces <- bind_rows(traces, tibble(
      Bout = b,
      Event = n,
      Time = t_n,
      Whale = sub_bout$WhaleID[n],
      Influence = influence_sum
    ))
  }
  
  # Save in list, using bout number as "key"
  influence_traces[[as.character(b)]] <- traces
}

# ---- Combine results and annotate call times ----
influence_df <- bind_rows(influence_traces)

# Before i filter for interesting bouts - get an overview and pick a few out after
call_lines_df <- g3_data %>%
  select(Bout = bout, Time = StartTime)

# Filter only the bouts of interest. A fair amount didn't have much happening, so filter only the bouts with lots of activity
interesting_bouts <- c( 9, 12, 21, 37, 38)

influence_df <- influence_df %>% filter(Bout %in% interesting_bouts)
call_lines_df <- call_lines_df %>% filter(Bout %in% interesting_bouts)


# ---- Plot ----
ggplot(influence_df, aes(x = Time, y = Influence, color = Whale)) +
  geom_point(aes(size = abs(Influence))) +
  geom_vline(data = call_lines_df, aes(xintercept = Time),
             color = "grey", linetype = "dotted", size = 0.4) +
  facet_wrap(~ Bout, scales = "free_x") +
  labs(
    title = "Group 3: Influence Across Bouts 9, 12, 21, 37, 38",
    x = "Time within Bout",
    y = "Summed Influence from Previous Calls"
  ) +
  theme_minimal()



#### ------------ Final Network plot, have to edit node strength text in inkscape, Edge_fan messes up alignment ---- ####
# Step 1: Format edge list for plotting
edges <- A_summary_g3 %>%
  mutate(
    Source = whale_index_to_id_g3[as.character(j)],  # j = influencer
    Target = whale_index_to_id_g3[as.character(i)],  # i = influenced
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
    title = "Whale Influence Network (Group 3)",
    fill = "Most Influential",
    size = "Total Influence Received"
  )
# save as svg to fix in inkscape
ggsave("dbwlambda_model_infl_whale_network_g3.svg", width = 10, height = 10, units = "in", dpi = 300) 


#### -------- Influence exerced and recieced test plot, most of these are extracted all over th place already-- ####

# Reextrcating everything, because it's probably a mess at this point. I should learn to make good functions

# Extract A[i,j] summary
A_summary_g3 <- fit_g3$draws("A_raw") %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = stan_data_g3$pair_i[k], j = stan_data_g3$pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop")

# Influence received (row sums)
influence_received_g3 <- A_summary_g3 %>%
  group_by(i) %>%
  summarise(total_received = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])

# Influence exerted (column sums)
influence_exerted_g3 <- A_summary_g3 %>%
  group_by(j) %>%
  summarise(total_exerted = sum(mean_A)) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(j)])

# omega summary
omega_summary_g3 <- fit_g3$draws("omega") %>%
  spread_draws(omega[i]) %>%
  group_by(i) %>%
  summarise(omega_median = median(omega)) %>%
  mutate(Whale = whale_index_to_id_g3[as.character(i)])

# Combine all metrics
influence_profile_g3 <- full_join(influence_received_g3, influence_exerted_g3, by = "Whale") %>%
  full_join(omega_summary_g3, by = "Whale") %>%
  replace_na(list(total_received = 0, total_exerted = 0))

# Plot: omega vs influence exerted and received
ggplot(influence_profile_g3, aes(x = total_exerted, y = total_received, label = Whale)) +
  geom_point(aes(size = omega_median), color = "steelblue", alpha = 0.7) +
  geom_text_repel(size = 4) +
  scale_size_continuous(name = expression(omega[i]), range = c(3, 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  labs(
    title = "Group 3: Whale Social Role: Influence vs. Baseline Calling",
    x = "Total Influence Exerted (Outgoing Aᵢⱼ)",
    y = "Total Influence Received (Incoming Aᵢⱼ)"
  ) +
  theme_minimal(base_size = 14)

ggsave("Lambda_model_Influence_give_receive_g3.png", width = 10, height = 10, units = "in", dpi = 300)



#### --------------- Influence Matrix — Aij - Heatmaps ----------------- ####

# --- Create pairs A[i,j]
pair_i <- stan_data_g1$pair_i
pair_j <- stan_data_g1$pair_j


fit_g1$draws() %>%
  spread_draws(A_raw[k]) %>%
  mutate(i = pair_i[k], j = pair_j[k]) %>%
  group_by(i, j) %>%
  summarise(mean_A = mean(A_raw), .groups = "drop") %>%
  left_join(whale_id_map_g1, by = c("i" = "WhaleIndex")) %>%
  rename(From = WhaleID) %>%
  left_join(whale_id_map_g1, by = c("j" = "WhaleIndex")) %>%
  rename(To = WhaleID) %>%
  ggplot(aes(x = From, y = To, fill = mean_A)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(title = "Group 1: Influence Matrix A[i,j]: Who Drives Whom?",
       x = "Influenced Whale (i)", y = "Influencing Whale (j)",
       fill = "Mean Aᵢⱼ") +
  theme_minimal()

#ggsave("Lambda_model_heatmap_g1.png", width = 10, height = 10, units = "in", dpi = 300)





