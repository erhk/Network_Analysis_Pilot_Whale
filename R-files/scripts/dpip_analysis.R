# # # ------- Analysis and Visualisation ------- # # #

pacman::p_load(ggplot2, posterior, reshape2, dplyr, igraph, ggraph, tidyr, tibble, ggrepel, svglite)
#source("prep_function.r")


# Melt the A matrix
edges <- as.data.frame(A_matrix) %>%
  rownames_to_column(var = "Target") %>%
  pivot_longer(-Target, names_to = "Source", values_to = "Influence") %>%
  filter(Target != Source)  # remove self-loops


# Load data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")

# G1, G2, G3 data objects
data_G1_full <- prepare_discrete_phase_with_influence("G1", prep, time_scale = 5)
data_G2_full <- prepare_discrete_phase_with_influence("G2", prep, time_scale = 5)
data_G3_full <- prepare_discrete_phase_with_influence("G3", prep, time_scale = 5)


# Fitted models
fit_G1 <- readRDS("fit_G1_pairwise_fix.rds")
fit_G2 <- readRDS("fit_G2_pairwise.rds")
fit_G3 <- readRDS("fit_G3_pairwise.rds")

# # # ------- A_raw matrix ------- # # #

# Extract posterior means for A_raw
a_means <- colMeans(as_draws_df(fit_G2$draws("A_raw")))

# Reconstruct A matrix
N <- data_G2_full$data$N_whales
A_matrix <- matrix(0, N, N)
for (k in seq_along(a_means)) {
  i <- data_G2_full$data$pair_i[k]
  j <- data_G2_full$data$pair_j[k]
  A_matrix[i, j] <- a_means[k]
}



# # # ---------------------- Map whale ID to omega and A-matrix for clarity --- # # # 
# Map whaleID to the index number
id_map_G2 <- data_G2_full$group_data %>%
  select(WhaleIndex, WhaleID) %>%
  distinct() %>%
  arrange(WhaleIndex)


# omega summary
omega_summary_named_G2 <- summarize_omega(fit_G2, data_G2_full) %>%
  distinct(whale, .keep_all = TRUE) %>%
  left_join(id_map_G2, by = c("whale" = "WhaleIndex"))

# Drop extra id coloumn if i accidentially duplicate in above
omega_summary_named_G2 <- omega_summary_named_G2 %>%
  select(-WhaleID.x) %>%
  rename(WhaleID = WhaleID.y)

# Create labels for ID
labels <- omega_summary_named_G2 %>%
  arrange(whale) %>%        # ensure 1-to-N order
  pull(WhaleID)

# Label sources in A-matrix
rownames(A_matrix) <- labels  # Influenced whale (target)
colnames(A_matrix) <- labels  # Calling whale (source)



# # # ------------ Plot Heatmap  --- # # # 


A_df <- melt(A_matrix, varnames = c("Target", "Source"), value.name = "Influence")

ggplot(A_df, aes(x = Source, y = Target, fill = Influence)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "green", mid = "white", midpoint = 0) +
  labs(
    title = "G2 Whale Influence Matrix (A[i,j])",
    x = "Calling Whale",
    y = "Influenced Whale"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Add edge labels

edges <- as.data.frame(A_matrix) %>%
  rownames_to_column("Target") %>%
  pivot_longer(-Target, names_to = "Source", values_to = "Influence") %>%
  filter(Target != Source) %>%
  mutate(
    Strength = abs(Influence),
    Color = case_when(
      Influence > 0 ~ "green",
      Influence < 0 ~ "red",
      TRUE ~ "gray80"
    ),
    Label = sprintf("%.2f", Influence)  # Format to 2 decimal places
  )



#### ------------------ Ugly graph -----------------------

# graph objects
g <- graph_from_data_frame(edges, directed = TRUE)

# Edge midpoint coordinates for labels - made them prettier in inkscape
edges <- edges %>%
  mutate(
    from_x = layout$x[match(Source, layout$name)],
    from_y = layout$y[match(Source, layout$name)],
    to_x   = layout$x[match(Target, layout$name)],
    to_y   = layout$y[match(Target, layout$name)],
    dx = to_x - from_x,
    dy = to_y - from_y,
    norm = sqrt(dx^2 + dy^2),
    offset_x = -dy / norm * 0.05,  # Perpendicular offset
    offset_y = dx / norm * 0.05,
    mid_x = (from_x + to_x) / 2 + offset_x,
    mid_y = (from_y + to_y) / 2 + offset_y
  )

# Compute node strength
node_strength <- edges %>%
  group_by(Source) %>%
  summarise(InfluenceStrength = sum(abs(Influence))) %>%
  ungroup()

# Plot
layout <- create_layout(graph_from_data_frame(edges), layout = 'linear', circular = TRUE) %>%
  left_join(node_strength, by = c("name" = "Source")) %>%
  mutate(
    InfluenceStrength = replace_na(InfluenceStrength, 0),
    Highlight = if_else(InfluenceStrength == max(InfluenceStrength), TRUE, FALSE)
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
  geom_node_text(aes(label = name, x = x * 1.15, y = y * 1.15), size = 4, fontface = "bold") +
  geom_text(data = edges,
            aes(x = mid_x, y = mid_y, label = sprintf("%.2f", Influence)),
            size = 3, color = "black", fontface = "italic") +
  scale_edge_color_gradient2(low = "red", mid = "gray95", high = "green", midpoint = 0) +
  theme_void() +
  labs(title = "Whale Influence Network (G2)",
       fill = "Most Influential",
       size = "Total Influence")
  
#ggsave("whale_network_G2.svg", width = 10, height = 10, units = "in", dpi = 300)

