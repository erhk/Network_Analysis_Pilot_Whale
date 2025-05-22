# Whale Influence Network Visualization (G1)

# ---------------------- Libraries
pacman::p_load(tidyr, cmdstanr, dplyr, reshape2, ggplot2, igraph, tidygraph, scales, ggraph)

# ---------------------- Load Data
fit <- readRDS("fit_kuramoto_G1_small.rds")
load("data_G1_small.RData")  # contains data_G1_small

# ---------------------- Influence Matrix A[i,j] 
posterior <- fit$draws(variables = "A_raw", format = "draws_matrix")
A_mean <- colMeans(posterior)

A_mat <- matrix(0, data_G1_small$data$N_whales, data_G1_small$data$N_whales)
for (k in seq_along(A_mean)) {
  i <- data_G1_small$data$pair_i[k]
  j <- data_G1_small$data$pair_j[k]
  A_mat[i, j] <- A_mean[k]
}

whale_lookup <- data_G1_small$group_data %>%
  select(WhaleID, WhaleIndex) %>%
  distinct() %>%
  arrange(WhaleIndex)

rownames(A_mat) <- whale_lookup$WhaleID
colnames(A_mat) <- whale_lookup$WhaleID

# ----------------------  Heatmap
A_long <- melt(A_mat, varnames = c("receiver", "source"))
heatmap_plot <- ggplot(A_long, aes(x = source, y = receiver, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Estimated Pairwise Influence Matrix A",
       x = "Source Whale (j)", y = "Receiver Whale (i)", fill = "A[i,j]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("whale_A_matrix_heatmap.png", heatmap_plot, width = 7, height = 6)

# ----------------------  Graph object and centrality 
g <- graph_from_adjacency_matrix(A_mat, mode = "directed", weighted = TRUE)
V(g)$name <- rownames(A_mat)
V(g)$in_strength <- rowSums(A_mat)
V(g)$out_strength <- colSums(A_mat)
V(g)$degree <- degree(g, mode = "all")
V(g)$eigen <- eigen_centrality(g, directed = TRUE, weights = abs(E(g)$weight))$vector
V(g)$betweenness <- betweenness(g, directed = TRUE, weights = abs(E(g)$weight))

if ("Source_subgroup" %in% colnames(data_G1_small$group_data)) {
  whale_meta <- data_G1_small$group_data %>%
    select(WhaleID, Source_subgroup, WhaleIndex) %>%
    distinct() %>%
    arrange(WhaleIndex)
  V(g)$subgroup <- whale_meta$Source_subgroup
}

# ---------------------- GGraph plot, cannot deal with negative values (weights)
tg <- as_tbl_graph(g) %>%
  mutate(
    total_strength = in_strength + out_strength,
    size = rescale(total_strength, to = c(4, 12))
  )

E(tg)$sign <- ifelse(E(tg)$weight > 0, "positive", "negative")
E(tg)$weight_abs <- abs(E(tg)$weight)

ggraph(tg, layout = "circle") +
  geom_edge_link(aes(color = sign, width = weight_abs),
                 arrow = arrow(length = unit(4, "mm")),
                 end_cap = circle(3, 'mm')) +
  scale_edge_color_manual(values = c("positive" = "#1f78b4", "negative" = "#e31a1c")) +
  scale_edge_width(range = c(0.2, 2.5)) +
  geom_node_point(aes(size = size), shape = 21, fill = "white", color = "black") +
  geom_node_text(aes(label = name), size = 4, repel = FALSE) +
  theme_graph() +
  labs(title = "Whale Influence Network (G1)",
       edge_color = "Influence Type", edge_width = "Influence Strength")

# ----------------------  Exporting: Matrix, Nodes_df, Influence_graph mm.
write.csv(A_mat, "whale_influence_adjacency_matrix.csv", row.names = TRUE)

node_df <- data.frame(
  WhaleID = V(g)$name,
  InStrength = V(g)$in_strength,
  OutStrength = V(g)$out_strength,
  Eigenvector = V(g)$eigen,
  Betweenness = V(g)$betweenness
)

write.csv(node_df, "whale_node_centrality.csv", row.names = FALSE)
saveRDS(g, "whale_influence_graph_full.rds")
