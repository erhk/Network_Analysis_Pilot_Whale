### Helper functions 

pacman::p_load(ggplot2, posterior, reshape2, dplyr)


# Load data
prep <- read.delim("../../Data/PrepData.csv", sep = ",")

# # # -------- Prep data for stan -------- # # #
prepare_discrete_phase_with_influence <- function(group_name, prep_data, time_scale = 1) {
  # Step 1: define unique WhaleIDs for this group only
  unique_ids <- prep_data %>%
    filter(Group == group_name) %>%
    distinct(WhaleID) %>%
    arrange(WhaleID) %>%
    pull(WhaleID)
  
  # Step 2: generate group_data with clean WhaleIndex
  group_data <- prep_data %>%
    filter(Group == group_name) %>%
    arrange(bout, StartTime) %>%
    group_by(bout) %>%
    arrange(StartTime, .by_group = TRUE) %>%
    mutate(
      CallTime   = StartTime * time_scale,
      WhaleIndex = as.integer(factor(WhaleID, levels = unique_ids)),
      bout_id    = as.integer(factor(bout))
    ) %>%
    ungroup()
  
  # Step 3: rebuild consistent unique_whales AFTER cleaning WhaleIndex
  unique_whales <- sort(unique(group_data$WhaleIndex))
  
  # Step 4: build pairwise data
  pair_data <- expand.grid(i = unique_whales, j = unique_whales) %>%
    filter(i != j)
  
  # Step 5: build Stan list
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





# # # -------- Summarize omega -------- # # #
summarize_omega <- function(fit, data_group) {
  omega_summary <- fit$summary(variables = "omega") %>%
    mutate(whale = as.integer(sub("omega\\[", "", sub("\\]", "", variable))))
  
  index_map <- data_group$group_data %>%
    distinct(WhaleIndex, WhaleID) %>%
    mutate(WhaleID = as.character(WhaleID))
  
  omega_summary_named <- omega_summary %>%
    left_join(index_map, by = c("whale" = "WhaleIndex")) %>%
    mutate(
      status = ifelse(is.na(WhaleID), "inactive", "active"),
      WhaleID = ifelse(is.na(WhaleID), paste0("inactive_", whale), WhaleID)
    )
  
  return(omega_summary_named)
}

# -------- continue --------

