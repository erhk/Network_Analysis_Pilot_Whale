# Analysis and visualisation of Cumulative influence model G1, G2, G3

# Load fitted models
fit_g1 <- readRDS("fit_G1_dpi_cum.rds")
fit_g2 <- readRDS("fit_G2_dpi_cum.rds")
fit_g3 <- readRDS("fit_G3_dpi_cum.rds")

# Load stan data list
stan_data_g1 <- readRDS("g1_stan_data.rds")
stan_data_g2 <- readRDS("g2_stan_data.rds")
stan_data_g3 <- readRDS("g3_stan_data.rds")

# Load subsetted data
g1_df <- readRDS("g1_data.rds")
g2_df <- readRDS("g2_data.rds")
g3_df<- readRDS("g3_data.rds")

# Id map
whale_id_map_g1 <- readRDS("g1_whale_id_map.rds")
whale_id_map_g2 <- readRDS("g2_whale_id_map.rds")
whale_id_map_g3 <- readRDS("g3_whale_id_map.rds")
