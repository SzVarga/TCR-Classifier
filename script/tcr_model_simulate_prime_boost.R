#!/usr/bin/env Rscript

# TCR Prime-Boost Model Simulation Example
# This script demonstrates the usage of the new prime-boost model
# with viral dynamics and T cell-virus interactions

# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# Load required libraries
library(adaptivetau)
library(sn)

# Source the TCR model functions
source("R/tcr_model.R")
source("R/tcr_model_prime_boost.R")

# Set simulation parameters
sim_times <- c("P10" = 10, "V2" = 10, "S10" = 10, "S68" = 58,
               "S210" = 152, "V3" = 20, "T10" = 10, "T108" = 98, "T189" = 79)
param_scale <- 100
num_generations <- 1
num_clones <- 3
clone_size <- 100
carry_cap <- num_clones * clone_size / 0.6
data_dir <- "data/tcr_data/"
data_name <- paste0("tcrColl_", num_clones, "clo_a", clone_size, "_x",
                    param_scale, ".rds")

# Example 1: Prime-boost vaccination schedule
cat("Example 1: Prime-boost vaccination schedule\n")
cat("=========================================\n")

# Define viral parameters for prime-boost
viral_params_prime_boost <- list(
  initial_load = 100,                           # Initial viral load
  replication_rate = param_scale * 0.008,       # High replication during active periods
  clearance_rate = param_scale * 0.001,         # Clearance rate
  replication_intervals = list(
    c(0, 14),                   # Prime: days 0-14
    c(20, 34),                  # Boost V2: days 20-34
    c(260, 274)                 # Boost V3: days 260-274
  )
)


tcr_collection <- list()
for (generation in 1:num_generations) {
  # Create TCR repertoire with viral parameters
  tcr_prime_boost <- new_tcr_primeBoostModel(sim_times, carry_cap, viral_params_prime_boost)

  # Add clones with different avidity levels
  # Persistent clones: high avidity (strong TCR-antigen binding)
  tcr_prime_boost <- add_clone_prime_boost_model(
    tcr = tcr_prime_boost,
    label = "persistent",
    init_size = clone_size,
    birth = rep(param_scale * 0.008, 9),
    death = rep(param_scale * 0.0032, 9),
    avidity = 0.5             # High avidity -> strong proliferation response
  )

  # Contracting clones: moderate avidity
  tcr_prime_boost <- add_clone_prime_boost_model(
    tcr = tcr_prime_boost,
    label = "contracting",
    init_size = clone_size,
    birth = rep(param_scale * 0.008, 9),
    death = rep(param_scale * 0.0042, 9),
    avidity = 0.1             # Moderate avidity -> moderate response
  )

  # Late emerging clones: very high avidity (strongest binders)
  tcr_prime_boost <- add_clone_prime_boost_model(
    tcr = tcr_prime_boost,
    label = "late_emerging",
    init_size = 3,
    birth = c(rep(0, 6), rep(param_scale * 0.010, 3)), # Late emerging pattern
    death = c(rep(0, 6), rep(param_scale * 0.0028, 3)),
    avidity = 0             # Very high avidity -> strongest response
  )

  # Run simulation
  cat("Running prime-boost simulation...\n")
  tcr_prime_boost$data <- tcr_simulate_prime_boost_model(tcr_prime_boost, benchmark = TRUE)

  # push tcr into tcr_collection
  tcr_collection <- append(tcr_collection, list(tcr_prime_boost))

  # print progress
  print(paste(round(generation * 100 / num_generations), "%"))
}

# save tcr object
saveRDS(tcr_collection, file.path(getwd(), data_dir, data_name))