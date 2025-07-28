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

get_logistic_death <- function(birth_rate, carry_cap, eff_cap = 0.95) {
  # Calculate the logistic death rate based on birth rate and carrying capacity
  death_rate <- birth_rate * (1 - (carry_cap * eff_cap) / carry_cap)
  return(death_rate)
}

# Set simulation parameters
bip <- 50
sim_times <- c("burn-in" = bip, "P10" = 10, "V2" = 10, "S10" = 10, "S68" = 58,
               "S210" = 152, "V3" = 20, "T10" = 10, "T108" = 98, "T189" = 79)
param_scale <- 100
num_generations <- 1
num_clones <- 100
clone_size <- 100
carry_cap <- num_clones * clone_size * 2 # Mainly persistent and late emerging clones
data_dir <- "data/tcr_data/"
data_name <- paste0("tcrColl_", num_clones, "clo_a", clone_size, "_x",
                    param_scale, ".rds")

# Example 1: Prime-boost vaccination schedule
cat("Example 1: Prime-boost vaccination schedule\n")
cat("=========================================\n")

# Define viral parameters for prime-boost
viral_params_prime_boost <- list(
  initial_load = 0,             # dV/dt = r(t) × V - c × T_total × V
  replication_rate = 4E-1,      # log(V) ~ 2
  clearance_rate = 8E-6,        # 100clones x 100cells x 1E2 viruses ~ 1E6
  replication_intervals = list(
    c(0, 7)+bip,                    # Prime: days 0-7 + burn-in-phase
    c(20, 27)+bip,                    # Boost V2: days 20-27
    c(260, 267)+bip                   # Boost V3: days 260-267
  ),
  viral_burden_sensitivity = 0  # Enable cumulative viral burden effect
)

tcr_collection <- list()
for (generation in 1:num_generations) {
  # Create TCR repertoire with viral parameters
  tcr_prime_boost <- new_tcr_primeBoostModel(sim_times, carry_cap, viral_params_prime_boost)

  for (i in 1:num_clones) {
    # Randomly select clone parameters
    birth_persistent <- sample_skewed_normal(param_scale * 0.0008, param_scale * 0.0008, 0)
    death_persistent <- get_logistic_death(birth_persistent, carry_cap)
    birth_contracting <- sample_skewed_normal(param_scale * 0.0008, param_scale * 0.00008, 0)
    death_contracting <- get_logistic_death(birth_contracting, carry_cap) * sample_skewed_normal(3, 0.5, 0)
    birth_late_emerging <- sample_skewed_normal(param_scale * 0.0008, param_scale * 0.0008 * 1e-2, 0) * sample_skewed_normal(6, 2, 0)
    death_late_emerging <- get_logistic_death(birth_late_emerging, carry_cap) * sample_skewed_normal(0.01, 0.001, 0)

    # Add clones with different avidity levels
    # Persistent clones: high avidity (strong TCR-antigen binding)
    tcr_prime_boost <- add_clone_prime_boost_model(
      tcr = tcr_prime_boost,
      label = "persistent",
      init_size = clone_size,
      birth = c(birth_persistent, rep(birth_persistent, 9)),
      death = c(get_logistic_death(birth_persistent, carry_cap), rep(death_persistent, 9)),
      avidity = 0
    )

    # Contracting clones: moderate avidity
    tcr_prime_boost <- add_clone_prime_boost_model(
      tcr = tcr_prime_boost,
      label = "contracting",
      init_size = clone_size,
      birth = c(birth_contracting, rep(birth_contracting, 9)),
      death = c(get_logistic_death(birth_contracting, carry_cap), rep(death_contracting, 9)),
      avidity = 0
    )

    # Late emerging clones: very high avidity (strongest binders)
    tcr_prime_boost <- add_clone_prime_boost_model(
      tcr = tcr_prime_boost,
      label = "late_emerging",
      init_size = 3,
      birth = c(rep(0, 7), rep(birth_late_emerging, 3)), # Late emerging pattern
      death = c(rep(0, 7), rep(death_late_emerging, 3)),
      avidity = 0
    )
  }

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