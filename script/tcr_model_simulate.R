# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# load libraries and functions
require(adaptivetau)    # for stochastic model simulation
source("R/tcr_model.R") # for tcr model handling

# get command line argument provided by the sbatch file
configs <- commandArgs(trailingOnly = TRUE)

num_generations <- as.integer(configs[1])
num_clones <- as.integer(configs[2])
clone_size <- as.integer(configs[3])
param_scale <- as.numeric(configs[4])
data_dir <- "data/tcr_data/"
data_name <- paste0("tcrColl_", num_clones, "clo_a", clone_size, "_x",
                    param_scale, ".rds")

# parameters for clonal dynamics
birth_pers <- param_scale * c(0.008, 0.008, 0.008)
death_pers <- param_scale * c(0.0032, 0.0032, 0.0032)
birth_cont <- param_scale * c(0.008, 0.008, 0.008)
death_cont <- param_scale * c(0.0042, 0.0042, 0.0042)
birth_late <- param_scale * c(0, 0, 0.010)
death_late <- param_scale * c(0, 0, 0.0028)

tcr_collection <- list()
for (generation in 1:num_generations) {
    # create a tcr object
    tcr <- new_tcr(sim_times = c("P1" = 10, "S1" = 10, "S2" = 10),
                  carry_cap = num_clones * clone_size / 0.6)


    # add clones to TCR-repertoire
    for (i in 1:num_clones) {
        tcr <- add_clone(tcr = tcr, label = "persistent",
                         init_size = clone_size,
                         birth = birth_pers,
                         death = death_pers)
        tcr <- add_clone(tcr = tcr, label = "contracting",
                         init_size = clone_size,
                         birth = birth_cont,
                         death = death_cont)
        tcr <- add_clone(tcr = tcr, label = "late_emerging",
                         init_size = 1,
                         birth = birth_late,
                         death = death_late)
  }

  # perform tcr-repertoire simulation
  tcr$data <- tcr_simulate(repertoire = tcr)

  # push tcr into tcr_collection
  tcr_collection <- append(tcr_collection, list(tcr))

  # print progress
  print(paste(round(generation * 100 / num_generations), "%"))
}

# save tcr object
saveRDS(tcr_collection, file.path(getwd(), data_dir, data_name))
