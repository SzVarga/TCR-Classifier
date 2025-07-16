# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# load libraries and functions
require(adaptivetau)    # for stochastic model simulation
require(sn)             # for skewed normal distribution
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

tcr_collection <- list()
for (generation in 1:num_generations) {
  # create a tcr object
  tcr <- new_tcr(sim_times = c("P10"=10, "V2"=10, "S10"=10, "S68"=58, "S210"=152, "V3"=20, "T10"=10, "T108"=98, "T189"=79),
                 carry_cap = num_clones * clone_size / 0.6)

  # add clones to TCR-repertoire
  for (i in 1:num_clones) {
    # sample birth and death rates from skewed normal distribution
    # i-th persistent clone
    birth_pers_val <- sample_skewed_normal(param_scale * 0.008, param_scale * 0.00008, 0)
    death_pers_val <- sample_skewed_normal(param_scale * 0.0032, param_scale * 0.000032, 0)

    # i-th contracting clone
    birth_cont_val <- sample_skewed_normal(param_scale * 0.008, param_scale * 0.00008, 0)
    death_cont_val <- sample_skewed_normal(param_scale * 0.0042, param_scale * 0.000042, 0)

    # i-th late-emerging clone
    birth_late_val <- sample_skewed_normal(param_scale * 0.010, param_scale * 0.0001, 0)
    death_late_val <- sample_skewed_normal(param_scale * 0.0028, param_scale * 0.000028, 0)

    tcr <- add_clone(tcr = tcr, label = "persistent",
                     init_size = clone_size,
                     birth = rep(birth_pers_val, 9),
                     death = rep(death_pers_val, 9))
    tcr <- add_clone(tcr = tcr, label = "contracting",
                     init_size = clone_size,
                     birth = rep(birth_cont_val, 9),
                     death = rep(death_cont_val, 9))
    tcr <- add_clone(tcr = tcr, label = "late_emerging",
                     init_size = 1,
                     birth = c(rep(0, 5), rep(birth_late_val, 4)),
                     death = c(rep(0, 5), rep(death_late_val, 4)))
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
