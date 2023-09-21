# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# load libraries and functions
source("R/tcr_model.R")
source("R/measures.R")

# load tcr_collection object
tcr_collection <- readRDS("data/tcr_data/tcrColl_10clo_a1000_x100.rds")

# retreive a single tcr object
tcr <- tcr_collection[[1]]

# CLI arguments
args <- commandArgs(trailingOnly = TRUE)

# set sample size
sample_size <- as.numeric(args[1])
iteration <- as.numeric(args[2])

# reference data
smpl <- list(P1 = sample_at(tcr, sample_size, "P1", detect_lim = 3),
             S1 = sample_at(tcr, sample_size, "S1", detect_lim = 3),
             S2 = sample_at(tcr, sample_size, "S2", detect_lim = 3))

# combine ref data
smpl_ref <- do.call(sample_bind, smpl)

# calculate measures
ref_measure <- get_measures(tcr, smpl_ref, draws = 10000, progress = TRUE)

# save output
dir_out <- "data/ref_measures"
file_name <- paste0("ref_measure_", sample_size, "_", iteration, ".rds")
saveRDS(ref_measure, file.path(dir_out, file_name))
