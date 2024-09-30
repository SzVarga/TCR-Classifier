# clean env
rm(list = ls())

# load libraries and functions
library(filelock)
source("R/tcr_model.R")
source("R/measures.R")

# set seed
set.seed(42)

# CLI arguments
args <- commandArgs(trailingOnly = TRUE)

DATA_IN <- args[1]                        # TCR collection data
TCR_ID <- 1                               # tcr id to sample from
SMPL_SIZE <- ifelse(
  args[2] == "all",
  args[2],
  as.numeric(args[2]))                    # sample size
N_SMPL_LOOP <- 30                         # number of replicates
DATA_OUT <- "data"                        # output directory
DETECT_LIM <- 3                           # clone detection limit
DRAWS <- 10000                            # resampling draws
PROGRESS <- FALSE                         # show progress bar

tryCatch({
  # load tcr_collection object
  tcr_collection <- readRDS(DATA_IN)

  # retreive a single tcr object
  tcr <- tcr_collection[[TCR_ID]]
}, error = function(e) {
  stop("Error reading data: ", e$message)
})

# main loop over all tcr objects
for (i in seq(1, N_SMPL_LOOP, 1)) {

  # reference data
  smpl <- list(P1 = sample_at(tcr, SMPL_SIZE, "P1", detect_lim = DETECT_LIM),
               S1 = sample_at(tcr, SMPL_SIZE, "S1", detect_lim = DETECT_LIM),
               S2 = sample_at(tcr, SMPL_SIZE, "S2", detect_lim = DETECT_LIM))

  # combine ref data
  smpl <- do.call(sample_bind, smpl)

  # calculate measures
  measures <- get_measures(tcr, smpl, draws = DRAWS, progress = PROGRESS)

  # create subdata by transposing the sample matrix
  subdata <- t(smpl)

  # check if the rownames of the measure matrix and the subdata match
  if (any(rownames(measures$measure_matrix) != rownames(subdata))) {
    stop("Error: rownames of measure matrix and of subdata do not match.")
  } else {
    #combine subdata and measure matrix
    subdata <- cbind(subdata, measures$measure_matrix)

    # add additional information
    clone <- as.numeric(rownames(subdata))
    sample_size <- rep(SMPL_SIZE, nrow(subdata))
    tcr_id <- rep(TCR_ID, nrow(subdata))
    label <- sapply(clone, function(x) {tcr$clonotypes[[x]]$label})
    subdata <- cbind(tcr_id, clone, sample_size, subdata, label)

    # remove rownames
    rownames(subdata) <- seq(1, nrow(subdata), 1)
  }

  # save output
  # Define the file path
  file_path <- file.path(DATA_OUT, "sample_table.csv")

  # Write or append the subdata to the CSV file
  if (!file.exists(file_path)) {
    # stop and throw an error saying that the file does not exist
    stop("Error: The output file sample_table.csv does not exist.")
  } else {
    # Lock the file
    lock <- filelock::lock(file_path)

    # check if the file is empty
    if (file.size(file_path) == 0) {
      # Write the header
      write.table(
        subdata[0,], file_path, sep = ",",
        col.names = TRUE, row.names = FALSE
      )
    }

    # Append the data
    write.table(
      subdata, file_path, sep = ",",
      col.names = FALSE, row.names = FALSE, append = TRUE
    )

    # Unlock the file
    filelock::unlock(lock)
  }
}