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

data_in <- args[1]                        # TCR collection data
smpl_size <- ifelse(
  args[2] == "all",
  args[2],
  as.numeric(args[2]))                    # sample size
data_out <- "data"                        # output directory
detect_lim <- 3                           # clone detection limit
draws <- 10000                            # resampling draws
progress <- FALSE                         # show progress bar


# load tcr_collection object
tryCatch({
  tcr_collection <- readRDS(data_in)
}, error = function(e) {
  stop("Error reading data: ", e$message)
})


# main loop over all tcr objects
for (i in seq_along(tcr_collection)) {

  # retreive a single tcr object
  tcr <- tcr_collection[[i]]

  # reference data
  smpl <- list(P1 = sample_at(tcr, smpl_size, "P1", detect_lim = detect_lim),
               S1 = sample_at(tcr, smpl_size, "S1", detect_lim = detect_lim),
               S2 = sample_at(tcr, smpl_size, "S2", detect_lim = detect_lim))

  # combine ref data
  smpl <- do.call(sample_bind, smpl)

  # calculate measures
  measures <- get_measures(tcr, smpl, draws = draws, progress = progress)

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
    sample_size <- rep(smpl_size, nrow(subdata))
    tcr_id <- rep(i, nrow(subdata))
    label <- sapply(clone, function(x) {tcr$clonotypes[[x]]$label})
    subdata <- cbind(tcr_id, clone, sample_size, subdata, label)

    # remove rownames
    rownames(subdata) <- seq(1, nrow(subdata), 1)
  }

  # save output
  # Define the file path
  file_path <- file.path(data_out, "sample_table.csv")

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