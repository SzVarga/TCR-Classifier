# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# load libraries and functions
require(caret)          # confusion matrix evaluation
source("R/tcr_model.R")
source("R/measures.R")
source("R/classifier.R")

# load tcr_collection object
tcr_collection <- readRDS("data/tcr_data/tcrColl_main.rds")

# CLI arguments
args <- commandArgs(trailingOnly = TRUE)

# settings
iterators <- 1:10
draw <- 1000
smpl_size <- as.numeric(args[1])

# data container
classif_name <- "naive"
clo_labels <- tcr_collection[[1]]$clone_labels
rep <- c()
sample_size <- c()
clo_type <- c()
classif <- c()
accuracy <- c()
tp_rate <- c()
tn_rate <- c()
precision <- c()


for (iterator in iterators) {
  # retreive a single tcr object
  tcr <- tcr_collection[[iterator]]

  #sample data
  starttime <- Sys.time()
  smpl <- list(
    P1 = sample_at(tcr, smpl_size, "P1", detect_lim = 3),
    S1 = sample_at(tcr, smpl_size, "S1", detect_lim = 3),
    S2 = sample_at(tcr, smpl_size, "S2", detect_lim = 3))

  # predict labels and generate confusion matrix
  binary_model <- c("persistent", "contracting")
  full_model <- c("persistent", "contracting", "late_emerging")
  if (identical(tcr$clone_labels, binary_model)) {
    # binary model
    estimation <- naive_classifier(smpl$S1, smpl$S2, binary = TRUE)
  } else if (identical(tcr$clone_labels, full_model)) {
    # full model
    estimation <- naive_classifier(smpl$S1, smpl$S2, binary = FALSE)
  } else {
    stop("unexpected clone labels")
  }
  conf_matrix <- evaluate_prediction(tcr, estimation)

  # expand conf matrix to be 3x3
  mat <- matrix(data = 0, nrow = 3, ncol = 3)
  mat[seq_len(dim(conf_matrix)[1]), seq_len(dim(conf_matrix)[2])] <- conf_matrix
  # name rows and columns
  colnames(mat) <- c(1, 2, 3)
  rownames(mat) <- c(1, 2, 3)
  # evaluate the confusion matrix
  conf_mat <- caret::confusionMatrix(mat)

  # append data vectors
  for (clo_label in seq_along(clo_labels)) {
    rep <- c(rep, iterator)
    sample_size <- c(sample_size, smpl_size)
    clo_type <- c(clo_type, clo_labels[clo_label])
    classif <- c(classif, classif_name)
    accuracy <- c(accuracy, conf_mat$overall[["Accuracy"]])
    tp_rate <- c(tp_rate, conf_mat$byClass[clo_label, "Sensitivity"])
    tn_rate <- c(tn_rate, conf_mat$byClass[clo_label, "Specificity"])
    precision <- c(precision, conf_mat$byClass[clo_label, "Precision"])
  }

  # print required simulation time
  endtime <- Sys.time()
  sim_time <- difftime(endtime, starttime, units = "secs")
  cat("simulation time:", sim_time, "s", "\n")
}

data <- data.frame(tcr = rep,
                   sample_size = sample_size,
                   clone_label = clo_type,
                   classifier = classif,
                   accuracy = accuracy,
                   sensitivity = tp_rate,
                   specificity = tn_rate,
                   precision = precision)


# save classifier performance data
saveRDS(data, paste0("data/naive/perf_naive_", smpl_size, ".rds"))
