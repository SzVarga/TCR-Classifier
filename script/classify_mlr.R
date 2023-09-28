# clean env
rm(list = ls())

# set seed
set.seed(20220314)

# load libraries and functions
require(nnet)        # multinomial logistic regression
library(caret)       # confMatrix eval
source("R/tcr_model.R")
source("R/measures.R")
source("R/classifier.R")

# load tcr_collection object
tcr_collection <- readRDS("data/tcr_data/tcrColl_main.rds")

# load reference data
measure_ref <- readRDS("data/ref_measures/main_measure_ref.rds")

# CLI arguments
args <- commandArgs(trailingOnly = TRUE)

# settings
iterators <- 1:10
draw <- 10
smpl_size <- as.numeric(args[1])

# data container
classif_name <- "mlr"
clo_labels <- tcr_collection[[1]]$clone_labels
rep <- c()
sample_size <- c()
clo_type <- c()
classif <- c()
accuracy <- c()
tp_rate <- c()
tn_rate <- c()
precision <- c()
jackknife_only <- TRUE

# Train the multinomial regression model
# keep only measures calculated using d1-jackknife
if (jackknife_only) {
  idx_jackknife <- grep("_d", colnames(measure_ref$measure_matrix))
  measure_ref$measure_matrix <- measure_ref$measure_matrix[, idx_jackknife]
}

# scale reference measure_matrix
measure_ref$measure_matrix <- measure_scale(measure_ref$measure_matrix)

# create reference data frame structure
ref_data <- as.data.frame(measure_ref$measure_matrix)
ref_data$label <- measure_ref$label_vec

# relevel labels
ref_data$label <- as.factor(ref_data$label)
ref_data$label <- relevel(ref_data$label, ref = 1)

# train model
model <- nnet::multinom(label ~ ., data = ref_data, maxit = 1000)

# iterate over tcr objects (for statistical power)
# and predict clone type from sample data
for (iterator in iterators) {
  # retreive a single tcr object
  tcr <- tcr_collection[[iterator]]

  #sample data
  starttime <- Sys.time()
  smpl <- list(
    P1 = sample_at(tcr, smpl_size, "P1", detect_lim = 3),
    S1 = sample_at(tcr, smpl_size, "S1", detect_lim = 3),
    S2 = sample_at(tcr, smpl_size, "S2", detect_lim = 3))

  # combine samples
  smpl_pred <- do.call(sample_bind, smpl)

  # calculate measures
  measure_pred <- get_measures(tcr, smpl_pred, draws = draw, progress = TRUE)

  # keep only measures calculated using d1-jackknife
  if (jackknife_only) {
    idx_jackknife <- grep("_d", colnames(measure_pred$measure_matrix))
    measure_pred$measure_matrix <- measure_pred$measure_matrix[, idx_jackknife]
  }

  # scale prediction measure_matrix
  measure_pred$measure_matrix <- measure_scale(measure_pred$measure_matrix)

  # create predicion data frame structure
  pred_data <- as.data.frame(measure_pred$measure_matrix)

  # predict clone type based on multinomial regression model
  prediction <- t(predict(model, newdata = pred_data, type = "prob"))

  # use most likely prediction label
  binary_model <- c("persistent", "contracting")
  full_model <- c("persistent", "contracting", "late_emerging")
  if (identical(tcr$clone_labels, binary_model)) {
    # binary model - estimation values are on [0,1]
    estimation <- apply(prediction, 2, function(p) round(p) + 1)
  } else if (identical(tcr$clone_labels, full_model)) {
    # full model - estimation values are prob vectors of length 3
    estimation <- apply(prediction, 2, which.max)
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

# collect classifier performance data
data <- data.frame(tcr = rep,
                   sample_size = sample_size,
                   clone_label = clo_type,
                   classifier = classif,
                   accuracy = accuracy,
                   sensitivity = tp_rate,
                   specificity = tn_rate,
                   precision = precision)


# save classifier performance data
saveRDS(data, paste0("data/mlr/perf_mlr_", smpl_size, ".rds"))
