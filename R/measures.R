# Simpson index
si <- function(smpl) {
  # avoid division by 0
  total_sum <- sum(smpl)
  if (total_sum <= 1) {
    return(0)
  }

  # Remove elements equal to 0
  smpl <- smpl[smpl != 0]

  index <- 1 - (sum(smpl / sum(smpl) * (smpl - 1) / (sum(smpl) - 1)))
  return(index)
}

# Morisita-Horn overlap
mh <- function(smpl_matrix) {
  result <- vegan::vegdist(smpl_matrix, method = "horn")
  result <- as.matrix(result)

  return(result)
}

# Chao overlap
ch <- function(smpl_matrix)  {
  result <- vegan::vegdist(smpl_matrix, method = "chao")
  result <- as.matrix(result)

  return(result)
}

# calculate measurements
get_measures <- function(tcr, smpl_table, draws = 10000, progress = FALSE) {
  # determine size nR (nR < nE) to which TCR set needs to be reduced
  smpl_min <- sample_min_size(smpl_table)

  # crete empty measure list
  measures <- list()

  if (progress) {
    # Create a progress bar
    pb <- utils::txtProgressBar(min = 0, max = ncol(smpl_table), style = 3)
  }

  for (i in seq_along(colnames(smpl_table))) {
    # clone id
    clo_id <- as.numeric(colnames(smpl_table)[i])

    # Simpson
    si_p1_d <- numeric(draws)
    si_s1_d <- numeric(draws)
    si_s2_d <- numeric(draws)
    # Morisita-Horn
    mh_s1p1 <- numeric(draws)
    mh_s2s1 <- numeric(draws)
    mh_s2p1 <- numeric(draws)
    mh_s1p1_d <- numeric(draws)
    mh_s2s1_d <- numeric(draws)
    mh_s2p1_d <- numeric(draws)
    # Chao-Jaccard
    ch_s1p1 <- numeric(draws)
    ch_s2s1 <- numeric(draws)
    ch_s2p1 <- numeric(draws)
    ch_s1p1_d <- numeric(draws)
    ch_s2s1_d <- numeric(draws)
    ch_s2p1_d <- numeric(draws)

    for (draw in seq(1, draws)) {
      # randomly draw, without replacement, nR TCRs from the set of nE TCRs
      draw_w <- sample_reduce(smpl_table, smpl_min, clone = clo_id, incl = T)
      draw_wo <- sample_reduce(smpl_table, smpl_min, clone = clo_id, incl = F)

      # calculate the div/ovlp measures for the reduced set
      # Simpson
      si_p1_d[draw] <- si(draw_w["P1", ]) - si(draw_wo["P1", ])
      si_s1_d[draw] <- si(draw_w["S1", ]) - si(draw_wo["S1", ])
      si_s2_d[draw] <- si(draw_w["S2", ]) - si(draw_wo["S2", ])

      # Morisita-Horn
      mh_w <- mh(draw_w)
      mh_wo <- mh(draw_wo)
      mh_s1p1[draw] <- mh_wo["S1", "P1"]
      mh_s2s1[draw] <- mh_wo["S2", "S1"]
      mh_s2p1[draw] <- mh_wo["S2", "P1"]
      mh_s1p1_d[draw] <- mh_w["S1", "P1"] - mh_wo["S1", "P1"]
      mh_s2s1_d[draw] <- mh_w["S2", "S1"] - mh_wo["S2", "S1"]
      mh_s2p1_d[draw] <- mh_w["S2", "P1"] - mh_wo["S2", "P1"]

      # Chao-Jaccard
      ch_w <- ch(draw_w)
      ch_wo <- ch(draw_wo)
      ch_s1p1[draw] <- ch_wo["S1", "P1"]
      ch_s2s1[draw] <- ch_wo["S2", "S1"]
      ch_s2p1[draw] <- ch_wo["S2", "P1"]
      ch_s1p1_d[draw] <- ch_w["S1", "P1"] - ch_wo["S1", "P1"]
      ch_s2s1_d[draw] <- ch_w["S2", "S1"] - ch_wo["S2", "S1"]
      ch_s2p1_d[draw] <- ch_w["S2", "P1"] - ch_wo["S2", "P1"]
    }

    # estimate the div/ovlp for the reduced TCR set from the median
    # Simpson
    measures$si_p1_d <- c(measures$si_p1_d, median(si_p1_d))
    measures$si_s1_d <- c(measures$si_s1_d, median(si_s1_d))
    measures$si_s2_d <- c(measures$si_s2_d, median(si_s2_d))

    # Morisita-Horn
    measures$mh_s1p1 <- c(measures$mh_s1p1, median(mh_s1p1))
    measures$mh_s2s1 <- c(measures$mh_s2s1, median(mh_s2s1))
    measures$mh_s2p1 <- c(measures$mh_s2p1, median(mh_s2p1))
    measures$mh_s1p1_d <- c(measures$mh_s1p1_d, median(mh_s1p1_d))
    measures$mh_s2s1_d <- c(measures$mh_s2s1_d, median(mh_s2s1_d))
    measures$mh_s2p1_d <- c(measures$mh_s2p1_d, median(mh_s2p1_d))

    # Chao-Jaccard
    measures$ch_s1p1 <- c(measures$ch_s1p1, median(ch_s1p1))
    measures$ch_s2s1 <- c(measures$ch_s2s1, median(ch_s2s1))
    measures$ch_s2p1 <- c(measures$ch_s2p1, median(ch_s2p1))
    measures$ch_s1p1_d <- c(measures$ch_s1p1_d, median(ch_s1p1_d))
    measures$ch_s2s1_d <- c(measures$ch_s2s1_d, median(ch_s2s1_d))
    measures$ch_s2p1_d <- c(measures$ch_s2p1_d, median(ch_s2p1_d))

    if (progress) {
      # update progress bar
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Convert measures list to a matrix for calculations
  measure_matrix <- do.call(cbind, measures)
  # add clone ids
  rownames(measure_matrix) <- colnames(smpl_table)
  # add clone labels
  label_vec <- sapply(colnames(smpl_table), function(id) {
    id <- as.numeric(id)
    return(as.numeric(tcr$clonotypes[[id]]$label))
  })

  measure_obj <- list(label_vec = label_vec,
                     measure_matrix = measure_matrix)

  if (progress) {
    # close progress bar
    close(pb)
  }

  # return measures
  return(measure_obj)
}

# scale measureMatrix
measure_scale <- function(measure, ...) {
  # store additional arguments
  args <- list(...)

  if (all(c("shift", "scale") %in% names(args))) {
    # shift and scale using user provided arguments
    args$scale <- ifelse(args$scale == 0, 1, args$scale) # avoid division by 0
    measure <- scale(measure, center = args$shift, scale = args$scale)

  } else if ("method" %in% names(args)) {
    # shift and scale based on method provided
    if (args$method == "minmax") {
      # scale to [0, 1]
      measure <- apply(measure, 2, function(x) {
        (x - min(x)) / (max(x) - min(x))
      })

    } else if (args$method == "robust") {
      # shift by median and scale by IQR
      median_measure <- apply(measure, 2, median)
      iqr_measure <- apply(measure, 2, IQR)

      # Avoid division by 0
      iqr_measure <- ifelse(iqr_measure == 0, 1, iqr_measure)

      measure <- scale(measure, center = median_measure, scale = iqr_measure)

    } else if (args$method == "z-score") {
      # shift and scale based on data provided
      sd_measure <- apply(measure, 2, sd)

      # Avoid division by 0
      sd_measure <- ifelse(sd_measure == 0, 1, sd_measure)

      # scale measures
      measure <- scale(measure, center = TRUE, scale = sd_measure)

    } else {
      # unknown method
      stop("Unknown method '", args$method, "'")
    }

  } else {
    # apply z-score normalization
    # shift and scale based on data provided
    sd_measure <- apply(measure, 2, sd)

    # Avoid division by 0
    sd_measure <- ifelse(sd_measure == 0, 1, sd_measure)

    # scale measures
    measure <- scale(measure, center = TRUE, scale = sd_measure)
  }

  return(measure)
}
