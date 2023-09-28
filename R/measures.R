#' Calculate the Simpson's Diversity Index for a sample.
#'
#' This function quantifies the probability that two individuals
#' randomly selected from a sample will belong to different species.
#' It ranges from 0 (no diversity) to 1 (maximum diversity).
#'
#' @param smpl A numeric vector representing the abundance of each species
#' in the sample.
#' @return The Simpson's Diversity Index for the sample, ranging from 0 to 1.
#'
#' @export
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

#' Calculate the Morisita-Horn dissimilarity matrix for a sample.
#'
#' This function calculates the Morisita-Horn dissimilarity matrix for a sample
#' using the package \code{vegan}.
#'
#' @param smpl_matrix A numeric matrix where each row represents a time point,
#'        and each column represents the abundance of a particular clone.
#' @return A dissimilarity matrix calculated using the Morisita-Horn method.
#'
#' @export
mh <- function(smpl_matrix) {
  result <- vegan::vegdist(smpl_matrix, method = "horn")
  result <- as.matrix(result)

  return(result)
}

#' Calculate the Chao-Jaccard dissimilarity matrix for a sample matrix.
#'
#' This function calculates the Chao-Jaccard dissimilarity matrix for a sample
#' using the package \code{vegan}.
#'
#' @param smpl_matrix A numeric matrix where each row represents a time point,
#'        and each column represents the abundance of a particular clone.
#' @return A dissimilarity matrix calculated using the Chao-Jaccard method.
#'
#' @export
ch <- function(smpl_matrix)  {
  result <- vegan::vegdist(smpl_matrix, method = "chao")
  result <- as.matrix(result)

  return(result)
}

#' Calculate diversity and overlap measures for a samples.
#'
#' This function calculates diversity and overlap measures,
#' including Simpson's Index, Morisita-Horn Index, and Chao-Jaccard Index,
#' for clones sampled at different time points. It estimates these measures
#' for each clone in the repertoire based on downsampling of the sample table,
#' to avoid sampling bias.
#'
#' @param tcr A TCR repertoire object.
#' @param smpl_table A sample table representing the abundance of each clone
#' in the samples. Rows correspond to different time points, and columns
#' correspond to clone IDs.
#' @param draws The number of random draws to generate the resampled
#' abundance data.
#' @param progress A logical value indicating whether to display a progress bar.
#'
#' @return A list containing diversity and overlap measures for each clone.
#'        The list includes the following components:
#'        - label_vec: A vector of clone labels corresponding to each clone.
#'        - measure_matrix: A matrix containing the calculated diversity and
#'        overlap measures for each clone.
#'
#' @export
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

#' Scale and shift the measures matrix.
#'
#' This function scales and shifts the measures matrix using various methods.
#' You can specify the method for scaling and shifting, or provide custom
#' values. The available methods are "minmax" (scale to [0, 1]),
#' "robust" (shift by median and scale by IQR), and "z-score"
#' (center at 0 and scale by standard deviation). If no method or
#' custom values are provided, the function applies z-score normalization
#' by default. If a column has 0 standard deviation or IQR, the function will
#' not scale the corresponding column.
#'
#' @param measure A matrix of measures where rows represent observations and
#'        columns represent measures (e.g., Morisita-Horn).
#' @param ... Additional arguments that can be provided to customize the scaling
#'        and shifting. You can use the following arguments:
#'        - method: A character string specifying the scaling and shifting
#'        method ("minmax", "robust", or "z-score"). Default is "z-score".
#'        - shift: Vector specifying the shift to be applied to each column.
#'        - scale: Vector specifying the scaling factor to be applied.
#'
#' @return A matrix of scaled and shifted measures.
#'
#' @export
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
