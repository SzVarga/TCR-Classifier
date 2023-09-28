#' Naive Classifier
#'
#' This function performs a naive classification of T cell clones
#' based on their presence at two time points. It can classify clones
#' into three categories: persisting, contracting, or late emerging.
#' Alternatively, you can use a binary classification that only distinguishes
#' between persisting and contracting clones.
#'
#' @param before Numeric vector representing the clone sizes at first time
#'               point. The names of the vector should correspond to clone IDs.
#' @param after Numeric vector representing the clone sizes at second time
#'              point. The names of the vector should correspond to clone IDs.
#' @param binary Logical value indicating whether to use binary classification.
#'        If TRUE, clones are classified as either persisting (1) or
#'        contracting (2).
#'        If FALSE, three categories are used: persisting (1), contracting (2),
#'        and late emerging (3).
#'
#' @return Named numeric vector of classification results. Clone IDs are
#'         represented by names, the classification results are represented
#'         by the values.
#'
#' @export
naive_classifier <- function(before, after, binary = FALSE) {

  all_ids <- union(names(before), names(after))
  result <- numeric(length(all_ids))

  # predict labels
  for (i in seq_along(all_ids)) {
    id <- all_ids[i]

    # has clone been sampled in before and after
    in_before <- id %in% names(before)
    in_after <- id %in% names(after)

    if (in_before && in_after) {
      # classify as persisting
      result[i] <- 1
    } else if (in_before) {
      # classify as contracting
      result[i] <- 2
    } else if (in_after) {
      if (!binary) {
        # classify as late emerging
        result[i] <- 3
      } else {
        # classify as persisting
        result[i] <- 1
      }
    }
  }

  names(result) <- all_ids
  return(result)
}


#' Perform Principal Component Analysis (PCA)
#'
#' This function performs Principal Component Analysis (PCA)
#' on a matrix of diversity or overlap measures to allow for
#' dimensionality reduction and feature extraction based on
#' principal components.
#'
#' @param measures A matrix of diversity or overlap measures where
#' rows represent clones, and columns represent indices.
#' 
#' @param center A logical value indicating whether to center
#' the input data before performing PCA. Default is FALSE.
#' @param scale A logical value indicating whether to scale the
#' input data before performing PCA. Default is FALSE.
#' 
#' @details Shifting and scaling of the data is performed by a
#' separate function.
#'
#' @return An object of class 'prcomp' representing the PCA results.
#'
#' @export
get_pca <- function(measures, center = FALSE, scale = FALSE) {
  pc <- prcomp(measures, center = center, scale = scale)
  return(pc)
}

#' Create a PCA Mapper Function
#'
#' This function creates a PCA (Principal Component Analysis) mapper function
#' based on a reference measures. The mapper function can be used to project
#' new datasets onto a PCA space based on the reference dataset's
#' principal components.
#'
#' @param measures_ref Reference measures, including a matrix of measures
#'        where rows represent clones, and columns represent different indices.
#' 
#' @param pc_x Index of first PC to use in the mapping function.
#' @param pc_y Index of second PC to use in the mapping function.
#' @param pcx_features Number of top loadings to consider for the first PC.
#' @param pcy_features Number of top loadings to consider for the second PC.
#'
#' @return A list containing the following elements:
#'  - transform : The mapping function that can be called to map new data.
#'  - pc : The results of the PCA performed on the reference dataset.
#'  - x_label : The formula (string) for the x-axis values.
#'  - y_label : The formula (string) for the y-axis values.
#'
#' @export
create_pca_mapper <- function(measures_ref,
                              pc_x, pc_y,
                              pcx_features, pcy_features) {

  # perform pca
  pc <- get_pca(measures_ref$measure_matrix)

  # define principal components to use
  pc_x <- paste0("PC", pc_x)
  pc_y <- paste0("PC", pc_y)

  # get the first n (biggest) loadings
  idx_x <- order(abs(pc$rotation[, pc_x]), decreasing = TRUE)[1:pcx_features]
  idx_y <- order(abs(pc$rotation[, pc_y]), decreasing = TRUE)[1:pcy_features]

  # extract loadings
  loadings <- list()
  loadings$x <- pc$rotation[idx_x, pc_x]
  loadings$y <- pc$rotation[idx_y, pc_y]

  # create transformation labels
  x_lab_comp <- character(length(loadings$x))
  y_lab_comp <- character(length(loadings$x))
  for (i in seq_along(loadings$x)) {
    x_lab_comp[i] <- paste0(round(loadings$x[[i]], digits = 3),
                            "*", names(loadings$x)[i])
  }
  for (i in seq_along(loadings$y)) {
    y_lab_comp[i] <- paste0(round(loadings$y[[i]], digits = 3),
                            "*", names(loadings$y)[i])
  }
  # concatenate label components
  x_lab <- paste(x_lab_comp, collapse = " + ")
  y_lab <- paste(y_lab_comp, collapse = " + ")

  # Define the mapping function
  pca_mapping_function <- function(measures) {
    # create data structure
    data <- list(id = rownames(measures$measure_matrix),
                 label = measures$label_vec,
                 x = numeric(nrow(measures$measure_matrix)),
                 y = numeric(nrow(measures$measure_matrix)))
    #convert measures to data frame
    measures <- as.data.frame(measures$measure_matrix)
    # map data
    for (i in seq_along(loadings$x)) {
      data$x <- data$x +
                loadings$x[[i]] * measures[, names(loadings$x)[i]]
    }
    for (i in seq_along(loadings$y)) {
      data$y <- data$y +
                loadings$y[[i]] * measures[, names(loadings$y)[i]]
    }
    return(data)
  }

  # return result object
  result <- list(transform = pca_mapping_function,
                 pc = pc,
                 x_label = x_lab,
                 y_label = y_lab)
  return(result)
}

##' Create a scatter plot of PCA-transformed data.
#'
#' This function generates a scatter plot of PCA-transformed data points. It can
#' also overlay predicted data points on the same plot if pred_data is provided.
#'
#' @param ref_data Reference data containing PCA-transformed values.
#' @param pred_data Optional predicted data to overlay on the plot.
#' @param legend Position of the legend in the plot (e.g., "bottomright").
#' @param ... Additional arguments to be passed to the plot.
#'
#' @export
plt_pca_features <- function(ref_data, pred_data = NULL,
                             legend = "bottomright", ...) {
  # picking color based on true clone label
  map_to_colors <- function(vector) {
    result <- sapply(vector, function(i) {
      switch(i,
             "1" = "#8DC955",
             "2" = "#DE5A43",
             "3" = "#33528F",
             "black"  # Default color when 'i' does not match any case
      )
    })
    return(result)
  }

  # plotting
  plot(ref_data$x, ref_data$y, type = "p",
       main = "Scatterplot", col = map_to_colors(ref_data$label),
       xlim = c(min(ref_data$x, pred_data$x), max(ref_data$x, pred_data$x)),
       ylim = c(min(ref_data$y, pred_data$y), max(ref_data$y, pred_data$y)),
       ...)
  if (!is.null(pred_data)) {
    points(pred_data$x, pred_data$y, type = "p", pch = 4,
           col = map_to_colors(pred_data$label))
  }
  if (legend != "none") {
    legend(legend, legend = c("Reference", "Sample",
                              "Persistent", "Contracting", "Late Emerging"),
           col = c("black", "black", "#8DC955", "#DE5A43", "#33528F"),
           pch = c(1, 4, 20, 20, 20))
  }
}

#' Create a heatmap of PCA loadings.
#'
#' This function generates a heatmap of the absolute values
#' of PCA loadings for the provided measures.
#'
#' @param measures A matrix of measures to compute PCA loadings from.
#' @param scale Logical, indicating whether to scale the measures before PCA.
#'
plt_pca_loadings <- function(measures, scale = FALSE) {
  # scale measures
  if (scale) {
    measures <- measure_scale(measures)
  }

  # perform pca
  pc <- get_pca(measures)

  # create heatmap of loadings
  plt <- heatmap(abs(pc$rotation), Rowv = NA, Colv = NA)

  return(plt)
}

#' Predict clone labels using k-Nearest Neighbors (k-NN).
#'
#' This function predicts clone labels for a given set of samples
#' using k-NN statistics.
#'
#' @param tcr The TCR object containing information about clones and labels.
#' @param reference_data The reference dataset used for k-NN classification.
#' @param prediction_data The dataset for which clone labels will be predicted.
#' @param k The number of nearest neighbors to consider for classification.
#'
#' @return A vector of predicted clone labels.
#'
#' @export
knn_predict <- function(tcr, reference_data, prediction_data, k = 1, ...) {
  # get indices for closest clones
  knn_index <- FNN::get.knnx(reference_data,
                             prediction_data,
                             k = 1, ...)$nn.index

  # name knn_index rows and columns
  rownames(knn_index) <- rownames(prediction_data)
  colnames(knn_index) <- seq(1, ncol(knn_index))

  # predict labels based on knn statistics
  estimation <- apply(knn_index, 1, function(x) {
    ids <- rownames(reference_data)[x]
    ids <- as.numeric(ids)
    labels <- sapply(ids, function(id) return(tcr$clonotypes[[id]]$label))
    label <- which.max(table(labels))
    return(label)
  })

  return(estimation)
}

#' Evaluate the predicted clone labels against the true labels.
#'
#' This function creates the confusion matrix by comparing the
#' predicted clone labels to the true labels.
#'
#' @param repertoire The TCR repertoire object containing true clone labels.
#' @param estimated_labels A named vector of predicted clone labels.
#'
#' @return A confusion matrix comparing estimated labels to true labels.
#'
#' @export
# return the confusion matrix
evaluate_prediction <- function(repertoire, estimated_labels) {
  ids <- as.integer(names(estimated_labels))
  true_labels <- sapply(ids, function(id) {
    as.integer(repertoire$clonotypes[[id]]$label)
  })
  estimated_labels <- as.numeric(estimated_labels)

  confusion_matrix <- table(estimated = estimated_labels, truth = true_labels)
  return(confusion_matrix)
}
