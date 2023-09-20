# naive classifier based on clonal incidence
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


# principle component analysis
get_pca <- function(measures, center = FALSE, scale = FALSE) {
  pc <- prcomp(measures, center = center, scale = scale)
  return(pc)
}


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


# classifier
get_pca_features <- function(measures_ref, measures_pred,
                             pc_x = 1, pc_y = 2,
                             pcx_features = 2, pcy_features = 2) {

  # perform pca
  pc <- get_pca(measures_ref$measure_matrix)

  # define principal components to use
  pc_x <- paste0("PC", pc_x)
  pc_y <- paste0("PC", pc_y)

  # get loadings for principal components
  loadings <- list()
  loadings$x <- pc$rotation[order(abs(pc$rotation[, pc_x]), decreasing = TRUE)[1:pcx_features], pc_x] 
  loadings$y <- pc$rotation[order(abs(pc$rotation[, pc_y]), decreasing = TRUE)[1:pcy_features], pc_y]


  # create reference data set
  ref_data <- list(id = rownames(measures_ref$measure_matrix),
                  label = measures_ref$label_vec,
                  x = numeric(nrow(measures_ref$measure_matrix)),
                  y = numeric(nrow(measures_ref$measure_matrix)))
  measures_ref <- as.data.frame(measures_ref$measure_matrix)
  for (i in seq_along(loadings$x)) {
    ref_data$x <- ref_data$x +
                  loadings$x[[i]] * measures_ref[, names(loadings$x)[i]]
  }
  for (i in seq_along(loadings$y)) {
    ref_data$y <- ref_data$y +
                  loadings$y[[i]] * measures_ref[, names(loadings$y)[i]]
  }

  # create prediction data set
  pred_data <- list(id = rownames(measures_pred$measure_matrix),
                   label = measures_pred$label_vec,
                   x = numeric(nrow(measures_pred$measure_matrix)),
                   y = numeric(nrow(measures_pred$measure_matrix)))
  measures_pred <- as.data.frame(measures_pred$measure_matrix)
  for (i in seq_along(loadings$x)) {
    pred_data$x <- pred_data$x +
                   loadings$x[[i]] * measures_pred[, names(loadings$x)[i]]
  }
  for (i in seq_along(loadings$y)) {
    pred_data$y <- pred_data$y +
                   loadings$y[[i]] * measures_pred[, names(loadings$y)[i]]
  }

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

  result <- list(ref_data = ref_data,
                 pred_data = pred_data,
                 pc_ref = pc,
                 x_label = x_lab,
                 y_label = y_lab)
  return(result)
}

# plot classification data reference data vs sample
plt_pca_features <- function(data, legend = "bottomright", ...) {
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
  plot(data$ref_data$x, data$ref_data$y, type = "p",
       main = "Scatterplot", col = map_to_colors(data$ref_data$label),
       xlab = data$x_label, ylab = data$y_label, ...)
  points(data$pred_data$x, data$pred_data$y, type = "p", pch = 4,
         col = map_to_colors(data$pred_data$label)
         )
  if (legend != "none") {
    legend(legend, legend = c("Reference", "Sample",
                              "Persistent", "Contracting", "Late Emerging"),
           col = c("black", "black", "#8DC955", "#DE5A43", "#33528F"),
           pch = c(1, 4, 20, 20, 20))
  }
}

# find the closest reference data point and return its label
knn_predict <- function(reference_data, prediction_data, k = 1) {
  # find closest labels for each observation
  closest_labels <- sapply(seq_along(prediction_data$x), function(i) {
    distances <- sqrt((reference_data$x - prediction_data$x[i])^2 +
                      (reference_data$y - prediction_data$y[i])^2)
    distances <- order(distances, decreasing = FALSE)

    # Get the labels of the k-nearest neighbors
    k_nearest_labels <- reference_data$label[distances[1:k]]

    # Perform majority vote with first occurrence tie-breaking
    majority_label <- k_nearest_labels[max(table(k_nearest_labels))]

    # Return the predicted label
    return(majority_label)
  })

  # Add clone ids
  names(closest_labels) <- prediction_data$id

  return(closest_labels)
}

# Evaluate classification
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
