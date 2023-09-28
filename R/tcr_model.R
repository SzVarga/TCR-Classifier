#' Get and increment the unique identifier for clones in the TCR object.
#'
#' This function retrieves a unique identifier (ID) from the global environment.
#' If the ID does not exist, it initializes it to 0 and then increments it by 1.
#'
#' @return The unique identifier.
#'
get_id <- function() {
  if (is.null(.GlobalEnv$id)) {
    .GlobalEnv$id <- 0
  }

  id  <<- id + 1

  return(id)
}

#' Reset the unique identifier counter to zero.
#'
#' This function sets the unique identifier (ID) in the global environment to
#' zero.
#'
reset_id_counter <- function() {
  .GlobalEnv$id <- 0
}

#' Get label from index or name.
#'
#' This function retrieves a label from either an index (1, 2, or 3) or a name
#' ("persistent", "contracting", or "late_emerging"). It maps indices to names
#' and vice versa based on a predefined mapping.
#'
#' @param x The index (1, 2, 3) or name ("persistent", "contracting", "late_emerging").
#' @return The corresponding label (factor).
#' @throws An error is thrown if the input is not valid.
#'
#' @examples
#' # Get the label for an index
#' get_label(2)
#'
#' # Get the label for a name
#' get_label("persistent")
#'
#' @export
get_label <- function(x) {
  # convert label to index
  if (x == 1 || x == 2 || x == 3) {
    # nothing to do
  } else if (x == "persistent") {
    x <- 1
  } else if (x == "contracting") {
    x <- 2
  }else if (x == "late_emerging") {
    x <- 3
  } else {
    stop("x must be one of the following \n 1 or \"persistent\", 
    \n 2 or \"contracting\", \n 3 or \"late_emerging\"")
  }

  # check if global label (factor) exists
  if (is.null(.GlobalEnv$Label)) {
    .GlobalEnv$Label <- factor(c("persistent", "contracting", "late_emerging"),
    levels = c("persistent", "contracting", "late_emerging"))
  }

  return(.GlobalEnv$Label[x])
}

#' Map simulation event name to event time.
#'
#' This function maps a simulation event name to its corresponding event time.
#' This function allows for the use of event names in sampling functions.
#'
#' @param simt_name The name of the simulation event.
#' @return The event time corresponding to the event name.
#'
#' @examples
#' # Map an event name to its time
#' sim_name2time("P1")
#'
#' @export
sim_name2time <- function(simt_name) {
  # map evnt name to evnt time
  if (simt_name == "begin") {
    time_point <- 0
  } else if (simt_name == "P1") {
    time_point <- 10
  } else if (simt_name == "S1") {
    time_point <- 20
  } else if (simt_name == "S2") {
    time_point <- 30
  } else {
    # unknown event name
    time_point <- simt_name
  }

  return(time_point)
}

#' Create a new TCR (T-cell receptor) object.
#'
#' This function creates a new TCR object with the specified simulation times,
#' carry capacity and the simulation data. This object also keeps track of the
#' unique clone labels within the TCR.
#'
#' @param sim_times A vector of simulation times for the TCR object.
#' @param carry_cap The carrying capacity of the TCR object.
#' @return A new TCR object.
#'
#' @examples
#' # Create a new TCR object with simulation times and carry capacity
#' tcr <- new_tcr(sim_times = c("P1" = 10, "S1" = 10, "S2" = 10),
#'                carry_cap = num_clones * clone_size / 0.6)
#'
#' @export
new_tcr <- function(sim_times, carry_cap) {
  # sim_times must have length 1 or more
  stopifnot(length(sim_times) >= 1)

  # reset id numbers
  reset_id_counter()

  # construct TCR struct
  tcr <- list(
    clonotypes = list(),
    sim_times = sim_times,
    carry_cap = carry_cap,
    clone_labels = c(),
    data = matrix()
  )
  return(tcr)
}

#' Add a new clone to the TCR (T-cell receptor) object.
#'
#' This function adds a new clone with specified characteristics to an
#' existing TCR object. Clones represent individual T-cell clones with
#' attributes such as label, initial size, birth, and death rates.
#'
#' @param tcr The TCR object to which the clonotype will be added.
#' @param label The label of the clonotype.
#' @param init_size The initial size of the clonotype.
#' @param birth A vector of birth rates corresponding to simulation times.
#' @param death A vector of death rates corresponding to simulation times.
#' @return The updated TCR object with the new clonotype.
#'
#' @examples
#' # Create a new TCR object
#' tcr <- new_tcr(c(0, 10, 20, 30), 100)
#'
#' # Add a new clonotype to the TCR object
#' tcr <- add_clone(tcr = tcr, label = "persistent",
#'                  init_size = 100,
#'                  birth = 100 * c(0.008, 0.008, 0.008),
#'                  death = 100 * c(0.0032, 0.0032, 0.0032))
#'
#' @export
add_clone <- function(tcr, label, init_size, birth, death) {
  # birth-death parameters must have same length as sim_times
  stopifnot(length(birth) == length(tcr$sim_times))
  stopifnot(length(death) == length(tcr$sim_times))

  # construct clonotype struct
  clonotype <- list(
    clone_id = get_id(),
    label = get_label(label),
    init_size = init_size,
    params = list(birth = birth, death = death)
  )

  # append clonotype to tcr
  tcr$clonotypes <- c(tcr$clonotypes, list(clonotype))
  # append label to tcr
  tcr$clone_labels <- unique(c(tcr$clone_labels, label))

  return(tcr)
}

#' Simulate stochastic clonal dynamics
#' of the TCR (T-cell receptor) over a time partition.
#'
#' This function simulates T-cell clonal expansion over time using the Adaptive
#' Tau-Leaping algorithm. It takes a TCR repertoire, initial values, and a
#' parameter index as input and returns simulation data.
#'
#' @param repertoire The TCR repertoire object containing clones
#' and simulation parameters.
#' @param init_values A list of initial values for each clone.
#' @param param_idx The index of simulation parameters to use.
#' @return Simulation data including clone populations over time.
#'
tcr_simulate_tpart <- function(repertoire, init_values, param_idx) {
    # generate required data for simulation
    init_values <- init_values
    transitions <- list()
    params <- list()

    # iterate through clonotypes
    for (clonotype in repertoire$clonotypes) {

      # generate transitions and params lists
      for (i in  1:2) {
        # transitions
        tr <- c(1, -1)[i]
        transition <- c(`names<-`(tr, clonotype$clone_id))
        transitions <- c(transitions, list(transition))

        # params
        pr <- c("birth", "death")[i]
        param <- c(`names<-`(clonotype$params[[pr]][param_idx],
                   paste(pr, clonotype$clone_id, sep = ".")))
        params <- c(params, param)
      }
    }

    # fetch TCR-repertoire carrying capacity
    carry_cap <- repertoire$carry_cap

    # function to return rate transitions
    rate_func <- function(vars, params, t) {

      # calculate total population
      pop_all <- Reduce(sum, vars)

      # calculate rate transitions
      vec <- c()
      for (i in seq_along(vars)) {
        birth <- paste("birth", i, sep = ".")
        death <- paste("death", i, sep = ".")
        # birth_rate of i-th clone
        vec <- c(vec, max(0, params[[birth]] *
                 vars[[i]] * (1 - pop_all / carry_cap)))
        # death rate of i-th clone
        vec <- c(vec, params[[death]] * vars[[i]])
      }

      # return rate transitions
      return(vec)
    }

    # perform the simulation
    data <- adaptivetau::ssa.adaptivetau(init_values, transitions,
                                         rate_func, params,
                                         repertoire$sim_times[[param_idx]])

    # calculate total population size
    total <- c()
    for (i in seq_len(nrow(data))) {
      total <- c(total, sum(data[i, 2:ncol(data)]))
    }

    # collect data
    data <- cbind(data, total)

    # return data
    return(data)
}

#' Simulate T-cell clonal dynamics over multiple time partitions.
#'
#' This function simulates T-cell clonal expansion over multiple time partitions
#' using the provided TCR repertoire and optional arguments. It iterates through
#' each time partition, calls the funtion that simulates the dynamics over the
#' time partition and combines the results to create a global simulation.
#'
#' @param repertoire The TCR repertoire object containing
#' clones and simulation parameters.
#' @param ... Optional arguments. Currently supports the following argument:
#' - benchmark : If set to TRUE, calculation time will be displayed.
#' @return Simulation data including clone populations over time.
#'
#' @export
tcr_simulate <- function(repertoire, ...) {
  # optional arguments
  args <- list(...)

  # start optional benchmarking
  if (!is.null(args$benchmark)) {
    starttime <- Sys.time()
  }

  # variables to track global simulation
  glob_time <- 0
  init_values <- c()
  results <- NULL

  # iterate over each simulation partition
  for (tPart in seq_along(repertoire$sim_times)){
    # fetch initial values
    if (is.null(results)) {
      # fetch from clonotypes
      for (clonotype in repertoire$clonotypes) {
        init_values <- c(init_values,
                         `names<-`(clonotype$init_size, clonotype$clone_id))
      }
    } else {
      #clear init.values
      init_values <- c()
      # fetch from results (columns 2 to second to last)
      for (col in 2:(ncol(data) - 1)) {
        init_values <- c(init_values,
                         `names<-`(results[nrow(results), col],
                                   colnames(results)[col]))
      }
    }

    # simulate time partition
    data <- tcr_simulate_tpart(repertoire, init_values, tPart)

    # combine results
    if (is.null(results)) {
      # override
      results <- data
    } else {
      # update data time
      data[, "time"] <- data[, "time"] + glob_time

      # append without duplicate
      results <- rbind(results, data[-1, ])
    }
    # update global simulation time
    glob_time <- glob_time + repertoire$sim_times[[tPart]]
  }

  # end optional benchmarking
  if (!is.null(args$benchmark)) {
    endtime <- Sys.time()
    cat("simulation time:",
        difftime(endtime, starttime, units = "secs"), "s", "\n")
  }

  # return results
  return(results)
}

#' Extract T-cell receptor (TCR) labels from a sample.
#'
#' This function extracts TCR labels from a sample based on the provided
#' TCR repertoire. It takes a sample (list or vector) of clone IDs and returns
#' a vector of corresponding TCR labels.
#'
#' @param repertoire The TCR repertoire object containing the clones and labels.
#' @param smpl A sample of clone IDs for which labels will be extracted.
#' @return A vector of TCR labels corresponding to the sample.
#'
tcr_labels <- function(repertoire, smpl) {
  ids <- as.integer(names(smpl))

  for (i in seq_along(smpl)) {
    smpl[[i]] <- get_label(repertoire$clonotypes[[ids[i]]]$label)
  }

  return(smpl)
}

# Sampling
#' Retrieve true clonal abundance data at a specific time point.
#'
#' This function retrieves clonal data from a repertoire at a specific
#' time point. You can specify additional arguments to filter the data based
#' on clone ID or other criteria.
#'
#' @param repertoire The TCR repertoire object containing clonotype data.
#' @param time_point The time point at which to retrieve data.
#' @param ... Additional optional arguments.
#' Currently supports the following argument:
#' -clone_id : Specify a clone ID to filter data for a specific clone or clones.
#' @return Clonotype data at the specified time point,
#' optionally filtered by clone ID.
#'
#' @export
data_at <- function(repertoire, time_point, ...) {
  # additional arguments
  args <- list(...)

  # convert event name to timePoint
  time_point <- sim_name2time(time_point)

  # return clonotype-data at given time_point
  filter <- repertoire$data[, "time"] == time_point
  if (!any(filter)) {
    stop("Time point has not been calculated!")
  }

  if ("clone_id" %in% names(args)) {
    # return only clones matching clone_id
    return(repertoire$data[filter, as.character(args$clone_id)])

  } else {
    # return all clone data
    return(repertoire$data[filter, ])
  }
}

#' Sample T-cell receptor (TCR) clones at a specific time point.
#'
#' This function samples TCR clones from a repertoire at a specific time point,
#' applying a detection limit and optionally specifying the sample size. It
#' returns a table with the sampled clone IDs and their counts.
#'
#' @param repertoire The TCR repertoire object containing clonotype data.
#' @param size The number of clones to sample. Defaults to 15.
#' Use "all" to sample all available clones.
#' @param time_point The time point at which to sample clones.
#' Supports event names such as "P1", "S1", "S2".
#' @param detect_lim The detection limit for clone counts.
#' Clones with counts less than this limit will be excluded.
#' @return A table with sampled clone IDs and their counts.
#'
#' @export
sample_at <- function(repertoire, size = 15, time_point, detect_lim = 3) {
  # fetch tcr data at time_point
  data <- data_at(repertoire = repertoire, time_point = time_point)

  # apply detection limit
  filter <- data[2:(length(data) - 1)] >= detect_lim
  filter <- c(TRUE, filter, TRUE)
  data <- data[filter]

  # check if pool size > size
  pool_size <- Reduce(sum, data[2:(length(data) - 1)])
  if (size == "all") {
    size <- pool_size
  }
  if (pool_size <= size) {
    return(data[2:(length(data) - 1)])
  }
  # create pool to sample from
  ids <- names(data)[2:(length(data) - 1)]

  # sampled tcr
  smpld <- numeric(size)
  # sampling sequentially
  for (i in 1:size){
    # calculate weights
    weights <- data[2:(length(data) - 1)] / data[length(data)]
    # sample clone id
    smpl <- sample(ids, size = 1, replace = FALSE, prob = weights)
    # append sampled clone to sampled vector
    smpld[i] <- smpl
    # decrement sampled clone
    data[[smpl]] <- data[[smpl]] - 1
  }
  return(table(smpld))
}

#' Bind vectors into a matrix by matching names.
#'
#' This function takes a variable number of input vectors and binds them into a
#' matrix by matching names. It returns a matrix where each row represents an
#' input vector, and columns are the union of all unique names from the input
#' vectors. Missing values are filled with zeros.
#'
#' @param ... Input vectors to be bound into a matrix.
#' @return A matrix with input vectors matched by names.
#'
#' @export
sample_bind <- function(...) {
  # Convert input vectors to a list
  vector_list <- list(...)

  # Get the union of all unique names
  all_names <- unique(unlist(lapply(vector_list, names)))
  vector_names <- names(list(...))

  # Initialize an empty matrix to store the results
  result_matrix <- matrix(0, nrow = length(vector_list),
                          ncol = length(all_names))

  # Fill the matrix with values from input vectors
  for (i in seq_along(vector_list)) {
    names_i <- names(vector_list[[i]])
    values_i <- as.vector(vector_list[[i]])
    col_indices <- match(names_i, all_names)
    result_matrix[i, col_indices] <- values_i
  }

  # name matrix columns & rows
  colnames(result_matrix) <- all_names
  rownames(result_matrix) <- vector_names

  return(result_matrix)
}

#' Calculate the minimum total count across all clones in a sample table.
#'
#' This function calculates the minimum total count across all clones in a
#' sample table if a single clone is not counted. The sample table is typically
#' the result of the `sample_at` function and represents sampled clone counts.
#' This functuin is used to determine the minimum sample size to wich all sample
#' are resampled to avoid sampling bias.
#'
#' @param smpl_table A sample table with clone counts.
#' @return The minimum total count across all clones in the sample table.
#'
#' @references{
#'  \insertRef{Venturi2007}{title = Methods for comparing the diversity
#'   of samples of the T cell receptor repertoire},
#'  \insertRef{Venturi2008}{title = Method for assessing the similarity between
#'   subsets of the T cell receptor repertoire
#' }
#'
#' @export
sample_min_size <- function(smpl_table) {
  smpl_colmin <- numeric(ncol(smpl_table))

  for (i in seq_len(ncol(smpl_table))) {
    smpl_colmin[i] <- min(rowSums(smpl_table[, -i]))
  }

  smpl_min <- min(smpl_colmin)
  return(smpl_min)
}

#' Reduce the counts of a specific clone in a sample table by resampling.
#'
#' This function reduces the counts of cells in a sample table by
#' resampling, optionally including or excluding a specific clone.
#'
#' @param smpl_table A sample table with clone counts.
#' @param size The number of cells to resample to.
#' @param clone The name of the specific clone.
#' @param incl Logical. If TRUE, the target clone is included in resampling.
#'             If FALSE, the target clone is excluded.
#' @return A modified sample table with reduced and resampled counts.
#'
#' @references{
#'  \insertRef{Venturi2007}{title = Methods for comparing the diversity
#'   of samples of the T cell receptor repertoire},
#'  \insertRef{Venturi2008}{title = Method for assessing the similarity between
#'   subsets of the T cell receptor repertoire
#'
#' @export
sample_reduce <- function(smpl_table, size, clone, incl) {
  # create empty copy of smpl_table
  re_sampled <- matrix(0, nrow = nrow(smpl_table), ncol = ncol(smpl_table))
  colnames(re_sampled) <- colnames(smpl_table)
  rownames(re_sampled) <- rownames(smpl_table)

  # clone index in smpl_table
  clo_idx <- which(colnames(smpl_table) == clone)

  # for each time point
  for (tp in seq_along(rownames(smpl_table))) {
    reps <- size
    # decide on excluding / including the clone
    if (smpl_table[tp, clo_idx] > 0 && incl == TRUE) {
      # add clone once to re.sampled manually
      re_sampled[tp, clo_idx] <- re_sampled[tp, clo_idx] + 1
      smpl_table[tp, clo_idx] <- smpl_table[tp, clo_idx] - 1
      reps <- reps - 1
    } else if (incl == FALSE) {
      # exclude clone from sampling
      smpl_table[tp, clo_idx] <- 0
    }

    # sample from smpl_table reps times
    # draw from pool
    pool <- rep(names(smpl_table[tp, ]), times = smpl_table[tp, ])
    smpld <- sample(pool, reps, replace = FALSE)
    smpld_names <- names(table(smpld))
    smpld_values <- table(smpld)
    # fill re.sampled vector
    re_sampled[tp, smpld_names] <- re_sampled[tp, smpld_names] + smpld_values
  }

  return(re_sampled)
}

# Plotting
#' Create a temporal composition plot for a T-cell receptor (TCR) repertoire.
#'
#' This function generates a temporal composition plot, showing clonal dynamics
#' over time. It uses the ggplot2 package for visualization.
#'
#' @param data A data frame or matrix with temporal TCR repertoire data.
#' Rows represent time points, and columns represent clones.
#' The last column might include the total cell count.
#' @param ... Additional optional arguments.
#' Currently supports the following arguments:
#' - carry_cap : Add a horizontal dashed line
#' to represent the carrying capacity of the TCR.
#' - legend_labels : Specify custom legend labels for clones.
#' Provide a vector of labels.
#' 
#' @return A ggplot2 plot object displaying
#' the temporal TCR repertoire composition.
#'
#' @export
tcr_temporal_comp <- function(data, ...) {
  # optional parameter
  args <- list(...)

  # convert data into a long format
  data <- tidyr::pivot_longer(
    data = as.data.frame(data),
    cols = 2:ncol(data),
    names_to = "clone",
    values_to = "size"
  )

  # convert clone id to factor
  data$clone <- as.factor(data$clone)

  # Plot using ggplot2
  plt <- ggplot(data, aes(x = time, y = size, color = clone)) +
    geom_line(linewidth = 1) +
    labs(x = "Time", y = "Size", color = "Clone",
         title = "Temporal TCR Repertoire Composition") +
    scale_x_continuous("Time", breaks = c(0, 10, 20, 30),
                       labels = c(0, "P1", "S1", "S2")) +
    theme_bw()

  # optional features
  if ("carry_cap" %in% names(args)) {
    plt <- plt + geom_hline(yintercept = args$carry_cap,
                            linetype = "dashed", color = "black")
  }
  if ("legend_labels" %in% names(args)) {
    plt <- plt + scale_color_discrete(labels = c(args$legend_labels, "total"))
  }

  return(plt)
}

#' Create a temporal box and evolution plot for TCR data.
#'
#' This function generates a plot combining a boxplot of clone size distribution
#' and a clone size evolution plot. The boxplot provides a snapshot of the
#' clone size distribution at different time points, while the evolution plot
#' shows how individual clone sizes change over time.
#'
#' @param repertoire A TCR repertoire object containing clonotype and time data.
#' @param init_smpl A named vector representing the initial clone sizes.
#' @param box Logical. If TRUE, display the boxplot of clone size distribution.
#' @param evo Logical. If TRUE, display the evolution plot of clone size changes
#' @param rel Logical. If TRUE, plot clone sizes relative to first time point;
#' if FALSE, plot them as absolute counts.
#' @param ... Additional optional arguments.
#' Currently supports the following arguments:
#' - main : The title of the plot. Default is "Temporal TCR Distribution."
#' - box_col : The color of the boxplot. Default is white.
#' @return A combined plot displaying the temporal TCR repertoire composition
#' and evolution.
#'
#' @export
tcr_box_evo <- function(repertoire, init_smpl, box = TRUE,
                        evo = TRUE, rel = TRUE, ...) {
  # aditional arguments
  args <- list(...)

  # Sample data
  init_ids <- as.integer(names(init_smpl))
  init_size <- unname(init_smpl)

  # combine the sampled data into a list
  data_list <- list()

  # calculate absolute sample sizes and set y_label accordingly
  if (rel) {
    text_yax <- "Relative Clone Size"
  } else {
    text_yax <- "Absolute Clone Size"
    init_size <- 1
  }

  # retrieve clone data of sample at diff timepoints
  tps <- c("begin", names(repertoire$sim_times))
  for (tp in tps) {
    data_list <- append(data_list,
                        list(data_at(repertoire, tp, clone_id = init_ids) / init_size)) # nolint: line_length_linter.
  }

  # show boxplot?
  if (box) {
    data_list_box <- data_list
  } else {
    data_list_box <- vector("list", length(data_list))
  }

  # modify plot according to additional args
  if ("main" %in% names(args)) {
    main <- args$main
  } else {
    main <- "Temporal TCR Distribution"
  }
  if ("box_col" %in% names(args)) {
    boxcol <- args$boxcol
  } else {
    boxcol <- rgb(1, 1, 1, 1)
  }

  # colorize based on clone label
  col_labels <- c()
  for (id in init_ids) {
    label <- repertoire$clonotypes[[id]]$label
    if (label == "persistent") {
      col_labels <- c(col_labels, "#8DC955")
    } else if (label == "contracting") {
      col_labels <- c(col_labels, "#DE5A43")
    } else if (label == "late_emerging") {
      col_labels <- c(col_labels, "#33528F")
    } else {
      # unknown label
      col_labels <- c(col_labels, rgb(0, 0, 0, 0.5))
    }
  }

  # create the boxplots
  plt <- boxplot(data_list_box, main = main, ylab = text_yax, xlab = "Time",
                 col = boxcol, xaxt = "n", outline = FALSE,
                 ylim = c(min(unlist(data_list)), max(unlist(data_list))))

  # add x-axis labels
  axis(1, at = seq_along(data_list), labels = tps)

  # show evolution of clones?
  if (evo) {
    # connect sizes of clones (lines)
    for (m in seq_along(data_list[[1]])) {
      vec <- c()
      for (n in seq_along(data_list)) {
        vec <- c(vec, data_list[[n]][m])
      }
      lines(vec, col = col_labels[m])
    }

    # add clone sizes @timePoint (points)
    for (i in seq_along(data_list)) {
      points(rep(i, length(data_list[[i]])), data_list[[i]],
             pch = 1, col = col_labels)
    }
  }

  # Create a legend
  legend("topleft", legend = c("Persistent", "Contracting", "Late Emerging"),
         fill = c("#8DC955", "#DE5A43", "#33528F"), title = "Clone Label")

  return(plt)
}
