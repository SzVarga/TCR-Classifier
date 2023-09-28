# Function to return unique id
get_id <- function() {
  if (is.null(.GlobalEnv$id)) {
    .GlobalEnv$id <- 0
  }

  id  <<- id + 1

  return(id)
}

# Function to reset unique id
reset_id_counter <- function() {
  .GlobalEnv$id <- 0
}

# Function to return a specific label
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

# Function to map sim_time name to value
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

# Constructor function for TCR struct
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

# Function to add a Clonotype to the TCR-repertoire
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

# Model simulation
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

# return sample labels
tcr_labels <- function(repertoire, smpl) {
  ids <- as.integer(names(smpl))

  for (i in seq_along(smpl)) {
    smpl[[i]] <- get_label(repertoire$clonotypes[[ids[i]]]$label)
  }

  return(smpl)
}

# Sampling
# read data from simulation data
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

# sample from simulation data
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

# rbind samples
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

# identify the smallest overall sample size
# if clone i is not included
sample_min_size <- function(smpl_table) {
  smpl_colmin <- numeric(ncol(smpl_table))

  for (i in seq_len(ncol(smpl_table))) {
    smpl_colmin[i] <- min(rowSums(smpl_table[, -i]))
  }

  smpl_min <- min(smpl_colmin)
  return(smpl_min)
}

# resample to smaller sample size
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
# plot temporal repertoire composition
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

# boxplot and simple clone evolution plot
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
