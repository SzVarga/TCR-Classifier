#' Create a new TCR (T-cell receptor) object with viral dynamics support.
#'
#' This function creates a new TCR object for prime-boost modeling that includes
#' viral parameters alongside the standard TCR repertoire components. It supports
#' viral replication intervals, allowing simulation of prime-boost vaccination
#' scenarios or chronic infections with periodic viral replication.
#'
#' @param sim_times A vector of simulation times for the TCR object.
#' @param carry_cap The carrying capacity of the TCR object.
#' @param viral_params A list containing viral parameters:
#'   - initial_load: Initial viral load (numeric)
#'   - replication_rate: Rate of viral replication when active (numeric)
#'   - clearance_rate: Immune-mediated viral clearance rate per T cell (numeric)
#'   - replication_intervals: List of [start, end] time intervals when virus replicates
#'   - viral_burden_sensitivity: System-wide sensitivity to cumulative viral burden (numeric, default: 0)
#' @return A new TCR object with viral dynamics support.
#'
#' @examples
#' # Prime-boost vaccination schedule
#' viral_params <- list(
#'   initial_load = 1e6,
#'   replication_rate = 0.8,
#'   clearance_rate = 0.2,
#'   replication_intervals = list(c(0, 14), c(180, 194)),
#'   viral_burden_sensitivity = 0.1  # System-wide memory effect
#' )
#' tcr <- new_tcr_primeBoostModel(sim_times, carry_cap, viral_params)
#'
#' @export
new_tcr_primeBoostModel <- function(sim_times, carry_cap, viral_params) {
  # sim_times must have length 1 or more
  stopifnot(length(sim_times) >= 1)
  
  # Basic validation of viral parameters
  stopifnot(!is.null(viral_params$initial_load))
  stopifnot(!is.null(viral_params$replication_rate))
  stopifnot(!is.null(viral_params$clearance_rate))
  stopifnot(!is.null(viral_params$replication_intervals))
  stopifnot(is.list(viral_params$replication_intervals))
  
  # Validate that initial_load is positive
  stopifnot(viral_params$initial_load >= 0)
  stopifnot(viral_params$replication_rate >= 0)
  stopifnot(viral_params$clearance_rate >= 0)
  
  # Set default viral burden sensitivity if not provided
  if (is.null(viral_params$viral_burden_sensitivity)) {
    viral_params$viral_burden_sensitivity <- 0
  }
  stopifnot(is.numeric(viral_params$viral_burden_sensitivity))
  
  # reset id numbers
  reset_id_counter()

  # construct TCR struct with viral parameters
  tcr <- list(
    clonotypes = list(),
    sim_times = sim_times,
    carry_cap = carry_cap,
    clone_labels = c(),
    data = matrix(),
    viral_params = viral_params
  )
  return(tcr)
}

#' Check if viral replication is active at a given time point.
#'
#' This indicator function determines whether viral replication should occur
#' at the specified time point based on the provided replication intervals.
#' Returns 1 if the time falls within any replication interval, 0 otherwise.
#'
#' @param current_time The current simulation time point
#' @param replication_intervals List of [start, end] time intervals for viral replication
#' @return 1 if viral replication is active, 0 otherwise
#'
#' @examples
#' # Define replication intervals for prime-boost
#' intervals <- list(c(0, 10), c(150, 170))
#' 
#' # Check if virus replicates at different time points
#' is_viral_replication_active(5, intervals)   # Returns 1 (within first interval)
#' is_viral_replication_active(50, intervals)  # Returns 0 (between intervals)
#' is_viral_replication_active(160, intervals) # Returns 1 (within second interval)
#'
#' @export
is_viral_replication_active <- function(current_time, replication_intervals) {
  # explicitly handle empty list intervals
  if (length(replication_intervals) == 0) {
    return(0)  # No intervals means virus never replicates
  }

  for (interval in replication_intervals) {
    if (current_time >= interval[1] && current_time <= interval[2]) {
      return(1)  # Virus replicates
    }
  }
  return(0)  # Virus does not replicate
}

#' Calculate cumulative viral burden from simulation results.
#'
#' This helper function calculates the cumulative viral load by summing
#' all viral particles seen at previous time points.
#'
#' @param results_data Matrix containing simulation results with virus column
#' @param current_time Current simulation time point
#' @return Sum of all viral particles from previous time points
#'
calculate_cumulative_viral_burden <- function(results_data, current_time) {
  if (is.null(results_data) || length(results_data) == 0 || 
      !is.matrix(results_data) || nrow(results_data) == 0) {
    return(0)  # No previous data or empty matrix
  }
  
  # Check if required columns exist
  if (is.null(colnames(results_data)) || 
      !("time" %in% colnames(results_data)) || 
      !("virus" %in% colnames(results_data))) {
    return(0)  # Missing required columns
  }
  
  # Sum all viral loads from previous time points
  previous_time_filter <- results_data[, "time"] < current_time
  if (!any(previous_time_filter)) {
    return(0)  # No previous time points
  }
  
  previous_viral_loads <- results_data[previous_time_filter, "virus"]
  return(sum(previous_viral_loads))
}

#' Simulate stochastic clonal dynamics with viral load
#' of the TCR (T-cell receptor) over a time partition.
#'
#' This function simulates T-cell clonal expansion and viral dynamics over time 
#' using the Adaptive Tau-Leaping algorithm. It extends the basic model to include
#' viral replication, clearance, and avidity-based T-cell proliferation responses.
#'
#' @param repertoire The TCR repertoire object containing clones and viral parameters
#' @param init_values A list of initial values for each clone and virus
#' @param param_idx The index of simulation parameters to use
#' @param global_time The current global simulation time (for viral intervals)
#' @return Simulation data including clone populations and viral load over time
#'
tcr_simulate_tpart_primeBoostModel <- function(repertoire, init_values, param_idx, global_time = 0) {
    # generate required data for simulation
    transitions <- list()
    params <- list()

    # iterate through clonotypes to set up clone transitions
    for (clonotype in repertoire$clonotypes) {
      # generate transitions and params lists for birth and death
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
        
        # Add avidity parameters for birth rates
        if (pr == "birth" && !is.null(clonotype$params$avidity)) {
          avidity_param <- c(`names<-`(clonotype$params$avidity,
                                paste("avidity", clonotype$clone_id, sep = ".")))
          params <- c(params, avidity_param)
        }
        
      }
    }
    
    # Add viral transitions: replication (+1) and clearance (-1)
    transitions <- c(transitions, list(c(virus = 1)), list(c(virus = -1)))
    
    # Add viral parameters
    viral_replication_param <- c(`names<-`(repertoire$viral_params$replication_rate, "viral_replication"))
    viral_clearance_param <- c(`names<-`(repertoire$viral_params$clearance_rate, "viral_clearance"))
    params <- c(params, viral_replication_param, viral_clearance_param)
    
    # Store viral replication intervals and global time for rate function
    params$viral_replication_intervals <- repertoire$viral_params$replication_intervals
    params$global_time <- global_time
    params$viral_burden_sensitivity <- repertoire$viral_params$viral_burden_sensitivity
    params$repertoire_data <- repertoire$data

    # fetch TCR-repertoire carrying capacity
    carry_cap <- repertoire$carry_cap

    # Enhanced rate function with viral dynamics
    rate_func <- function(vars, params, t) {
      # Extract populations (last element is virus)
      clone_vars <- vars[1:(length(vars)-1)]
      viral_load <- vars[length(vars)]
      
      # Calculate total clone population
      pop_all <- sum(clone_vars)
      
      # Calculate cumulative viral burden for this time point
      current_absolute_time <- params$global_time + t
      cumulative_burden <- calculate_cumulative_viral_burden(params$repertoire_data, current_absolute_time)

      # Calculate rate transitions for clones
      vec <- c()
      for (i in seq_along(clone_vars)) {
        birth_param <- paste("birth", i, sep = ".")
        death_param <- paste("death", i, sep = ".")
        avidity_param <- paste("avidity", i, sep = ".")
        
        # Base birth rate with logistic growth and viral burden effect
        base_birth_rate <- params[[birth_param]] * clone_vars[[i]] * (1 - pop_all / carry_cap)
        
        # Modify base birth rate based on cumulative viral burden
        if (params$viral_burden_sensitivity > 0 && cumulative_burden > 0) {
          viral_burden_effect <- params$viral_burden_sensitivity * log10(cumulative_burden + 1)
          base_birth_rate <- base_birth_rate * (1 + viral_burden_effect)
        }
        
        # Apply avidity-based proliferation boost if parameter exists
        if (avidity_param %in% names(params) && viral_load > 0) {
          avidity_boost <- params[[avidity_param]] * viral_load
          birth_rate <- max(0, base_birth_rate + avidity_boost)
        } else {
          birth_rate <- max(0, base_birth_rate)
        }
        
        # Death rate
        death_rate <- params[[death_param]] * clone_vars[[i]]
        
        vec <- c(vec, birth_rate, death_rate)
      }
      
      # Calculate viral dynamics
      current_absolute_time <- params$global_time + t
      replication_active <- is_viral_replication_active(
        current_absolute_time, 
        params$viral_replication_intervals
      )
      
      # Viral replication (only when active)
      viral_replication_rate <- replication_active * params$viral_replication * viral_load
      
      # Viral clearance (immune-mediated only, based on total population)
      viral_clearance <- params$viral_clearance * pop_all * viral_load
      
      # Add viral rates to the vector
      vec <- c(vec, viral_replication_rate, viral_clearance)

      # return rate transitions
      return(vec)
    }

    # perform the simulation
    data <- adaptivetau::ssa.adaptivetau(init_values, transitions,
                                         rate_func, params,
                                         repertoire$sim_times[[param_idx]])

    # calculate total clone population size (excluding virus)
    total_clones <- c()
    for (i in seq_len(nrow(data))) {
      total_clones <- c(total_clones, sum(data[i, 2:(ncol(data)-1)]))
    }

    # collect data with virus and total_clones columns
    data <- cbind(data, total_clones)

    # return data
    return(data)
}

#' Simulate T-cell clonal dynamics with viral load over multiple time partitions.
#'
#' This function simulates T-cell clonal expansion and viral dynamics over
#' multiple time partitions using the provided TCR repertoire and viral
#' parameters. It tracks global time to properly handle viral replication
#' intervals that may span multiple time partitions.
#'
#' @param repertoire The TCR repertoire object containing clones and viral
#'   parameters
#' @param ... Optional arguments. Currently supports the following argument:
#'   - benchmark : If set to TRUE, calculation time will be displayed.
#' @return Simulation data including clone populations and viral load over time
#'
#' @export
tcr_simulate_prime_boost_model <- function(repertoire, ...) {
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
      # add initial viral load
      init_values <- c(init_values,
                       `names<-`(repertoire$viral_params$initial_load, "virus"))
    } else {
      #clear init.values
      init_values <- c()
      # fetch from results (all columns except time and total_clones)
      for (col in 2:(ncol(data) - 1)) {
        init_values <- c(init_values,
                         `names<-`(results[nrow(results), col],
                                   colnames(results)[col]))
      }
    }

    # simulate time partition with current global time
    data <- tcr_simulate_tpart_primeBoostModel(repertoire, init_values, tPart,
                                               glob_time)

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

#' Add a new clone to the TCR (T-cell receptor) object for prime-boost modeling.
#'
#' This function adds a new clone with avidity parameters to an existing
#' TCR object for prime-boost modeling. Clones with higher avidity show
#' stronger proliferative responses to viral load. Viral clearance is handled 
#' at the population level.
#'
#' @param tcr The TCR object to which the clonotype will be added.
#' @param label The label of the clonotype.
#' @param init_size The initial size of the clonotype.
#' @param birth A vector of birth rates corresponding to simulation times.
#' @param death A vector of death rates corresponding to simulation times.
#' @param avidity TCR avidity parameter determining proliferation response to
#'   viral load (default: 0). Higher values indicate stronger TCR-antigen
#'   binding and greater proliferative response.
#' @return The updated TCR object with the new clonotype.
#'
#' @examples
#' # Create a new TCR object with viral parameters
#' viral_params <- list(
#'   initial_load = 1e6,
#'   replication_rate = 0.8,
#'   clearance_rate = 0.2,
#'   replication_intervals = list(c(0, 14))
#' )
#' tcr <- new_tcr_primeBoostModel(sim_times, carry_cap, viral_params)
#'
#' # Add a clone with high avidity (strong response to viral load)
#' tcr <- add_clone_prime_boost_model(tcr = tcr, label = "persistent",
#'                                    init_size = 100,
#'                                    birth = rep(0.008, 9),
#'                                    death = rep(0.0032, 9),
#'                                    avidity = 0.5)
#'
#' @export
add_clone_prime_boost_model <- function(tcr, label, init_size, birth, death,
                                        avidity = 0) {
  # birth-death parameters must have same length as sim_times
  stopifnot(length(birth) == length(tcr$sim_times))
  stopifnot(length(death) == length(tcr$sim_times))

  # construct clonotype struct with viral parameters
  clonotype <- list(
    clone_id = get_id(),
    label = get_label(label),
    init_size = init_size,
    params = list(
      birth = birth,
      death = death,
      avidity = avidity
    )
  )

  # append clonotype to tcr
  tcr$clonotypes <- c(tcr$clonotypes, list(clonotype))
  # append label to tcr
  tcr$clone_labels <- unique(c(tcr$clone_labels, label))

  return(tcr)
}