# Visualization script for skewed normal distribution sampling
# This script demonstrates the sample_skewed_normal function with different parameter settings

# clean environment
rm(list = ls())

# set seed for reproducibility
set.seed(42)

# load required libraries
require(sn)
require(ggplot2)
require(gridExtra)
require(dplyr)

# source the TCR model functions
source("R/tcr_model.R")

# create plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

#' Generate samples and create visualization plot
#'
#' This function generates a specified number of samples from the skewed normal
#' distribution and creates a histogram with density overlay.
#'
#' @param mean Mean parameter for the distribution
#' @param sd Standard deviation parameter
#' @param shape Shape parameter (skewness)
#' @param n_samples Number of samples to generate
#' @param title_suffix Additional text for plot title
#' @return A ggplot2 object
create_distribution_plot <- function(mean, sd, shape, n_samples = 10000, title_suffix = "") {
  # generate samples
  samples <- replicate(n_samples, sample_skewed_normal(mean, sd, shape))
  
  # calculate statistics
  sample_mean <- mean(samples)
  sample_sd <- sd(samples)
  
  # create data frame
  data <- data.frame(x = samples)
  
  # create plot
  p <- ggplot(data, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.7, 
                   fill = "lightblue", color = "black") +
    geom_density(color = "red", linewidth = 1) +
    labs(
      title = paste0("Skewed Normal Distribution", title_suffix),
      subtitle = paste0("Parameters: mean=", round(mean, 4), 
                       ", sd=", round(sd, 4), 
                       ", shape=", shape),
      x = "Value", 
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = paste0("Sample mean: ", round(sample_mean, 4), "\n",
                           "Sample sd: ", round(sample_sd, 4), "\n"),
             size = 3, color = "darkblue")
  
  return(p)
}

# define parameter settings to test
parameter_settings <- list(
  list(mean = 0.008, sd = 0.0008, shape = 0, title = "\n(Normal - Shape = 0)"),
  list(mean = 0.008, sd = 0.0008, shape = 2, title = "\n(Right Skewed - Shape = 2)"),
  list(mean = 0.008, sd = 0.0008, shape = -2, title = "\n(Left Skewed - Shape = -2)"),
  list(mean = 0.008, sd = 0.0008, shape = 5, title = "\n(Heavy Right Skew - Shape = 5)"),
  list(mean = 0.008, sd = 0.0016, shape = 0, title = "\n(Normal, Higher SD)"),
  list(mean = 0.008, sd = 0.0004, shape = 0, title = "\n(Normal, Lower SD)")
)

# create individual plots
plots <- list()
for (i in seq_along(parameter_settings)) {
  params <- parameter_settings[[i]]
  plots[[i]] <- create_distribution_plot(
    mean = params$mean,
    sd = params$sd, 
    shape = params$shape,
    title_suffix = params$title
  )
}

# save to PDF
pdf("plots/skewed_normal_visualization.pdf", width = 12, height = 15)
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

# also save as PNG for easier viewing
png("plots/skewed_normal_visualization.png", width = 1200, height = 1500, res = 100)
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

# print summary message
cat("Visualization complete!\n")
cat("Files saved:\n")
cat("- plots/skewed_normal_visualization.pdf\n")
cat("- plots/skewed_normal_visualization.png\n")
cat("\nParameter settings tested:\n")
for (i in seq_along(parameter_settings)) {
  params <- parameter_settings[[i]]
  cat(sprintf("  %d. mean=%.4f, sd=%.4f, shape=%d%s\n", 
              i, params$mean, params$sd, params$shape, params$title))
}