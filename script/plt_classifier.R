# clean env
rm(list = ls())

# load libraries
require(dplyr)
require(ggplot2)
require(ggsignif)

# set path to data
path_naive <- "data/classifier/full/naive"
path_pca_knn <- "data/classifier/full/pca_knn"
path_pca_mlr <- "data/classifier/full/pca_mlr"
path_mlr <- "data/classifier/full/mlr"
path_knn <- "data/classifier/full/knn"

# List all RDS files in a directory (change the path to your directory)
file_naive <- list.files(path_naive, "*.rds", full.names = TRUE)
file_pca_knn <- list.files(path_pca_knn, "*.rds", full.names = TRUE)
file_pca_mlr <- list.files(path_pca_mlr, "*.rds", full.names = TRUE)
file_mlr <- list.files(path_mlr, "*.rds", full.names = TRUE)
file_knn <- list.files(path_knn, "*.rds", full.names = TRUE)

# combine files
file_list <- c(file_naive, file_pca_knn, file_pca_mlr, file_mlr, file_knn)

# Initialize an empty data frame to store the combined data
data <- data.frame()

# Loop through each RDS file and append its data to combined_data
for (file_path in file_list) {
  # Load the data from the RDS file
  loaded_data <- readRDS(file_path)

  # Append the loaded data to combined_data
  data <- rbind(data, loaded_data)
}

# create list of classifier combinations for geom_signif
combinations <- combn(unique(data$classifier), 2)
# Convert columns of the matrix to a list of vectors
combinations <- lapply(seq_len(ncol(combinations)), function(col) {
                                                      combinations[, col]})

# set the plotting order of classifiers
cls_order <- c("naive", "knn", "pca_knn", "mlr", "pca_mlr")


# Plotting
# boxplot comparison and statistical testing
data_box <- data %>%
  group_by(classifier) %>%
  filter(tcr != 1) %>%
  filter(sample_size < 100)

filter <- data_box$clone_label == "persistent"
tp_box_accuracy <- ggplot(data_box[filter, ],
                          aes(x = classifier, y = accuracy,
                              fill = classifier)) +
  geom_hline(yintercept = 1 / 3, linetype = "dashed", color = "gray80") +
  geom_violin() +
  ggsignif::geom_signif(test = "wilcox.test",
                        comparisons = combinations,
                        margin_top = 0.05, tip_length = 0.01, vjust = 0.4,
                        map_signif_level = TRUE, textsize = 2.5,
                        step_increase = 0.05) +
  xlab("Classifier") +
  ylab("Accuracy") +
  labs(title = "Classifier performance \nsample size > 100 cells",
       fill = "Classifier") +
  scale_x_discrete(limits = cls_order) +
  theme_bw()
show(tp_box_accuracy)

tp_box_sensitivity <- ggplot(data_box,
                             aes(x = classifier, y = sensitivity,
                                 fill = classifier)) +
  geom_boxplot() +
  ggsignif::geom_signif(test = "wilcox.test",
                        comparisons = combinations,
                        margin_top = 0.05, tip_length = 0.01, vjust = 0.4,
                        map_signif_level = TRUE, textsize = 2.5,
                        step_increase = 0.05) +
  facet_wrap(~ clone_label, ncol = 3) +
  xlab("Classifier") +
  ylab("Sensitivity (True positive rate)") +
  labs(title = "Classifier performance \nsample size > 100 cells",
       fill = "Classifier") +
  scale_x_discrete(limits = cls_order) +
  theme_bw()
show(tp_box_sensitivity)

tp_box_specificity <- ggplot(data_box,
                             aes(x = classifier, y = specificity,
                                 fill = classifier)) +
  geom_boxplot() +
  ggsignif::geom_signif(test = "wilcox.test",
                        comparisons = combinations,
                        margin_top = 0.05, tip_length = 0.01, vjust = 0.4,
                        map_signif_level = TRUE, textsize = 2.5,
                        step_increase = 0.05) +
  facet_wrap(~ clone_label, ncol = 3) +
  xlab("Classifier") +
  ylab("Specificity (True negative rate)") +
  labs(title = "Classifier performance \nsample size > 100 cells",
       fill = "Classifier") +
  scale_x_discrete(limits = cls_order) +
  theme_bw()
show(tp_box_specificity)


# accuracy of classifier against sample size
robustness_data <- data %>%
  filter(sample_size > 100) %>%
  dplyr::group_by(sample_size, classifier) %>%
  dplyr::summarise(mean_ac = mean(accuracy),
                   sd_ac = sd(accuracy),
                   sem_ac = sd(accuracy) / sqrt(n()))

# plot robustnes data
rob_ac <- ggplot(data = robustness_data, aes(x = sample_size)) +
  geom_ribbon(aes(ymin = mean_ac - sd_ac, ymax = mean_ac + sd_ac,
                  color = classifier, fill = classifier), alpha = 0.2) +
  facet_wrap(~ classifier, ncol = 2) +
  geom_line(aes(y = mean_ac, color = classifier)) +
  xlab("Sample size") +
  ylab("Accuracy") +
  labs(title = "Classifier robustness",
       color = "Classifier", fill = "Classifier") +
  theme_bw()
show(rob_ac)


# robustness of classifier against sample size
robustness_data_by_label <- data %>%
  filter(sample_size > 100) %>%
  dplyr::group_by(sample_size, clone_label, classifier) %>%
  dplyr::summarise(mean_sens = mean(sensitivity),
                   sd_sens = sd(sensitivity),
                   mean_spec = mean(specificity),
                   sd_spec = sd(specificity))

# plot robustness data
rob_sens <- ggplot(data = robustness_data_by_label, aes(x = sample_size)) +
  geom_ribbon(aes(ymin = mean_sens - sd_sens, ymax = mean_sens + sd_sens,
                  color = clone_label, fill = clone_label), alpha = 0.5) +
  facet_wrap(~ classifier, ncol = 2) +
  geom_line(aes(y = mean_sens, color = clone_label)) +
  xlab("Sample size") +
  ylab("Sensitivity") +
  labs(title = "Classifier robustness",
       color = "clone label", fill = "clone label") +
  theme_bw()
show(rob_sens)

# plot robustness data
rob_spec <- ggplot(data = robustness_data_by_label, aes(x = sample_size)) +
  geom_ribbon(aes(ymin = mean_spec - sd_spec, ymax = mean_spec + sd_spec,
                  color = clone_label, fill = clone_label), alpha = 0.5) +
  facet_wrap(~ classifier, ncol = 2) +
  geom_line(aes(y = mean_spec, color = clone_label)) +
  xlab("Sample size") +
  ylab("Specificity") +
  labs(title = "Classifier robustness",
       color = "clone label", fill = "clone label") +
  theme_bw()
show(rob_spec)
