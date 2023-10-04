# TCR-Classifier
Lab rotation project to model the T-cell receptor (TCR) repertoire and classify clones based on their dyamic behaviour from sampled data
![Workflow diagram](/example/workflow.svg)

## Dependencies
The following R packages are utilisied in this project:
- adaptivetau (adaptivetau function for stochastic model computation)
- vegan (vegdist function for computation of dissimilarity indices)
- nnet (multinom function for multinomial logistic regression model)
- FNN (fast k-nearest neighbour statistic)
- caret (evaluation of classification models)
- utils (progress bar)
- tidyr (data manipulation)
- dplyr (data manipulation)
- ggplot2 (data visualisation)
- ggsignif (statistical significance annotation)

To install the required dependencies, run the following command in R:
```R
# install package
install.packages(c("adaptivetau", "vegan", "nnet", "FNN", "caret", "utils", "tidyr", "dplyr", "ggplot2", "ggsignif"))
```
## Usage
This collection of functions and scripts can be run from the command line without installation.

On a slurm cluster the scripts can be called using the `cluster/submit.sh` script: \
Make sure to edit the script according to your needs prior submition.
```bash
sbatch cluster/submit.sh <script_arguments>
```
For iterative job submition, the `cluster/iterations.job` script can be used: \
The script iterates through a loop and calls submit.sh passing to it the iteration index as argument.
```bash
./cluster/iterations.job
```

## Model simulation
The `tcr_model_simulate.R` script can be used to simulate the TCR repertoire model. \
The script simulates all TCR repertoire models and saves them as `.rds` files in the `data` folder. \
The script takes the following arguments:
- `num_generation`: number of identical simulations to run (for statistics)
- `num_clones`: number of clones (per label) in the repertoire
- `clone_size`: number of cells per clone
- `param_scale`: scaling factor for the model's time scale

The script outputs a list of TCR objects (saved as `.rds` files) containing the simulated data for each TCR simulation run.
Each TCR object contains the following data:
- `clonotype`: list of `clone` objects of the TCR repertoire
- `sim_times`: named vector of time points at which the TCR repertoire must be calculated (sampling times)
- `carry_cap`: carrying capacity of the TCR repertoire
- `clone_labels` : character vector of unique clone labels contained in the TCR repertoire
- `data`: matrix containing the simulated number of cells within each clone (columns) at each time point (rows). Also the `total` number of cells is calculated (last column).

Adittionally, each `clone` object has the following data:
- `clone_id`: unique identifier of the clone (index starting at 1)
- `label`: factor indicating the clone's label (e.g. `persisting`, `contracting`, `late_emerging`)
- `init_size`: initial size of the clone
- `params`: list of birth death rates for each time partition $\left[ t_{n-1}, t_n \right]$

## Model visualisation
The temporal composition of the TCR repertoire can be visualised using the function `tcr_temporal_comp`. \
![TCR temporal composition](/example/CloneEvo.png)

## Reference measurement calculation
Refence measurements for model training can be calculated using the `ref_measure_calc.R` script. \
This script samples cells from the simulated TCR repertoire, calculates the reference measurements and saves them as `.rds` files in the `data` folder. \
The script takes the following arguments:
- `sample_size`: number of cells to sample from the TCR repertoire at each time point
- `iteration`: number (index) of identical simulations to run (for statistics)

## Classification scripts
The scripts starting with `classify_` can be used to classify sampled clones. \
These scripts load a previously simulated TCR repertoire, the corresponding reference measurements, classify the clones using the method specified in the script's name and save the classification performance based on the confusion matrix. \
Currently there are four classification methods implemented:
- `naive`: naive classification based on the clone's incidence (last two time points)
- `mlr`: multinomial logistic regression model based on the entire measurement calculated from the sampled cells
- `pca_mlr`: multinomial logistic regression model based on a subset of the calculated measurements (pca)
- `knn`: k-nearest neighbour classification based on the entire measurement calculated from the sampled cells
- `pca_knn`: k-nearest neighbour classification based on a subset of the calculated measurements (pca)or 

By setting the `jackknife_only` parameter to `TRUE`, the measurements utilsed for classification can be reduced to the jackknife measurements only.

## Dimensionality reduction and feature extraction (pca)
Scripts using the `pca` method for dimensionality reduction and feature extraction are prefixed with `pca_`. \
These scripts use the `create_pca_mapper` closure to apply pca to the reference measurements and generate an object that can be used to map data based on a user defined linear combination of the measured variables. \
Hereby, it is up to the user to select two principle components and a number of features per principle component that yields the desired clustering. \
Subsequently, the function generated can be used to transform data onto the same linear combination calling the `$transform` method of the created mapper.

## Variable loadings and Reference measurement visualisation (after pca)
To help feature selection, the matrix of variable loadings can be visualised using `plt_pca_loadings` function. Further, the reference measurements (after the dimension reduction step) can also be visualised using the `plt_pca_features` function. \
![PCA loadings](/example/loadings.svg)
![PCA features](/example/CloneClustering2.svg)

Axis labels are stored in the mapper object generated by the `create_pca_mapper` closure.

## Classification performance visualisation
The classification performance can be visualised using the `plt_classifier.R` script. \
Currently, the following plots are implemented:
- classification accuracy, specificity and sensitivity by classifcation method
- classification accuracy, specificity and sensitivity by sample size and classification method

![Classifier performance](/example/cls_perf.svg)
![Classifier robustness](/example/cls_robust.svg)
