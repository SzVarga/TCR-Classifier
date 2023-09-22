#!/bin/sh

#SBATCH --job-name=TCR
#SBATCH -t 3-00:00:00
#SBATCH --partition=single

# load R module
# module load math/R/4.0.0

# submit job with specific config file
# Rscript ../script/tcr_model_simulate.R $1 $2 $3 $4
# Rscript ../script/ref_measure_calc.R $1 $2
# Rscript ../script/classify_naive.R $1
# Rscript ../script/classify_pca_knn.R $1
# Rscript ../script/classify_pca_mlr.R $1
# Rscript ../script/classify_mlr.R $1 
