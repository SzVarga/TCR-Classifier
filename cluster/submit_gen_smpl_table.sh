#!/bin/bash

# Check if exactly two arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <path_to_tcr_collection> <sample_size>"
  exit 1
fi

# Path to the Project on host
TCR_PATH="/net/home.isilon/ag-pahle/Szilard/TCR"

# Path to the Project on container
CTR_PATH="/home/bq_szvarga"

# Assign arguments
TCR_DATA=$1
SMPL_SIZE=$2

# Path to prpoject related folders
DATA_PATH="${TCR_PATH}/data"
R_PATH="${TCR_PATH}/R"
SCRIPT_PATH="${TCR_PATH}/script"

# Path to the Singularity container
SIF_PATH="${TCR_PATH}/cluster/container_TCR.sif"

# Slurm job settings (modify according to your cluster)
SBATCH_JOB_NAME="gen_sample_table"   # Name for the Slurm job
SBATCH_PARTITION="single"            # Slurm partition to submit to
SBATCH_TIME="14-00:00:00"            # Job run time limit (adjust as necessary)
SBATCH_CPUS=1                        # Number of CPUs to allocate
SBATCH_MEM="4G"                      # Memory allocation

# Submit NUM_JOBS parallel jobs using a job array in Slurm
sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SBATCH_JOB_NAME}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --time=${SBATCH_TIME}
#SBATCH --cpus-per-task=${SBATCH_CPUS}
#SBATCH --mem=${SBATCH_MEM}

# Load any necessary modules (if required)
module load system/singularity/3.10.4

# Run the Singularity container with folder bindings, script argument, and index
singularity exec \
    --bind "${DATA_PATH}:${CTR_PATH}/data" \
    --bind "${R_PATH}:${CTR_PATH}/R" \
    --bind "${SCRIPT_PATH}:${CTR_PATH}/script" \
    "${SIF_PATH}" \
    Rscript ${CTR_PATH}/script/sample_population.R "${TCR_DATA}" "${SMPL_SIZE}"
EOF
