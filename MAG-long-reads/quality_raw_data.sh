#!/bin/bash

#SBATCH --job-name=qc_raw_data_bats
#SBATCH --output=%x.job%j.out
#SBATCH --error=%x.job%j.err
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=40G
#SBATCH --time=04:00:00

# --- LOADING MODULES ---

module purge
module load singularity
module load fastqc/0.12.1


SINGULARITY_CONTAINER_HOME=$HOME/singularity_containers


INPUT_FILE_FULL_PATH=$(realpath "$1")
INPUT_FILENAME=$(basename "$INPUT_FILE_FULL_PATH")
INPUT_DIR=$(dirname "$INPUT_FILE_FULL_PATH")
LOCAL_BASE_DIR=$(dirname "$INPUT_DIR")

FASTQC_OUTPUT_DIR="$LOCAL_BASE_DIR/fastqc_results"
FILTERED_FILENAME="filtred_$INPUT_FILENAME"
LOCAL_OUTPUT_PATH="$INPUT_DIR/$FILTERED_FILENAME"

mkdir -p "$FASTQC_OUTPUT_DIR"

echo "================================================="
echo "Starting workflow: Quality control of raw data"
echo "Sample: $INPUT_FILENAME"
echo "================================================="

# --- EXECUTION ---

# [1/3] FastQC 
echo "[1/3] Running FastQC (Pre-filter)..."
fastqc --nano --memory 1000 -t 15 -o "$FASTQC_OUTPUT_DIR" "$INPUT_FILE_FULL_PATH"

# [2/3] NanoFilt and Gzip (Singularity)
echo "[2/3] Running NanoFilt..."
gunzip -c "$INPUT_FILE_FULL_PATH" | \
    singularity run $SINGULARITY_CONTAINER_HOME/nanofilt.sif NanoFilt -q 7 | \
    gzip > "$LOCAL_OUTPUT_PATH"

# [3/3] FastQC (Post-filter)
echo "[3/3] Running FastQC (Post-filter)..."
fastqc --nano --memory 1000 -t 15 -o "$FASTQC_OUTPUT_DIR" "$LOCAL_OUTPUT_PATH"

# --- REPORTING RESOURCES ---
echo -e "\n================================================="
echo "End Time: $(date)"
sh -c "seff $SLURM_JOBID"
echo "================================================="
