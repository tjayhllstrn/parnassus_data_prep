#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --array=0-6
#SBATCH --mem-per-cpu=5000
#SBATCH --nodes=2
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=processing_parnassus_v0
#SBATCH --time=4:00:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/parnassus_data_prep

INPUT_DIR=$1
DATA_VERSION=$2
OUTPUT_DIR="/volatile/clas12/users/tjhellst/parnassus_data/$2/"
INPUT_FILE_BASE="*-030*.hipo"

# Enable nullglob to get empty array if no matches (instead of literal pattern)
shopt -s nullglob
FILES=($INPUT_DIR/$INPUT_FILE_BASE) 
shopt -u nullglob

# Debug output
echo "INPUT_DIR: $INPUT_DIR"
echo "INPUT_FILE_BASE: $INPUT_FILE_BASE"
echo "Pattern: $INPUT_DIR/$INPUT_FILE_BASE"
echo "Number of files found: ${#FILES[@]}"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Check if files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No files matched pattern $INPUT_DIR/$INPUT_FILE_BASE"
    exit 1
fi

echo "Running job ${SLURM_ARRAY_TASK_ID[@]} on file: ${FILES[$SLURM_ARRAY_TASK_ID]}"

root -l -b -q "macros/hipo2tree_parnassus.C(\"${FILES[$SLURM_ARRAY_TASK_ID]}\",\"$OUTPUT_DIR\",-1)"
