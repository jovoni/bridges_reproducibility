#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=50GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -o logs/bridges_fit_%A_%a.out
#SBATCH -e logs/bridges_fit_%A_%a.err
#SBATCH --job-name=bridges_fit
#SBATCH --array=1-35

mkdir -p logs results/bridges_trees

# activate your env
source ~/.bashrc
conda activate process

# run the indexed fit
Rscript run_fit.R "${SLURM_ARRAY_TASK_ID}"
