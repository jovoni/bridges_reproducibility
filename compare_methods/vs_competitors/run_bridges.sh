#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --array=1-60%10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/bridges_out_%A_%a.log
#SBATCH --error=logs/bridges_err_%A_%a.log
#SBATCH --job-name=bridges_fit

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules (if needed)
#source ~/.bashrc
#conda activate R_env

source ~/.bashrc
conda activate process

Rscript run_bridges.R