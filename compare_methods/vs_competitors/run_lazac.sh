#!/bin/bash
#SBATCH --job-name=bfbsim_array
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/lazac_%A_%a.log
#SBATCH --error=logs/lazac_%A_%a.err
#SBATCH --array=1-60%24

source ~/.bashrc
conda activate process

mkdir -p logs

# Get the current job's parameters from the CSV file
# SLURM_ARRAY_TASK_ID corresponds to the line number (1-indexed, excluding header)
line_num=$((SLURM_ARRAY_TASK_ID + 1))  # +1 to skip header
params=$(sed -n "${line_num}p" ../data/param_grid.csv)

# Parse parameters
IFS=',' read -r sim_id seed bfb_rate n_events max_cells <<< "$params"

echo "Processing job $SLURM_ARRAY_TASK_ID: sim_id=$sim_id"

# Set up directories and files
dir_name="../data/${sim_id}"
input_file_medicc="${dir_name}/medicc_input.csv"
input_file_lazac="${dir_name}/lazac_input.csv"
output_dir="results/${sim_id}/lazac/"
mkdir -p "$output_dir"

Rscript lazac_prep.R -i "$input_file_medicc" -o "$input_file_lazac"

start_time=$(date +%s.%N)
/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/lazac-copy-number/build/src/lazac nni "$input_file_lazac" -a 2 -o "$output_dir"
end_time=$(date +%s.%N)

execution_time=$(echo "$end_time - $start_time" | bc)
echo "$execution_time" > "${output_dir}/time.txt"
echo "Execution time: ${execution_time} seconds"
