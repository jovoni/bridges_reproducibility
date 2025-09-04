#!/bin/bash
#SBATCH --job-name=bfbsim_array
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/medicc_%A_%a.log
#SBATCH --error=logs/medicc_%A_%a.err
#SBATCH --array=1-60%24

# Load required modules
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate medicc_env

# Create logs directory if it doesn't exist
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
input_file="${dir_name}/medicc_input.csv"
output_dir="results/${sim_id}/medicc2"
mkdir -p "$output_dir"

echo "Running medicc2 on $input_file -> $output_dir"

# Run medicc2
start_time=$(date +%s.%N)
medicc2 "$input_file" "$output_dir"
end_time=$(date +%s.%N)

execution_time=$(echo "$end_time - $start_time" | bc)
echo "$execution_time" > "$output_dir/medicc_time.txt"
echo "Execution time: ${execution_time} seconds (saved to $output_dir/medicc_time.txt)"
