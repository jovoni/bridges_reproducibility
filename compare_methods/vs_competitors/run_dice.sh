#!/bin/bash
#SBATCH --job-name=bfbsim_array
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/dice.log

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate medicc_env

# 5. RUN MEDICC2 WITH PROPER ENVIRONMENT
tail -n +2 ../data/param_grid.csv | while IFS=',' read -r sim_id seed bfb_rate n_events max_cells; do
    dir_name="../data/${sim_id}"
    input_file="${dir_name}/dice_input.tsv"
    output_dir="results/${sim_id}/dice"

    mkdir -p "$output_dir"
    echo "Running DICE on $input_file -> $output_dir"

    # Run dice
    start_time=$(date +%s.%N)
    dice -i $input_file -o $output_dir -m balME -s
    end_time=$(date +%s.%N)

    execution_time=$(echo "$end_time - $start_time" | bc)
    echo "$execution_time" > "$output_dir/dice_time.txt"
    echo "Execution time: ${execution_time} seconds (saved to $output_dir)"
done
