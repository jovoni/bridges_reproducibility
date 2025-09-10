#!/bin/bash
#SBATCH --job-name=bfbsim_array
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=logs/sitka_%A_%a.log
#SBATCH --error=logs/sitka_%A_%a.err
#SBATCH --array=1-60%24

source ~/.bashrc
conda activate process
module load java/8-8u402b06

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
input_file_sitka="${dir_name}/sitka_input.csv"
output_dir="results/${sim_id}/sitka"
mkdir -p "$output_dir"

Rscript sitka_prep.R -i "$input_file_medicc" -o "$input_file_sitka"

cd "$output_dir"
input_file_sikta="../../../${dir_name}/sitka_input.csv"

start_time=$(date +%s.%N)
/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/sitkatree/sitka/build/install/nowellpack/bin/corrupt-straighten --input "$input_file_sikta" --neighborhoodSize 2
cp results/latest/output.csv ./

/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/sitkatree/sitka/build/install/nowellpack/bin/corrupt-filter --input output.csv --lowerFraction 0.05
cp results/latest/filtered.csv ./

/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/sitkatree/sitka/build/install/nowellpack/bin/corrupt-infer-with-noisy-params \
    --model.globalParameterization true \
    --model.binaryMatrix filtered.csv \
    --model.fprBound 0.1 \
    --model.fnrBound 0.5 \
    --engine PT \
    --engine.initialization FORWARD \
    --engine.nScans 1000 \
    --engine.nPassesPerScan 1 \
    --engine.nChains 1
cp results/latest/samples/phylo.csv ./

/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/sitkatree/sitka/build/install/nowellpack/bin/corrupt-average --csvFile phylo.csv --logisticTransform false
cp results/latest/average.csv ./

/orfeo/cephfs/scratch/cdslab/gsantacatterina/bridges_reproducibility/compare_methods/vs_competitors/sitkatree/sitka/build/install/nowellpack/bin/corrupt-greedy --tipInclusionProbabilities ReadOnlyCLMatrix average.csv
cp results/latest/tree.newick ./
end_time=$(date +%s.%N)

execution_time=$(echo "$end_time - $start_time" | bc)
echo "$execution_time" > "time.txt"
echo "Execution time: ${execution_time} seconds"
