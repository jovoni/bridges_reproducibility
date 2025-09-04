
# Iterate through each row in the parameter grid (excluding header)
tail -n +2 ../data/param_grid.csv | while IFS=',' read -r sim_id seed bfb_rate n_events max_cells; do
    dir_name="../data/${sim_id}"
    input_file="${dir_name}/medalt_input.tsv"
    output_dir="results/${sim_id}"

    mkdir -p "$output_dir"
    echo "Running Hamming and Euclidean building on $input_file -> $output_dir"

    # Run the R script
    start_time=$(date +%s.%N)
    Rscript run_hamming_and_euclidean.R "$input_file" "$output_dir"
    end_time=$(date +%s.%N)

    # Record runtime
    execution_time=$(echo "$end_time - $start_time" | bc)
    echo "Execution time: ${execution_time} seconds"
done
