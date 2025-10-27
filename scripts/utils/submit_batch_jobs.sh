#!/usr/bin/env bash

# Shared utility for submitting batch jobs with parallel/series processing
# Usage: submit_batch_jobs.sh <config_file> <batch_script> <job_name> [additional_params...]

set -o xtrace -o nounset -o errexit

# Parse arguments
CONFIG_FILE="${1}"
BATCH_SCRIPT="${2}"
JOB_NAME="${3}"
shift 3  # Remove first 3 arguments, rest are additional parameters

# Source the config file
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Get all files to process (this will be overridden by each script)
# Default implementation - scripts should override this
get_files_to_process() {
    echo "Error: get_files_to_process function not implemented in config"
    exit 1
}

# Check if files_to_process is empty
if [ -z "$(get_files_to_process)" ]; then
    echo "To be processed: 0"
    echo "All files processed"
    exit 0
else
    # Count files to process
    to_process_count=$(get_files_to_process | wc -l)
    echo "To be processed: $to_process_count"
    
    # Calculate batch configuration to fit within max_array_size constraint
    # Each array job processes: num_parallel files per step
    
    # Calculate minimum files per job to stay within max_array_size limit
    num_series=$(( (to_process_count + max_array_size - 1) / max_array_size + num_parallel - 1 ) / num_parallel ))
    
    # Calculate number of sequential steps per job
    num_batches=$(( (to_process_count + num_parallel * num_series - 1) / (num_parallel * num_series) ))
    
    echo "Configured array jobs: ${num_batches}"
    echo "Configured parallel files per step: ${num_parallel}"
    echo "Configured series steps per job: ${num_series}"
    echo "Total files per array job: $((num_parallel * num_series))"
fi

# Create batch files for SLURM array processing
file_list_folder="$output_dir/file_lists/file_lists_${JOB_NAME}"
rm -rf "${file_list_folder}"
mkdir -p "${file_list_folder}"

# Calculate files per batch file to distribute all files evenly
final_files_per_batch=$(( (to_process_count + num_batches - 1) / num_batches ))

# Split files into batch files
get_files_to_process | split -l "${final_files_per_batch}" --additional-suffix=".txt" - "${file_list_folder}/file_list_"

# Create list of batch file paths for SLURM array
file_list_paths="${output_dir}/file_lists/file_list_paths_${JOB_NAME}.txt"
rm -rf "${file_list_paths}"
printf "%s\n" "${file_list_folder}"/* > "${file_list_paths}"

echo "Batches created: ${num_batches}"
echo "Already completed: $((original_count - to_process_count))"
echo "To be processed: ${to_process_count}"
echo "Ideal array jobs needed: $(( (to_process_count + num_parallel * num_series - 1) / (num_parallel * num_series) ))"
echo "Actual array jobs created: ${num_batches}"
echo "Final files per batch file: ${final_files_per_batch}"

# Configure SLURM job parameters
sbatch_params=(
    --output "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --error "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --array "1-${num_batches}%250"
    --cpus-per-task "${num_parallel}"
    --time "${job_time:-40:00:00}"
    --mem "${job_mem:-256G}"
    --job-name "${JOB_NAME}_batch"
    "${BATCH_SCRIPT}"
    --file_list_paths "${file_list_paths}"
    --num_parallel "${num_parallel}"
    "$@"  # Pass through additional parameters
)

# Add any additional SLURM parameters from config
if [ -n "${job_partition:-}" ]; then
    sbatch_params+=(--partition "${job_partition}")
fi
if [ -n "${job_account:-}" ]; then
    sbatch_params+=(--account "${job_account}")
fi
if [ -n "${job_constraint:-}" ]; then
    sbatch_params+=(--constraint "${job_constraint}")
fi
if [ -n "${job_tmp:-}" ]; then
    sbatch_params+=(--tmp "${job_tmp}")
fi

# Submit job array to cluster
if [ "${submit_on}" = 'sherlock' ]; then
    # Additional parameters for sherlock (if not specified in config)
    if [ -z "${job_partition:-}" ]; then
        sbatch_params+=(--partition normal,owners)
    fi
    if [ -z "${job_tmp:-}" ]; then
        sbatch_params+=(--tmp 200G)
    fi
    sbatch "${sbatch_params[@]}"
elif [ "${submit_on}" = 'scg' ]; then
    # Additional parameters for scg (if not specified in config)
    if [ -z "${job_account:-}" ]; then
        sbatch_params+=(--account smontgom)
    fi
    if [ -z "${job_partition:-}" ]; then
        sbatch_params+=(--partition batch)
    fi
    if [ -z "${job_constraint:-}" ]; then
        sbatch_params+=(--constraint "nvme")
    fi
    sbatch "${sbatch_params[@]}"
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
