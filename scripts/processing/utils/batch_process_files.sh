#!/usr/bin/env bash

# Shared batch processing script for parallel/series processing
# Usage: batch_process_files.sh <file_list_paths> <num_parallel> <processing_script> [additional_params...]

set -o xtrace -o nounset -o pipefail -o errexit

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

# Define all possible parameters
options_array=(
    file_list_paths
    num_parallel
    processing_script
    file_param
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Save original arguments to preserve additional parameters with their option names
ORIGINAL_ARGS=("$@")

# Parse command line arguments with getopt
# Suppress stderr warnings about unrecognized options (they're passed through after --)
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'batch_process_files' -- "$@" 2>/dev/null)
eval set -- "${arguments}"

# Initialize optional variables
file_param=${file_param-}

while true; do
    case "${1}" in
        --file_list_paths )
            file_list_paths="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --num_parallel )
            num_parallel="${2}"; shift 2 ;;
        --processing_script )
            processing_script="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --file_param )
            file_param="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

# Extract additional parameters from original arguments, preserving option names
# Find the -- separator and take everything after it
ADDITIONAL_PARAMS=()
found_sep=false
for arg in "${ORIGINAL_ARGS[@]}"; do
    if [ "$arg" = "--" ]; then
        found_sep=true
        continue
    fi
    if [ "$found_sep" = true ]; then
        ADDITIONAL_PARAMS+=("$arg")
    fi
done

# Verify that we found the separator and have additional params
if [ "$found_sep" = false ]; then
    echo "Error: Expected -- separator before additional parameters" >&2
    exit 1
fi

# Get the paths to the files for this batch
line_number="${SLURM_ARRAY_TASK_ID}"
file_list="$(sed "${line_number}q; d" "${file_list_paths}")"

# Extract and copy reference directory to local temp if provided
# Convert --reference_dir (global) to --local_reference_dir (local copy)
for i in "${!ADDITIONAL_PARAMS[@]}"; do
    if [ "${ADDITIONAL_PARAMS[$i]}" = "--reference_dir" ] && [ -n "${ADDITIONAL_PARAMS[$((i+1))]:-}" ]; then
        original_ref_dir="${ADDITIONAL_PARAMS[$((i+1))]}"
        local_ref_dir="${TMPDIR}/references_${line_number}"
        mkdir -p "${local_ref_dir}"
        rsync -PrhLtv "${original_ref_dir}"/* "${local_ref_dir}/"
        ADDITIONAL_PARAMS[$i]="--local_reference_dir"
        ADDITIONAL_PARAMS[$((i+1))]="${local_ref_dir}"
        break
    fi
done

echo "Parallel jobs: ${num_parallel}"
echo "Requested CPUs: ${SLURM_CPUS_PER_TASK}"

# Process all files in batch file continuously with num_parallel jobs running
echo "Processing all files with ${num_parallel} parallel jobs"

# Use file_param if provided, otherwise default to --bam_file
file_param="${file_param:---bam_file}"
cat "${file_list}" | parallel -j"${num_parallel}" --ungroup --verbose \
    "${processing_script}" \
    "${ADDITIONAL_PARAMS[@]}" \
    "${file_param}" {}

echo "Batch finished"

