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

# Parse command line arguments
file_list_paths="${1}"
num_parallel="${2}"
processing_script="${3}"
shift 3  # Remove first 3 arguments, rest are additional parameters

check_for_file "file_list_paths" "${file_list_paths}"
check_for_file "processing_script" "${processing_script}"

# Get the paths to the files for this batch
line_number="${SLURM_ARRAY_TASK_ID}"
file_list="$(sed "${line_number}q; d" "${file_list_paths}")"

# Copy over references to this node (if reference_dir is provided)
if [ -n "${reference_dir:-}" ]; then
    reference_dir_prefix="${TMPDIR}/references_${line_number}"
    mkdir -p "${reference_dir_prefix}"
    rsync -PrhLtv "${reference_dir}"/* "${reference_dir_prefix}/"
fi

echo "Parallel jobs: ${num_parallel}"
echo "Requested CPUs: ${SLURM_CPUS_PER_TASK}"

# Process all files in batch file continuously with num_parallel jobs running
echo "Processing all files with ${num_parallel} parallel jobs"

# Determine parameter name based on processing script or file_type
# For coloc script (05_colocalize_regions.sh), use --ld_region parameter
# For other scripts, use --bam_file (default) or FILE_PARAM if set
if [[ "${processing_script}" == *"05_colocalize_regions.sh" ]] || [ "${file_type:-}" = "ld_regions_string" ]; then
    # For coloc, the file list contains region strings (e.g., "chr1:1000-2000")
    # Pass each region string with --ld_region parameter
    cat "${file_list}" | parallel -j"${num_parallel}" --ungroup --verbose \
        "${processing_script}" \
        "$@" \
        --ld_region {}
else
    # Default behavior: pass files with --bam_file parameter
    file_param="${FILE_PARAM:---bam_file}"
    cat "${file_list}" | parallel -j"${num_parallel}" --ungroup --verbose \
        "${processing_script}" \
        "$@" \
        "${file_param}" {}
fi

echo "Batch finished"

