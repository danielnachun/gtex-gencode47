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

# General get_files_to_process function based on file_type
get_files_to_process() {
    case "${file_type:-bam_files}" in
        "bam_files")
            # Get all bam files from input directory
            input_path="${input_dir:-${bam_dir:-${realign_bam_dir}}}"
            file_pattern="${file_pattern:-Aligned.sortedByCoord.out.patched.v11md.bam$}"
            
            # Get all matching files using ls + grep
            bam_dir_bam_list=$(ls "$input_path" | grep "${file_pattern}")
            
            # Filter to only those in gtex_ids if specified
            if [ -n "${gtex_ids:-}" ]; then
                full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
            else
                full_bam_list="$bam_dir_bam_list"
            fi
            original_count=$(echo "$full_bam_list" | wc -l)
            echo "Original sample count: ${original_count}" >&2
            
            # Determine which files to process based on regenerate_all setting
            if [ "${regenerate_all:-false}" = true ]; then
                # Process all files
                echo "$full_bam_list" | sed "s|^|${input_path}/|"
            else
                # Only process files that don't have completion markers
                completion_dir="${output_dir}/completed/${completion_subdir:-${JOB_NAME}}"
                # Use find + grep -v for bulk file checking
                echo "$full_bam_list" | grep -v -F -f <(find "$completion_dir" -name "*.completed" -exec basename {} \; | sed 's|\.completed$||') | sed "s|^|${input_path}/|"
            fi
            ;;
        "ld_regions")
            # Process LD regions from BED file
            if [ -n "${ld_blocks_bed:-}" ]; then
                # Write a padded TSV (1-based inclusive) from the BED
                padded_ld_blocks_tsv="${output_dir:-${out_dir}}/padded_ld_blocks.tsv"
                awk -v PAD="${padding:-0}" -F '\t' '
                  NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
                    chrom=$1; sub(/^chr/,"",chrom);
                    from1=$2+1; to1=$3;
                    pf=(from1>PAD? from1-PAD:1);
                    pt=to1+PAD;
                    printf "%s\t%d\t%d\n", chrom, pf, pt;
                  }
                ' "$ld_blocks_bed" > "$padded_ld_blocks_tsv"
                
                original_count=$(wc -l < "${padded_ld_blocks_tsv}")
                echo "Original region count: ${original_count}" >&2
                
                if [ "${regenerate_all:-false}" = true ]; then
                    cat "$padded_ld_blocks_tsv"
                else
                    # Find existing files and use grep -v for bulk filtering
                    regions_to_process_tsv="${output_dir:-${out_dir}}/regions_to_process.tsv"
                    rm -f "$regions_to_process_tsv"
                    
                    # Use find + grep -v for bulk file checking
                    cat "$padded_ld_blocks_tsv" | grep -v -F -f <(find "${output_dir:-${out_dir}}" -name "LD_chr*.ld.gz" -exec basename {} \; | sed 's|LD_chr\([0-9]*\)_\([0-9]*\)_\([0-9]*\)\.ld\.gz|\1\t\2\t\3|') > "$regions_to_process_tsv"
                    cat "$regions_to_process_tsv" 2>/dev/null || true
                fi
            else
                echo "Error: ld_blocks_bed not specified for ld_regions file_type"
                exit 1
            fi
            ;;
        "participants")
            # Process participants from participant list
            if [ -n "${participant_id_list:-}" ]; then
                # Check for duplicates
                if [ $(sort "${participant_id_list}" | uniq -d | wc -l) -gt 0 ]; then
                    echo "Error: Duplicate entries found in ${participant_id_list}:"
                    sort "${participant_id_list}" | uniq -d
                    exit 1
                fi
                
                if [ "${regenerate_all:-false}" = true ]; then
                    cat "${participant_id_list}"
                else
                    # Find existing files and use grep -v for bulk filtering
                    new_participant_list="${output_dir}/participants_to_process.txt"
                    > "${new_participant_list}"
                    
                    # Find all existing participant files and extract participant IDs
                    existing_participants=$(find "$output_dir" -name "*.snps.vcf.gz" -exec basename {} \; | sed 's|\.snps\.vcf\.gz$||' | sort -u)
                    
                    # Use grep -v for bulk filtering
                    cat "${participant_id_list}" | grep -v -F -f <(echo "$existing_participants") > "${new_participant_list}"
                    cat "${new_participant_list}"
                fi
            else
                echo "Error: participant_id_list not specified for participants file_type"
                exit 1
            fi
            ;;
        "tissues")
            # Process tissues from tissue list
            if [ -n "${tissue_id_list:-}" ]; then
                cat "${tissue_id_list}"
            else
                echo "Error: tissue_id_list not specified for tissues file_type"
                exit 1
            fi
            ;;
        "ld_regions_string")
            # Process LD regions as region strings (for coloc)
            if [ -n "${ld_region_list:-}" ]; then
                ld_regions=$(awk 'NR>1 {printf "%s:%d-%d\n", $1, $2+1, $3}' "${ld_region_list}")
                total_count=$(wc -l < "${ld_region_list}")
                
                if [ "${regenerate_all:-FALSE}" = "TRUE" ] || [ "${regenerate_all:-FALSE}" = "true" ]; then
                    echo "$ld_regions"
                else
                    # Find existing completion files and use grep -v for bulk filtering
                    completion_dir="${output_dir}/completed"
                    # Use find + grep -v for bulk file checking
                    echo "$ld_regions" | grep -v -F -f <(find "$completion_dir" -name "*.completed" -exec basename {} \; | sed 's|\.completed$||')
                fi
            else
                echo "Error: ld_region_list not specified for ld_regions_string file_type"
                exit 1
            fi
            ;;
        *)
            echo "Error: Unknown file_type '${file_type}'. Supported types: bam_files, ld_regions, participants, tissues, ld_regions_string"
            exit 1
            ;;
    esac
}

# Check if files_to_process is empty
files_to_process=$(get_files_to_process)
if [ -z "$files_to_process" ]; then
    echo "To be processed: 0"
    echo "All files processed"
    exit 0
else
    # Count files to process
    to_process_count=$(echo "$files_to_process" | wc -l)
    
    # Calculate original count for display (this is a bit hacky but works)
    case "${file_type:-bam_files}" in
        "bam_files")
            # Reuse the same logic as in get_files_to_process for consistency
            input_path="${input_dir:-${bam_dir:-${realign_bam_dir}}}"
            file_pattern="${file_pattern:-Aligned.sortedByCoord.out.patched.v11md.bam$}"
            
            bam_dir_bam_list=$(ls "$input_path" | grep "${file_pattern}")
            if [ -n "${gtex_ids:-}" ]; then
                full_bam_list=$(grep -F -f "$gtex_ids" <<< "$bam_dir_bam_list")
            else
                full_bam_list="$bam_dir_bam_list"
            fi
            original_count=$(echo "$full_bam_list" | wc -l)
            ;;
        *)
            # For other types, we can't easily calculate original count
            original_count=$to_process_count
            ;;
    esac
    
    completed_count=$((original_count - to_process_count))
    echo "Already completed: ${completed_count}"
    echo "To be processed: ${to_process_count}"
    
    # Calculate batch configuration to fit within max_array_size constraint
    # Each array job processes: num_parallel files per step
    
    # Calculate minimum files per job to stay within max_array_size limit
    num_series=$(( (to_process_count + max_array_size - 1) / max_array_size ))
    num_series=$(( (num_series + num_parallel - 1) / num_parallel ))
    echo "Configured series steps per job: ${num_series}"
    
    # Calculate number of sequential steps per job
    num_batches=$(( (to_process_count + num_parallel * num_series - 1) / (num_parallel * num_series) ))
    echo "Configured array jobs: ${num_batches}"
    echo "Configured parallel files per step: ${num_parallel}"
    echo "Total files per array job: $((num_parallel * num_series))"
fi

# Create batch files for SLURM array processing
file_list_folder="$output_dir/file_lists/file_lists_${JOB_NAME}"
rm -rf "${file_list_folder}"
mkdir -p "${file_list_folder}"

# Calculate files per batch file to distribute all files evenly
final_files_per_batch=$(( (to_process_count + num_batches - 1) / num_batches ))
echo "Final files per batch file: ${final_files_per_batch}"

# Split files into batch files
echo "$files_to_process" | split -l "${final_files_per_batch}" --additional-suffix=".txt" - "${file_list_folder}/file_list_"

# Create list of batch file paths for SLURM array
file_list_paths="${output_dir}/file_lists/file_list_paths_${JOB_NAME}.txt"
rm -rf "${file_list_paths}"
printf "%s\n" "${file_list_folder}"/* > "${file_list_paths}"

echo "Batches created: ${num_batches}"
ideal_array_jobs=$(( (to_process_count + num_parallel * num_series - 1) / (num_parallel * num_series) ))
echo "Ideal array jobs needed: ${ideal_array_jobs}"
echo "Actual array jobs created: ${num_batches}"

# Configure SLURM job parameters
sbatch_params=(
    --output "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --error "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --array "1-${num_batches}%250"
    --cpus-per-task "${num_parallel}"
    --time "${job_time:-12:00:00}"
    --mem "${job_mem:-256G}"
    --job-name "${JOB_NAME}_batch"
    "${BATCH_SCRIPT}"
    "${file_list_paths}"
    "${num_parallel}"
    "${code_dir}/02.A_realign_bam_sj.sh"
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
    # Additional parameters for sherlock
    sbatch \
        --partition normal,owners \
        --tmp 200G \
        "${sbatch_params[@]}"
elif [ "${submit_on}" = 'scg' ]; then
    # Additional parameters for scg
    sbatch \
        --account smontgom \
        --partition batch \
        --constraint "nvme" \
        "${sbatch_params[@]}"
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
