#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Generate list of all sample IDs to process
sample_list_file="${output_dir}/file_lists/sample_list.txt"
mkdir -p "$(dirname "${sample_list_file}")"
bam_file_end="${bam_file_end:-Aligned.sortedByCoord.out.patched.v11md.bam}"

# Get list of all BAM filenames, filtering by gtex_ids if specified, and extract sample IDs
if [ -n "${gtex_ids:-}" ] && [ -f "${gtex_ids}" ]; then
    ls "${bam_dir}" | grep "${bam_file_end}$" | grep -F -f "${gtex_ids}" | sed "s/.${bam_file_end}$//" > "${sample_list_file}"
else
    ls "${bam_dir}" | grep "${bam_file_end}$" | sed "s/.${bam_file_end}$//" > "${sample_list_file}"
fi
echo "Created sample ID list: ${sample_list_file}"

# Create completion directory
mkdir -p "${output_dir}/completed/call_mutect"

# Submit batch jobs: helper will filter completed items and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "call_mutect" \
    --items_list_file "${sample_list_file}" \
    --processing_script "${code_dir}/03.B_call_mutect.sh" \
    --file_param "--sample_id" \
    --completion_dir "${output_dir}/completed/call_mutect" \
    --regenerate_all "${regenerate_all:-false}" \
    -- \
    --reference_dir "${reference_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --bam_dir "${bam_dir}" \
    --completion_dir "${output_dir}/completed/call_mutect" \
    --vcf_dir "${vcf_dir}" \
    --reference_fasta "${reference_fasta}" \
    --dbsnp "${dbsnp}" \
    --indels_mills "${indels_mills}" \
    --indels_decoy "${indels_decoy}" \
    --gene_intervals_bed "${gene_intervals_bed}" \
    --full_vcf_file "${full_vcf_file}" \
    --exac_reference "${exac_reference}" \
    --bam_file_end "${bam_file_end}"
