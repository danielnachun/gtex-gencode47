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
mkdir -p "${output_dir}/completed/quantify_bam"

# Submit batch jobs: helper will filter completed items and divide into batches
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    --config_file "${CONFIG_FILE}" \
    --batch_script "${code_dir}/utils/batch_process_files.sh" \
    --job_name "quantify_bam" \
    --items_list_file "${sample_list_file}" \
    --processing_script "${code_dir}/03_quantify_bam.sh" \
    --file_param "--sample_id" \
    --completion_dir "${output_dir}/completed/quantify_bam" \
    --regenerate_all "${regenerate_all:-false}" \
    -- \
    --reference_dir "${reference_dir}" \
    --vcf_dir "${vcf_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --bam_dir "${bam_dir}" \
    --completion_dir "${output_dir}/completed/quantify_bam" \
    --reference_fasta "${reference_fasta}" \
    --genes_gtf "${genes_gtf}" \
    --chr_sizes "${chr_sizes}" \
    --intervals_bed "${intervals_bed}" \
    --ipa_annotation "${ipa_annotation}" \
    --editing_bed "${editing_bed}" \
    --bam_file_end "${bam_file_end}"
