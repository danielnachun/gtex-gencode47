#!/usr/bin/env bash

set -o nounset -o errexit

# Source the config file
CONFIG_FILE="${1:-}"
if [[ -z "${CONFIG_FILE}" ]]; then
    echo "Usage: $(basename "$0") <config_file>"
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

echo "Submitting job for realigning STAR outputs"
# Call the shared batch submission utility
exec "${code_dir}/utils/submit_batch_jobs.sh" \
    "${CONFIG_FILE}" \
    "${code_dir}/utils/batch_process_files.sh" \
    "realign_bam_sj" \
    "${code_dir}/02.A_realign_bam_sj.sh" \
    --local_reference_dir "${reference_dir}" \
    --vcf_dir "${vcf_dir}" \
    --output_dir "${output_dir}" \
    --code_dir "${code_dir}" \
    --reference_fasta "${reference_fasta}" \
    --star_index "${star_index}"

echo "Done"
