#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# source the config file
CONFIG_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/config/04_ld_gtex_eur.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

mkdir -p ${out_dir}/logs


# Prepare PLINK 1.9-compatible keep file (FID IID)
reformated_sample_ids="$out_dir/european_subject_ids.keep.txt"
awk 'NF{print $1"\t"$1}' "$sample_ids" > "$reformated_sample_ids"

# Write a padded TSV (1-based inclusive) from the BED without converting to region strings
# these are the regions to process
padded_ld_blocks_tsv="$out_dir/padded_ld_blocks.tsv"
awk -v PAD="$padding" -F '\t' '
  NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    chrom=$1; sub(/^chr/,"",chrom);
    from1=$2+1; to1=$3;
    pf=(from1>PAD? from1-PAD:1);
    pt=to1+PAD;
    printf "%s\t%d\t%d\n", chrom, pf, pt;
  }
' "$ld_blocks_bed" > "$padded_ld_blocks_tsv"

# Get total number of regions
original_count=$(wc -l < "${padded_ld_blocks_tsv}")
echo "Original region count: ${original_count}"

# if true do all regions
# if false, do all regions that do not already have a completed LD file
regenerate_all=${regenerate_all:-false}
if [ "${regenerate_all}" = true ]; then
    # run all regions
    regions_to_process_tsv="$padded_ld_blocks_tsv"
    regions_count=${original_count}
else
    # only process regions that don't already have completed LD files
    regions_to_process_tsv="$out_dir/regions_to_process.tsv"
    # Remove old regions_to_process file to ensure fresh filtering
    rm -f "$regions_to_process_tsv"
    # Filter out regions that already have completed LD files
    while IFS=$'\t' read -r chr_id from_bp to_bp; do
        out_file="${out_dir}/LD_chr${chr_id}_${from_bp}_${to_bp}.ld.gz"
        if [[ ! -f "$out_file" ]]; then
            echo -e "${chr_id}\t${from_bp}\t${to_bp}" >> "$regions_to_process_tsv"
        fi
    done < "$padded_ld_blocks_tsv"
    regions_count=$(wc -l < "$regions_to_process_tsv" 2>/dev/null || echo "0")
fi

# Check if regions_to_process is empty
if [ "${regions_count}" -eq 0 ]; then
    echo "To be processed: 0"
    echo "All LD regions completed"
    exit 0  # Exit the script with a success status
else
    echo "To be processed: ${regions_count}"
fi

completed_count=$((original_count - regions_count))

if [[ ${regions_count} -gt ${max_array_size} ]]; then
  regions_count=${max_array_size}
fi

echo "Already completed: ${completed_count}"
echo "To be processed: ${regions_count}"


sbatch_params=(
    --output "${out_dir}/logs/generate_ld/%A/%A_%a.log"
    --error "${out_dir}/logs/generate_ld/%A/%A_%a.log"
    --array "1-${regions_count}"
    --time 0:15:00
    --cpus-per-task 4
    --mem 4G
    --job-name generate_ld
    ${code_dir}/04_generate_ld.sh \
        --genotype_prefix ${genotype_prefix} \
        --sample_ids ${reformated_sample_ids} \
        --output_dir ${out_dir} \
        --code_dir ${code_dir} \
        --regions_bed ${regions_to_process_tsv} \
        --regenerate ${regenerate_all} \
)


# Submit on either sherlock or scg
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
        --constraint="nvme" \
        "${sbatch_params[@]}" 
else
    echo "must submit on either 'sherlock' or 'scg'"
fi