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
  NR==1 && $1=="chr" && $2=="start" { print "chr\tstart\tend"; next }
  NR==1 && $1!="chr" { print "chr\tstart\tend" }
  NF>=3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
    chrom=$1; sub(/^chr/,"",chrom);
    from1=$2+1; to1=$3;
    pf=(from1>PAD? from1-PAD:1);
    pt=to1+PAD;
    printf "%s\t%d\t%d\n", chrom, pf, pt;
  }
' "$ld_blocks_bed" > "$padded_ld_blocks_tsv"


regions_count=$(wc -l < "${padded_ld_blocks_tsv}")
echo "Regions to process: ${regions_count}"
if [[ ${regions_count} -gt ${max_array_size} ]]; then
  regions_count=${max_array_size}
fi


sbatch_params=(
    --output "${out_dir}/logs/generate_ld/%A/%A_%a.log"
    --error "${out_dir}/logs/generate_ld/%A/%A_%a.log"
    --array "1-${regions_count}"
    --time 2:00:00
    --cpus-per-task 8
    --mem 50G
    --job-name generate_ld
    ${code_dir}/04_generate_ld.sh \
        --genotype_prefix ${genotype_prefix} \
        --sample_ids ${reformated_sample_ids} \
        --output_dir ${out_dir} \
        --code_dir ${code_dir} \
        --regions_bed ${padded_ld_blocks_tsv} \
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