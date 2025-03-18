#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# sample args
full_vcf=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
output_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs
participant_id_list=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/test_participants.txt
code_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow/code


mkdir -p ${output_dir}
mkdir -p ${output_dir}/logs

participant_id_length=$(wc -l < ${participant_id_list})
echo "BAM array length: ${participant_id_length}"


sbatch --output "${output_dir}/logs/%A_%a.log" \
    --error "${output_dir}/logs/%A_%a.log" \
    --array "1-${participant_id_length}%250" \
    --time 4:00:00 \
    --partition normal,owners \
    --cpus-per-task 1 \
    --mem 64G \
    --job-name participant_vcf \
    ${code_dir}/extract_participant_vcfs.sh \
        --participant_id_list ${participant_id_list} \
        --full_vcf ${full_vcf} \
        --output_dir ${output_dir}
