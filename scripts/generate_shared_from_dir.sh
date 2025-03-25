#!/usr/bin/env bash

vcf_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs"
bam_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/raw/GTEx_Analysis_2022-06-06_v10_RNAseq_BAM_files"
output_file="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/all_tissues/all_shared_samples.txt"

# clear output file
> "$output_file"

# check each sample with expression data
for bam_file in "$bam_dir"/*.Aligned.sortedByCoord.out.patched.md.bam; do
    # extract the participant ID and sample ID 
    participant_id=$(basename "$bam_file" | cut -d '-' -f 1-2) 
    sample_id=$(basename "$bam_file" .Aligned.sortedByCoord.out.patched.md.bam)
    # each sample needs to have genotype data (in the vcf dir)
    vcf_file="$vcf_dir/${participant_id}.snps.het.vcf.gz"
    if [[ -f "$vcf_file" ]]; then
        if [[ -n "$sample_id" ]]; then
            echo "$sample_id" >> "$output_file"
        fi
    fi
done

echo "Matching IDs have been written to $output_file."