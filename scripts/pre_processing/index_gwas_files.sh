#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit
module load tabix

sumstat_dir="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/imputed_gwas_hg38_1.1"

echo $(date +"[%b %d %H:%M:%S] Starting bgzip conversion, chr prefix removal, and tabix indexing for all files")

cd ${sumstat_dir}

# Count total files
total_files=$(ls imputed_*.txt.gz 2>/dev/null | wc -l)
echo "Found ${total_files} summary statistics files to process"

if [[ ${total_files} -eq 0 ]]; then
    echo "No summary statistics files found in ${sumstat_dir}"
    exit 1
fi

current_file=0
for sumstat_file in imputed_*.txt.gz; do
    if [[ -f ${sumstat_file} ]]; then
        current_file=$((current_file + 1))
        
        # Check if already indexed
        if [[ -f ${sumstat_file}.tbi ]]; then
            echo $(date +"[%b %d %H:%M:%S] Skipping ${sumstat_file} - already indexed (${current_file}/${total_files})")
            continue
        fi
        
        echo $(date +"[%b %d %H:%M:%S] Processing ${sumstat_file} (${current_file}/${total_files})")
        
        # Get the base name without .gz extension
        base_name=${sumstat_file%.gz}
        
        # Convert from gzip to bgzip, remove chr prefix from chromosome column only, and create index
        gunzip ${sumstat_file}           # Creates .txt file
        
        # Remove chr prefix from chromosome column (column 3) only, preserve header and variant IDs
        awk 'BEGIN{FS=OFS="\t"} NR==1{print} NR>1{gsub(/^chr/, "", $3); print}' ${base_name} > ${base_name}.tmp
        mv ${base_name}.tmp ${base_name}
        
        bgzip ${base_name}               # Creates .txt.gz file (bgzipped)
        tabix -s 3 -b 4 -e 4 -S 1 ${sumstat_file}
        
        if [[ -f ${sumstat_file}.tbi ]]; then
            echo $(date +"[%b %d %H:%M:%S] Successfully processed ${sumstat_file}")
        else
            echo $(date +"[%b %d %H:%M:%S] ERROR: Failed to process ${sumstat_file}")
        fi
    fi
done

echo $(date +"[%b %d %H:%M:%S] Processing completed for ${total_files} files")

# Verify all files are now indexed
missing_indices=0
for file in imputed_*.txt.gz; do
    if [[ ! -f ${file}.tbi ]]; then
        echo "Missing index: ${file}"
        missing_indices=$((missing_indices + 1))
    fi
done

if [[ ${missing_indices} -gt 0 ]]; then
    echo "WARNING: ${missing_indices} files still missing indices"
else
    echo "SUCCESS: All files now have tabix indices"
fi