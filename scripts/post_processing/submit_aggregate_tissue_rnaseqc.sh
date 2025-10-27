#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

# Configuration
OUTPUT_DIR="/home/klawren/oak/gtex/output/all_tissues_quantifications"
SAMPLE_META_FILE="/home/klawren/oak/gtex/data/other_references/v10/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
TISSUE_LIST_FILE="/home/klawren/oak/gtex/data/pecotmr_references/qtl_tissue_ids.txt"
CODE_DIR="/home/klawren/oak/gtex/scripts/post_processing"
LOG_DIR="${OUTPUT_DIR}/logs/tissue_aggregation"

# Create directories
mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}/rnaseqc_agg"

# Extract tissue list from sample metadata using Python
echo "Extracting tissue list from sample metadata..."

python3 -c "
import pandas as pd

# Load sample metadata
sample_meta = pd.read_csv('${SAMPLE_META_FILE}', sep='\t')
passed_samples = sample_meta[sample_meta['SMAFRZE']=='RNASEQ']

# Apply tissue cutoff filter
tissue_cutoff = 30
large_sample_size_tissues = passed_samples.groupby('SMTSD').size()[passed_samples.groupby('SMTSD').size() > tissue_cutoff].index
print('continuing with {} tissues'.format(len(large_sample_size_tissues)))

# Get tissue sample lists for large sample size tissues only
tissue_sample_lists = passed_samples.groupby('SMTSD').agg({'SAMPID':'unique'})
tissue_sample_lists = tissue_sample_lists[tissue_sample_lists.index.isin(large_sample_size_tissues)]

# Write tissue list to file
with open('${TISSUE_LIST_FILE}', 'w') as f:
    for tissue in tissue_sample_lists.index:
        f.write(f'{tissue}\n')

print(f'Extracted {len(tissue_sample_lists)} tissues with >30 samples to ${TISSUE_LIST_FILE}')
"

# Count tissues
num_tissues=$(wc -l < "${TISSUE_LIST_FILE}")
echo "Total tissues to process: ${num_tissues}"

# Submit sbatch job arrays
echo "Submitting tissue aggregation job arrays..."

# Create array of tissue IDs for the job array
tissue_array_file="${TISSUE_LIST_FILE}.array"
cp "${TISSUE_LIST_FILE}" "${tissue_array_file}"

# Calculate total number of jobs (2 per tissue: tpm and reads)
total_jobs=$((num_tissues * 2))
echo "Total jobs to submit: ${total_jobs} (${num_tissues} tissues × 2 file types)"

# Submit job array for gene_tpm
echo "Submitting gene_tpm job array..."
job_id_tpm=$(sbatch \
    --account smontgom \
    --output "${LOG_DIR}/%A_%a_tpm.log" \
    --error "${LOG_DIR}/%A_%a_tpm.log" \
    --array "1-${num_tissues}" \
    --time 4:00:00 \
    --cpus-per-task 8 \
    --mem 32G \
    --job-name "agg_tissue_tpm" \
    "${CODE_DIR}/aggregate_tissue_rnaseqc.py" \
        --tissue_list_file "${tissue_array_file}" \
        --output_dir "${OUTPUT_DIR}" \
        --sample_meta_file "${SAMPLE_META_FILE}" \
        --file_end "gene_tpm" \
        --tissue_cutoff 0 | awk '{print $4}')

# Submit job array for gene_reads  
echo "Submitting gene_reads job array..."
job_id_reads=$(sbatch \
    --account smontgom \
    --output "${LOG_DIR}/%A_%a_reads.log" \
    --error "${LOG_DIR}/%A_%a_reads.log" \
    --array "1-${num_tissues}" \
    --time 4:00:00 \
    --cpus-per-task 8 \
    --mem 32G \
    --job-name "agg_tissue_reads" \
    "${CODE_DIR}/aggregate_tissue_rnaseqc.py" \
        --tissue_list_file "${tissue_array_file}" \
        --output_dir "${OUTPUT_DIR}" \
        --sample_meta_file "${SAMPLE_META_FILE}" \
        --file_end "gene_reads" \
        --tissue_cutoff 0 | awk '{print $4}')

echo "Submitted job array ${job_id_tpm} for gene_tpm processing"
echo "Submitted job array ${job_id_reads} for gene_reads processing"
echo "All tissue aggregation job arrays submitted!"
echo "Monitor jobs with: squeue -u \$USER"
echo "Check logs in: ${LOG_DIR}"
