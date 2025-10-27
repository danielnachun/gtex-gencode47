#!/usr/bin/env python3
"""
Process tissue-level RNA-seq aggregation for a single tissue using SLURM job arrays.
This script is designed to be run as part of a job array where each task processes one tissue.
"""

import pandas as pd
import os
import sys
import argparse
from tqdm import tqdm

def read_sample_file(sample_id, output_dir, file_end):
    """Read a single sample file and return the data."""
    expected_filename = f"{sample_id}.{file_end}.gct.gz"
    file_path = os.path.join(os.path.join(output_dir, 'rnaseqc'), expected_filename)
    if os.path.isfile(file_path):
        try:
            return pd.read_csv(file_path, sep='\t', skiprows=2).set_index(['Name', 'Description'])
        except Exception as e:
            print(f"Failed reading {file_path}: {e}")
            return None
    else:
        print(f"File not found: {file_path}")
        return None

def agg_rnaseqc(output_dir, tissue_sample_ids, file_end='gene_reads'):
    """Aggregate RNA-seq quantification data for a tissue."""
    sample_tpms = []
    for sample_id in tqdm(tissue_sample_ids, desc="Processing samples"):
        result = read_sample_file(sample_id, output_dir, file_end)
        if result is not None:
            sample_tpms.append(result)
    if sample_tpms:
        tpm_agg = pd.concat(sample_tpms, axis=1).reset_index()
        return tpm_agg
    else:
        raise FileNotFoundError("No sample files were found for aggregation.")

def get_aggregation_tissue(output_dir, tissue_sample_lists, tissue_id, file_end='gene_reads'):
    """Process aggregation for a single tissue."""
    print(f'Processing samples for {tissue_id}')
    
    # Use the tissue_id directly since it's now in GTEx format
    sample_ids = tissue_sample_lists.loc[tissue_id]['SAMPID']
    tissue_id_clean = tissue_id.replace(' - ', '_').replace('-', '_').replace(' ', '_').replace('(', '').replace(')', '')

    gene_agg_path = os.path.join(os.path.join(output_dir, 'rnaseqc_agg'), f"{tissue_id_clean}.v11.{file_end}.gct.gz")
    print(f'\t{gene_agg_path}')
    
    # Check if file already exists
    if os.path.exists(gene_agg_path):
        print(f'Aggregated file already exists: {gene_agg_path}')
        return gene_agg_path
    
    try:
        gene_agg = agg_rnaseqc(output_dir, sample_ids, file_end=file_end)
        os.makedirs(os.path.join(output_dir, 'rnaseqc_agg'), exist_ok=True)
        gene_agg.to_csv(gene_agg_path, sep='\t', index=False, compression='gzip')
        print(f'Successfully created aggregated file: {gene_agg_path}')
        return gene_agg_path
    except Exception as e:
        print(f'Error processing tissue {tissue_id}: {e}')
        raise

def main():
    parser = argparse.ArgumentParser(description='Process tissue-level RNA-seq aggregation using job arrays')
    parser.add_argument('--tissue_list_file', required=True, help='Path to tissue list file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--sample_meta_file', required=True, help='Path to sample metadata file')
    parser.add_argument('--file_end', default='gene_tpm', help='File ending (default: gene_tpm)')
    parser.add_argument('--tissue_cutoff', type=int, default=0, help='Minimum sample size cutoff')
    
    args = parser.parse_args()
    
    # Get the SLURM array task ID
    slurm_array_task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if not slurm_array_task_id:
        print("Error: SLURM_ARRAY_TASK_ID not found. This script must be run as part of a SLURM job array.")
        sys.exit(1)
    
    task_id = int(slurm_array_task_id)
    print(f"Processing array task {task_id}")
    
    # Load tissue list
    with open(args.tissue_list_file, 'r') as f:
        tissue_list = [line.strip() for line in f if line.strip()]
    
    if task_id < 1 or task_id > len(tissue_list):
        print(f"Error: Task ID {task_id} is out of range. Tissue list has {len(tissue_list)} entries.")
        sys.exit(1)
    
    # Get the tissue for this task
    tissue_id = tissue_list[task_id - 1]  # Convert to 0-based indexing
    print(f"Processing tissue: {tissue_id}")
    
    # Load sample metadata
    print(f'Loading sample metadata from {args.sample_meta_file}')
    sample_meta = pd.read_csv(args.sample_meta_file, sep='\t')
    passed_samples = sample_meta[sample_meta['SMAFRZE'] == 'RNASEQ']
    
    # Apply tissue cutoff filter
    tissue_cutoff = 30
    large_sample_size_tissues = passed_samples.groupby('SMTSD').size()[passed_samples.groupby('SMTSD').size() > tissue_cutoff].index
    print('continuing with {} tissues'.format(len(large_sample_size_tissues)))
    
    # Get tissue sample lists for large sample size tissues only
    tissue_sample_lists = passed_samples.groupby('SMTSD').agg({'SAMPID': 'unique'})
    tissue_sample_lists = tissue_sample_lists[tissue_sample_lists.index.isin(large_sample_size_tissues)]
    
    # Check if tissue exists
    if tissue_id not in tissue_sample_lists.index:
        print(f'Error: Tissue {tissue_id} not found in tissue list')
        print(f'Available tissues: {list(tissue_sample_lists.index)}')
        sys.exit(1)
    
    # Process the tissue
    try:
        result_path = get_aggregation_tissue(
            args.output_dir, 
            tissue_sample_lists, 
            tissue_id, 
            file_end=args.file_end
        )
        print(f'Successfully processed tissue {tissue_id}')
        print(f'Output file: {result_path}')
    except Exception as e:
        print(f'Failed to process tissue {tissue_id}: {e}')
        sys.exit(1)

if __name__ == '__main__':
    main()
