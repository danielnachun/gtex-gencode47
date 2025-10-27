#!/usr/bin/env python3

"""
combine_colocboost.py

Combine formatted colocboost outputs into aggregated files by category.

- Scans a coloc output directory for files matching patterns and concatenates them.
- Handles both regular and robust versions of colocboost outputs.
- Writes these optional outputs if matching files exist:
  Regular files:
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.v39.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.all.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.individual.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.v39.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.all.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.individual.txt
  Robust files:
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.v39.robust.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.all.robust.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.individual.robust.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.v39.robust.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.all.robust.txt
    <aggregated_output_dir>/<tissue_id>.separate_gwas.individual.robust.txt
"""

import argparse
import glob
import os
import sys

import pandas as pd
from tqdm.auto import tqdm


def read_and_tag_generator(files, desc=None):
    """Yield DataFrames read from TSVs with an added source_file column.

    Skips empty or malformed files gracefully and logs the skipped paths.
    """
    for f in tqdm(files, desc=desc):
        try:
            # Quick size check avoids parsing known-empty files
            if os.path.getsize(f) == 0:
                print(f"[combine_colocboost] Skipping empty file: {f}", file=sys.stderr)
                continue
            df = pd.read_table(f)
        except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError, ValueError) as e:
            print(f"[combine_colocboost] Skipping unreadable file: {f} ({e})", file=sys.stderr)
            continue
        df["source_file"] = f
        yield df


def concat_or_none(gen):
    """Consume a generator of DataFrames; return concatenated DataFrame or None if empty."""
    collected = []
    for df in gen:
        collected.append(df)
    if not collected:
        return None
    return pd.concat(collected, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Combine formatted colocboost outputs")
    parser.add_argument("--coloc_output_dir", required=True, help="Directory containing formatted outputs")
    parser.add_argument("--aggregated_output_dir", default=None, help="Directory to write aggregated outputs")
    parser.add_argument("--tissue_id", default=None, help="Tissue identifier for output filenames")
    parser.add_argument("--quiet", action="store_true", help="Suppress progress messages")
    args = parser.parse_args()

    coloc_output_dir = args.coloc_output_dir
    aggregated_output_dir = args.aggregated_output_dir or f"{coloc_output_dir}/aggregated"
    tissue_id = args.tissue_id or os.path.basename(os.path.normpath(coloc_output_dir))

    os.makedirs(aggregated_output_dir, exist_ok=True)

    if not args.quiet:
        print(f"Scanning dir: {coloc_output_dir}")
        print(f"Aggregated out: {aggregated_output_dir}")
        print(f"Tissue ID: {tissue_id}")

    # All xqtl_coloc text files (both regular and robust)
    coloc_files = glob.glob(f"{coloc_output_dir}/*.xqtl_coloc.txt")
    coloc_robust_files = glob.glob(f"{coloc_output_dir}/*.xqtl_coloc.robust.txt")

    # All separate_gwas text files (both regular and robust)
    separate_gwas_files = glob.glob(f"{coloc_output_dir}/*.separate_gwas.*.txt")
    separate_gwas_robust_files = glob.glob(f"{coloc_output_dir}/*.separate_gwas.*.robust.txt")

    # Helper function to filter out empty files
    def filter_non_empty_files(file_list):
        return [f for f in file_list if os.path.getsize(f) > 0]

    # Regular xqtl_coloc files
    v39_files = filter_non_empty_files([f for f in coloc_files if "v39_genes" in f])
    all_files = filter_non_empty_files([f for f in coloc_files if "all_genes" in f])
    individual_files = filter_non_empty_files(glob.glob(f"{coloc_output_dir}/individual_genes/*.xqtl_coloc.txt"))
    
    # Robust xqtl_coloc files
    v39_robust_files = filter_non_empty_files([f for f in coloc_robust_files if "v39_genes" in f])
    all_robust_files = filter_non_empty_files([f for f in coloc_robust_files if "all_genes" in f])
    individual_robust_files = filter_non_empty_files(glob.glob(f"{coloc_output_dir}/individual_genes/*.xqtl_coloc.robust.txt"))

    # Regular separate_gwas files
    v39_separate_files = filter_non_empty_files([f for f in separate_gwas_files if "v39_genes" in f])
    all_separate_files = filter_non_empty_files([f for f in separate_gwas_files if "all_genes" in f])
    individual_separate_files = filter_non_empty_files(glob.glob(f"{coloc_output_dir}/individual_genes/*.separate_gwas.*.txt"))
    
    # Robust separate_gwas files
    v39_separate_robust_files = filter_non_empty_files([f for f in separate_gwas_robust_files if "v39_genes" in f])
    all_separate_robust_files = filter_non_empty_files([f for f in separate_gwas_robust_files if "all_genes" in f])
    individual_separate_robust_files = filter_non_empty_files(glob.glob(f"{coloc_output_dir}/individual_genes/*.separate_gwas.*.robust.txt"))

    # --- Process both all, then both v39, then both individual ---

    # 1. all_genes (regular and robust)
    if all_files:
        all_gen = read_and_tag_generator(all_files, desc=None if args.quiet else "Reading all genes files")
        all_df = concat_or_none(all_gen)
        if all_df is not None:
            all_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.all.txt"
            all_df.to_csv(all_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {all_out}")
        else:
            print("[combine_colocboost] No valid all-genes files to write.", file=sys.stderr)

    if all_robust_files:
        all_robust_gen = read_and_tag_generator(all_robust_files, desc=None if args.quiet else "Reading all genes robust files")
        all_robust_df = concat_or_none(all_robust_gen)
        if all_robust_df is not None:
            all_robust_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.all.robust.txt"
            all_robust_df.to_csv(all_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {all_robust_out}")
        else:
            print("[combine_colocboost] No valid all-genes robust files to write.", file=sys.stderr)

    # 2. v39_genes (regular and robust)
    if v39_files:
        v39_gen = read_and_tag_generator(v39_files, desc=None if args.quiet else "Reading v39 files")
        v39_df = concat_or_none(v39_gen)
        if v39_df is not None:
            v39_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.v39.txt"
            v39_df.to_csv(v39_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {v39_out}")
        else:
            print("[combine_colocboost] No valid v39 files to write.", file=sys.stderr)

    if v39_robust_files:
        v39_robust_gen = read_and_tag_generator(v39_robust_files, desc=None if args.quiet else "Reading v39 robust files")
        v39_robust_df = concat_or_none(v39_robust_gen)
        if v39_robust_df is not None:
            v39_robust_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.v39.robust.txt"
            v39_robust_df.to_csv(v39_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {v39_robust_out}")
        else:
            print("[combine_colocboost] No valid v39 robust files to write.", file=sys.stderr)

    # 3. individual (regular and robust)
    if individual_files:
        ind_gen = read_and_tag_generator(individual_files, desc=None if args.quiet else "Reading individual gene files")
        ind_df = concat_or_none(ind_gen)
        if ind_df is not None:
            ind_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.individual.txt"
            ind_df.to_csv(ind_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {ind_out}")
        else:
            print("[combine_colocboost] No valid individual files to write.", file=sys.stderr)

    if individual_robust_files:
        ind_robust_gen = read_and_tag_generator(individual_robust_files, desc=None if args.quiet else "Reading individual gene robust files")
        ind_robust_df = concat_or_none(ind_robust_gen)
        if ind_robust_df is not None:
            ind_robust_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.individual.robust.txt"
            ind_robust_df.to_csv(ind_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {ind_robust_out}")
        else:
            print("[combine_colocboost] No valid individual robust files to write.", file=sys.stderr)

    # --- Process separate_gwas files ---

    # 1. all_genes separate_gwas (regular and robust)
    if all_separate_files:
        all_sep_gen = read_and_tag_generator(all_separate_files, desc=None if args.quiet else "Reading all genes separate_gwas files")
        all_sep_df = concat_or_none(all_sep_gen)
        if all_sep_df is not None:
            all_sep_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.all.txt"
            all_sep_df.to_csv(all_sep_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {all_sep_out}")
        else:
            print("[combine_colocboost] No valid all-genes separate_gwas files to write.", file=sys.stderr)

    if all_separate_robust_files:
        all_sep_robust_gen = read_and_tag_generator(all_separate_robust_files, desc=None if args.quiet else "Reading all genes separate_gwas robust files")
        all_sep_robust_df = concat_or_none(all_sep_robust_gen)
        if all_sep_robust_df is not None:
            all_sep_robust_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.all.robust.txt"
            all_sep_robust_df.to_csv(all_sep_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {all_sep_robust_out}")
        else:
            print("[combine_colocboost] No valid all-genes separate_gwas robust files to write.", file=sys.stderr)

    # 2. v39_genes separate_gwas (regular and robust)
    if v39_separate_files:
        v39_sep_gen = read_and_tag_generator(v39_separate_files, desc=None if args.quiet else "Reading v39 separate_gwas files")
        v39_sep_df = concat_or_none(v39_sep_gen)
        if v39_sep_df is not None:
            v39_sep_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.v39.txt"
            v39_sep_df.to_csv(v39_sep_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {v39_sep_out}")
        else:
            print("[combine_colocboost] No valid v39 separate_gwas files to write.", file=sys.stderr)

    if v39_separate_robust_files:
        v39_sep_robust_gen = read_and_tag_generator(v39_separate_robust_files, desc=None if args.quiet else "Reading v39 separate_gwas robust files")
        v39_sep_robust_df = concat_or_none(v39_sep_robust_gen)
        if v39_sep_robust_df is not None:
            v39_sep_robust_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.v39.robust.txt"
            v39_sep_robust_df.to_csv(v39_sep_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {v39_sep_robust_out}")
        else:
            print("[combine_colocboost] No valid v39 separate_gwas robust files to write.", file=sys.stderr)

    # 3. individual separate_gwas (regular and robust)
    if individual_separate_files:
        ind_sep_gen = read_and_tag_generator(individual_separate_files, desc=None if args.quiet else "Reading individual gene separate_gwas files")
        ind_sep_df = concat_or_none(ind_sep_gen)
        if ind_sep_df is not None:
            ind_sep_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.individual.txt"
            ind_sep_df.to_csv(ind_sep_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {ind_sep_out}")
        else:
            print("[combine_colocboost] No valid individual separate_gwas files to write.", file=sys.stderr)

    if individual_separate_robust_files:
        ind_sep_robust_gen = read_and_tag_generator(individual_separate_robust_files, desc=None if args.quiet else "Reading individual gene separate_gwas robust files")
        ind_sep_robust_df = concat_or_none(ind_sep_robust_gen)
        if ind_sep_robust_df is not None:
            ind_sep_robust_out = f"{aggregated_output_dir}/{tissue_id}.separate_gwas.individual.robust.txt"
            ind_sep_robust_df.to_csv(ind_sep_robust_out, sep="\t", index=False)
            if not args.quiet:
                print(f"Wrote {ind_sep_robust_out}")
        else:
            print("[combine_colocboost] No valid individual separate_gwas robust files to write.", file=sys.stderr)

    return 0

def split_robust_results(robust_gwas_df, cos_threshold=0.5, npc_threshold=0.75):
    """Split robust results with low cos_npc and filter by npc_threshold."""
    results_per_tissue = []

    for tissue_id in robust_gwas_df['tissue_id'].unique():
        robust_gwas_df_tissue = robust_gwas_df[robust_gwas_df['tissue_id'] == tissue_id].copy()

        print(f'\n=== Processing {tissue_id} ===')
        print(f'Total robust results: {len(robust_gwas_df_tissue)}')

        # Find low cos_npc results to split
        low_cos = robust_gwas_df_tissue[
            (robust_gwas_df_tissue['cos_npc'] < cos_threshold) &
            (robust_gwas_df_tissue['cs_type'] == 'trait_shared')
        ].copy()

        print(f'Found {len(low_cos)} robust results with cos_npc < {cos_threshold}')
        
        if not low_cos.empty:
            print(f'Converting {len(low_cos)} results to trait_specific_from_shared')

            # Convert to trait_specific_from_shared
            low_cos['cs_type'] = 'trait_specific_from_shared'
            low_cos['converted_from_shared'] = True
            low_cos['original_cos_npc'] = low_cos['cos_npc']
            low_cos['cos_npc'] = np.nan  # Set cos_npc to NaN for trait_specific_from_shared
            low_cos['original_cos_id'] = low_cos['cs_id']

            # Split columns and explode - fix dtype issues by converting to object first
            for col in ['neg_log10_p_value', 'npc_outcome']:
                low_cos[col] = low_cos[col].astype(str).str.split(', ')
            low_cos['phenotype_id'] = low_cos['phenotype_id'].astype(str).str.split(',')

            # Handle cs_id splitting
            cs_id_split = low_cos['cs_id'].astype(str).str.split(':')
            low_cos['cs_id'] = [x[1] if isinstance(x, list) and len(x) > 1 else np.nan for x in cs_id_split]
            low_cos['cs_id'] = low_cos['cs_id'].astype(str).str.split('_')

            # Explode all split columns
            low_cos = low_cos.explode(['phenotype_id', 'npc_outcome', 'neg_log10_p_value', 'cs_id'])
            low_cos['cs_id'] = 'ucos_new:' + low_cos['cs_id'].astype(str)
            low_cos = low_cos.reset_index(drop=True)

            print(f'After exploding: {len(low_cos)} individual results')

        # Get existing trait_specific_from_shared results (they already have original_npc_outcome)
        trait_specific_from_shared = robust_gwas_df_tissue[
            robust_gwas_df_tissue['cs_type'] == 'trait_specific_from_shared'
        ].copy()

        print(f'Found {len(trait_specific_from_shared)} existing trait_specific_from_shared results')

        # Ensure all trait_specific_from_shared have cos_npc = NaN
        if not trait_specific_from_shared.empty:
            trait_specific_from_shared['cos_npc'] = np.nan

        # Combine all split results
        split_robust_results_df = pd.concat([
            df for df in [trait_specific_from_shared, low_cos] if not df.empty
        ], ignore_index=True) if any(not df.empty for df in [trait_specific_from_shared, low_cos]) else pd.DataFrame()

        print(f'Combined split results: {len(split_robust_results_df)} total')

        # Filter by npc_threshold
        if not split_robust_results_df.empty:
            before_filter = len(split_robust_results_df)
            split_robust_results_df = split_robust_results_df[
                pd.to_numeric(split_robust_results_df['original_npc_outcome'], errors='coerce') > npc_threshold
            ]
            after_filter = len(split_robust_results_df)
            print(f'Filtered by npc_threshold > {npc_threshold}: {before_filter} -> {after_filter} results')
        else:
            print(f'No split results to filter')

        # Add high cos_npc results
        keep_mask = (
            (robust_gwas_df_tissue['cos_npc'] > cos_threshold) & (robust_gwas_df_tissue['cs_type'] == 'trait_shared') |
            (robust_gwas_df_tissue['cs_type'] == 'trait_specific')
        )
        high_cos_results = robust_gwas_df_tissue[keep_mask]
        print(f'High cos_npc results to keep: {len(high_cos_results)}')

        if not split_robust_results_df.empty:
            final_results = pd.concat([split_robust_results_df, high_cos_results], ignore_index=True)
            print(f'Final combined results: {len(final_results)} total')
        else:
            final_results = high_cos_results
            print(f'Final results (only high cos_npc): {len(final_results)} total')

        results_per_tissue.append(final_results)

    final_df = pd.concat(results_per_tissue, ignore_index=True) if results_per_tissue else pd.DataFrame()
    print(f'\n=== FINAL SUMMARY ===')
    print(f'Total results across all tissues: {len(final_df)}')
    print(f'Results by tissue: {[len(df) for df in results_per_tissue]}')

    return final_df


if __name__ == "__main__":
    sys.exit(main())

