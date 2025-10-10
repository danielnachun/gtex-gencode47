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
  Robust files:
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.v39.robust.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.all.robust.txt
    <aggregated_output_dir>/<tissue_id>.xqtl_coloc.individual.robust.txt
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

    # Regular files
    v39_files = [f for f in coloc_files if "v39_genes" in f]
    all_files = [f for f in coloc_files if "all_genes" in f]
    individual_files = glob.glob(f"{coloc_output_dir}/individual_genes/*.xqtl_coloc.txt")
    
    # Robust files
    v39_robust_files = [f for f in coloc_robust_files if "v39_genes" in f]
    all_robust_files = [f for f in coloc_robust_files if "all_genes" in f]
    individual_robust_files = glob.glob(f"{coloc_output_dir}/individual_genes/*.xqtl_coloc.robust.txt")

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

    return 0


if __name__ == "__main__":
    sys.exit(main())

