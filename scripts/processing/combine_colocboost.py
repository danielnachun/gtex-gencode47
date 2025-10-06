#!/usr/bin/env python3

"""
combine_colocboost.py

Combine formatted colocboost outputs into aggregated files by category.

- Scans a coloc output directory for files matching patterns and concatenates them.
- Writes these optional outputs if matching files exist:
  <aggregated_output_dir>/<tissue_id>.xqtl_coloc.v39.txt
  <aggregated_output_dir>/<tissue_id>.xqtl_coloc.all.txt
  <aggregated_output_dir>/<tissue_id>.xqtl_coloc.individual.txt
"""

import argparse
import glob
import os
import sys

import pandas as pd
from tqdm.auto import tqdm


def read_and_tag_generator(files, desc=None):
    """Yield DataFrames read from TSVs with an added source_file column."""
    for f in tqdm(files, desc=desc):
        df = pd.read_table(f)
        df["source_file"] = f
        yield df


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

    # All xqtl_coloc text files
    coloc_files = glob.glob(f"{coloc_output_dir}/*.xqtl_coloc.txt")

    v39_files = [f for f in coloc_files if "v39_genes" in f]
    all_files = [f for f in coloc_files if "all_genes" in f]
    individual_files = glob.glob(f"{coloc_output_dir}/individual_genes/*.xqtl_coloc.txt")

    if v39_files:
        v39_df = pd.concat(read_and_tag_generator(v39_files, desc=None if args.quiet else "Reading v39 files"), ignore_index=True)
        v39_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.v39.txt"
        v39_df.to_csv(v39_out, sep="\t", index=False)
        if not args.quiet:
            print(f"Wrote {v39_out}")

    if all_files:
        all_df = pd.concat(read_and_tag_generator(all_files, desc=None if args.quiet else "Reading all genes files"), ignore_index=True)
        all_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.all.txt"
        all_df.to_csv(all_out, sep="\t", index=False)
        if not args.quiet:
            print(f"Wrote {all_out}")

    if individual_files:
        ind_df = pd.concat(read_and_tag_generator(individual_files, desc=None if args.quiet else "Reading individual gene files"), ignore_index=True)
        ind_out = f"{aggregated_output_dir}/{tissue_id}.xqtl_coloc.individual.txt"
        ind_df.to_csv(ind_out, sep="\t", index=False)
        if not args.quiet:
            print(f"Wrote {ind_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

