#!/usr/bin/env python3
"""
Create a modified GWAS metadata file with specific columns.
Extracts Tag, Sample_Size, Cases, and calculates Controls.
"""

import pandas as pd
import sys

def create_modified_gwas_metadata(input_file, output_file):
    """
    Create a modified GWAS metadata file with Tag, Sample_Size, Cases, and Controls columns.
    
    Args:
        input_file (str): Path to the input metadata file
        output_file (str): Path to the output modified metadata file
    """
    # Read the input file
    df = pd.read_csv(input_file, sep='\t')
    
    # Create the modified dataframe with required columns
    modified_df = pd.DataFrame()
    modified_df['Tag'] = df['Tag']
    modified_df['Sample_Size'] = df['Sample_Size']
    modified_df['Cases'] = df['Cases']
    
    # Calculate Controls: if Cases is not NA, then Controls = Sample_Size - Cases
    # If Cases is NA, then both Cases and Controls should be 0
    modified_df['Controls'] = modified_df.apply(
        lambda row: (row['Sample_Size'] - row['Cases']) if pd.notna(row['Cases']) else 0, 
        axis=1
    )
    
    # Set Cases to 0 where it was originally NA
    modified_df['Cases'] = modified_df['Cases'].fillna(0)
    
    # Save the modified file
    modified_df.to_csv(output_file, sep='\t', index=False)
    print(f"Modified GWAS metadata saved to: {output_file}")
    print(f"Total rows: {len(modified_df)}")
    print(f"Columns: {list(modified_df.columns)}")

if __name__ == "__main__":
    input_file = "/home/klawren/oak/gtex/data/pecotmr_references/barberia_full_gwas_metadata.txt"
    output_file = "/home/klawren/oak/gtex/data/pecotmr_references/barberia_modified_gwas_metadata.txt"
    
    create_modified_gwas_metadata(input_file, output_file)
