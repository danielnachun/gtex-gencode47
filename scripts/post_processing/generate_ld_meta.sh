#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

CONFIG_FILE="${1:-/home/klawren/oak/gtex/config/04_ld_gtex_eur.sh}"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Determine output path: prefer OUTPUT_TSV; else derive from out_dir
if [[ -z "${OUTPUT_TSV:-}" ]]; then
  : "${out_dir:?Error: set OUTPUT_TSV or out_dir in the config}"
  OUTPUT_TSV="${out_dir%/}/ld_metadata.tsv"
fi

output_path="${out_dir%/}/ld_metadata.tsv"
# clear the file if it exists
if [[ -f "$output_path" ]]; then
  echo "Clearing $output_path"
  > "$output_path"
fi

padded_ld_blocks_path="${out_dir%/}/padded_ld_blocks.tsv"

# Build the metadata TSV
# - Accepts header with at least the first three columns = chr, start, end (any names)
# - Normalizes chromosome to include 'chr' prefix for file naming
# - Constructs path: <LD_DIR>/LD_chr<chrom>_<start>_<end>.ld.gz,<LD_DIR>/LD_chr<chrom>_<start>_<end>.bim
# - Skips rows where either the LD file or bim file is missing (warns to stderr)
awk -v dir="$out_dir" 'BEGIN{
  OFS = "\t";
}
BEGIN{
  # Print standardized header
  print "chrom","start","end","path";
}
{
  chr_raw = $1;
  start = $2;
  end = $3;
  if (chr_raw == "" || start == "" || end == "") next;
  
  ld_path = dir "/LD_chr" chr_raw "_" start "_" end ".ld.gz";
  bim_path = dir "/LD_chr" chr_raw "_" start "_" end ".bim";
  
  # Check if both files exist
  ld_cmd = "[ -f \"" ld_path "\" ]";
  bim_cmd = "[ -f \"" bim_path "\" ]";
  
  if (system(ld_cmd) == 0 && system(bim_cmd) == 0) {
    # Both files exist, create comma-separated path
    print chr_raw, start, end, ld_path "," bim_path;
  } else {
    # At least one file is missing, print warning
    if (system(ld_cmd) != 0) {
      printf "Warning: missing LD file for %s:%s-%s at %s\n", chr_raw, start, end, ld_path > "/dev/stderr";
    }
    if (system(bim_cmd) != 0) {
      printf "Warning: missing bim file for %s:%s-%s at %s\n", chr_raw, start, end, bim_path > "/dev/stderr";
    }
  }
}' "$padded_ld_blocks_path" > "$output_path"

echo "$output_path"