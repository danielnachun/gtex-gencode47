#!/bin/bash
#SBATCH --job-name=format_colocboosts
#SBATCH --account=smontgom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --time=12:00:00
#SBATCH --output=format_colocboosts_%j.out
#SBATCH --error=format_colocboosts_%j.err

#!/usr/bin/env bash
set -euo pipefail

# Wrapper to run `format_colocboost.R` on every .rds file under the given root (recursively).
# Adjust COLOC_DIR if needed.
COLOC_DIR="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/coloc_qtl_only"
SCRIPT="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing/format_colocboost.R"

source <(pixi shell-hook --environment pecotmr-dev --manifest-path /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing/pixi.toml)

count=0
while IFS= read -r -d '' file; do
  echo "Processing: $file"
  out_prefix="${file%.rds}"
  Rscript "$SCRIPT" --pecotmr_colocboost "$file" --output_prefix "$out_prefix"
  ((count++)) || true
done < <(find "$COLOC_DIR" -type f -name '*.rds' -print0)

echo "Done. Processed $count file(s)."