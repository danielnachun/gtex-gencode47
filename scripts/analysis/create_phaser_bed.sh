
#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 input.gtf output.bed" >&2
  exit 1
fi

GTF=$1
OUT=$2

awk '
  BEGIN { FS=OFS="\t" }
  # skip header/comment lines; only process gene features
  $0 !~ /^#/ && $3=="gene" {
    # extract gene_id
    if ( match($9, /gene_id "([^"]+)"/, m) ) {
      id = m[1]
      # GTF is 1-based inclusive; BED is 0-based half-open
      start = $4 - 1
      end   = $5
      print $1, start, end, id
    }
  }
' "$GTF" 
| sort -k1,1 -k2,2n \
> "$OUT"
