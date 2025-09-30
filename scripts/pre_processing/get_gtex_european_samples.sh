#!/usr/bin/env bash
set -euo pipefail

PHENOTYPES="/home/klawren/oak/gtex/data/other_references/v10/GTEx_Analysis_2022-06-06_v10_Annotations_SubjectPhenotypesDS.txt"
OUT="/home/klawren/oak/gtex/data/other_references/v10/race3_subject_ids.txt"

awk -F '\t' 'NR==1{
  for(i=1;i<=NF;i++){if($i=="RACE") rc=i; if($i=="SUBJID") sc=i}
  next
}
$rc==3 {print $sc}' "$PHENOTYPES" > "$OUT"

echo "Wrote $(wc -l < "$OUT") subject IDs to $OUT"