#!/bin/bash

# This script processes the Barberia full GWAS metadata file to create a simplified version
OUTPUT_PATH="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/barberia_simplified_gwas_metadata.txt"

awk -F'\t' 'BEGIN{OFS="\t"; print "Tag", "Sample_Size", "Cases", "Controls"} NR>1{if($19=="NA" || $19=="") {cases=0; controls=0} else {cases=$19; controls=$17-$19} print $2, $17, cases, controls}' /home/klawren/oak/gtex/data/pecotmr_references/barberia_full_gwas_metadata.txt > ${OUTPUT_PATH}