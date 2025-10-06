#!/bin/bash

# This script processes the Barberia full GWAS metadata file to create a simplified version
# If cases is not 0/NA, set sample size to 0, controls = sample size - cases
# If cases is NA or empty, sample size is original, cases and controls are 0
OUTPUT_PATH="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/barberia_simplified_gwas_metadata.txt"
META_PATH="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/barbeira_gtex_imputed/barberia_full_gwas_metadata.txt"

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "Tag", "Sample_Size", "Cases", "Controls"
}
NR > 1 {
    if ($19 == "NA" || $19 == "") {
        sample_size = $17;
        cases = 0;
        controls = 0;
    } else {
        sample_size = 0;
        cases = $19;
        controls = $17 - $19;
    }
    print $2, sample_size, cases, controls
}
' ${META_PATH} > ${OUTPUT_PATH}