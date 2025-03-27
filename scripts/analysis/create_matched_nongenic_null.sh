#!/bin/bash
module load bedtools

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

# Input parameters
GENE_GTF="/home/klawren/oak/gtex/data/realign_references/gencode.v47.genes.gtf"
ALL_GTF="/home/klawren/oak/gtex/data/realign_references/gencode.v47.annotation.gtf"
CHROM_SIZES="/home/klawren/oak/gtex/data/realign_references/GRCh38.chrsizes"
TMP_DIR="/home/klawren/oak/gtex/data/other_references/nongenic_null_with_introns"

# pull out exons
awk -v OFS='\t' '$3=="exon" {
    # extract exon_id and gene_id
    match($0, /exon_id "([^"]+)"/, eid)
    match($0, /gene_id "([^"]+)"/, gid)
    print $1, $4-1, $5, eid[1], gid[1]
}' "${ALL_GTF}" | \
    bedtools sort -i - > "${TMP_DIR}/all_exons.bed"

# pull out exons
awk -v OFS='\t' '$3=="exon" {
    # extract exon_id and gene_id
    match($0, /exon_id "([^"]+)"/, eid)
    match($0, /gene_id "([^"]+)"/, gid)
    print $1, $4-1, $5, eid[1], gid[1]
}' "${GENE_GTF}" | \
    bedtools sort -i - > "${TMP_DIR}/exons_to_match.bed"

# pull out genes 
awk -v OFS='\t' '$3=="gene" {print $1, $4-1, $5}' "${ALL_GTF}" | \
    bedtools sort -i - | \
    bedtools merge -i - > "${TMP_DIR}/merged_genes.bed"

# create intronic regions (genes minus exons + 100)
bedtools slop -i "${TMP_DIR}/all_exons.bed" -g "${CHROM_SIZES}" -b 100 > "${TMP_DIR}/exons_buffer.bed"
bedtools subtract -a "${TMP_DIR}/merged_genes.bed" \
    -b "${TMP_DIR}/exons_buffer.bed" > "${TMP_DIR}/intronic_regions.bed"

# create intergenic regions (genome minus genes + 1000)
bedtools slop -i "${TMP_DIR}/merged_genes.bed" -g "${CHROM_SIZES}" -b 100 > "${TMP_DIR}/genes_buffer.bed"
awk -v OFS='\t' '{print $1, "0", $2}' "${CHROM_SIZES}" | \
    bedtools subtract -a - -b "${TMP_DIR}/genes_buffer.bed" > "${TMP_DIR}/intergenic_regions.bed"

# combine intronic and intergenic regions
cat "${TMP_DIR}/intergenic_regions.bed" "${TMP_DIR}/intronic_regions.bed" | \
    bedtools sort -i - > "${TMP_DIR}/background_regions.bed"


# this works but is too slow, python implementiaion in notebooks/matched_nongenic_null.ipynb
# # find the nearest region that can accomedate each exon 
# > "${TMP_DIR}/nongenic_exon_matches.bed"
# TOTAL_EXONS=$(wc -l < "${TMP_DIR}/all_exons.bed")

# find_nearest_non_genic_region() {
#     local exon="$1"
#     local exon_chr=$(echo "$exon" | cut -f1)
#     local exon_start=$(echo "$exon" | cut -f2)
#     local exon_end=$(echo "$exon" | cut -f3)
#     local exon_id=$(echo "$exon" | cut -f4)
#     local gene_id=$(echo "$exon" | cut -f5)
#     local exon_length=$((exon_end - exon_start))

#     local nearest_region=$(awk -v chr="$exon_chr" -v start="$exon_start" -v exon_length="$exon_length" '
#         $1 == chr && ($3 - $2) >= exon_length {
#             distance = ($2 > start) ? $2 - start : start - $2
#             if (min_distance == "" || distance < min_distance) {
#                 min_distance = distance
#                 nearest_region = $0
#             }
#         }
#         END { print nearest_region }
#     ' "${TMP_DIR}/background_regions.bed")

#     if [ -n "$nearest_region" ]; then
#         local region_start=$(echo "$nearest_region" | cut -f2)
#         local region_end=$((region_start + exon_length))
#         echo -e "$exon_chr\t$region_start\t$region_end\t${exon_id}_null\t${gene_id}_null"
#     fi
# }

# # Use pv to show progress while reading the input file
# pv -l -s "$TOTAL_EXONS" "${TMP_DIR}/all_exons.bed" | while read -r exon; do
#     find_nearest_non_genic_region "$exon" >> "${TMP_DIR}/nongenic_exon_matches.bed"
# done


