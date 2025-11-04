#!/bin/bash
module load bedtools

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

# Input parameters
GENE_GTF="/oak/stanford/groups/smontgom/dnachun/data/dgtex/dGTEx_Consortium_resources/NHP-dGTEx_references/Enhanced_GTFs/Marmosets/v1/Marmoset_mCalJa1.2.pat.X_enhanced_v1_genes_collapsed_stranded.gtf"
ALL_GTF="/oak/stanford/groups/smontgom/dnachun/data/dgtex/dGTEx_Consortium_resources/NHP-dGTEx_references/Enhanced_GTFs/Marmosets/v1/Marmoset_mCalJa1.2.pat.X_enhanced_v1.gtf"
CHROM_SIZES="/oak/stanford/groups/smontgom/dnachun/data/dgtex/dGTEx_Consortium_resources/NHP-dGTEx_references/FASTA_files/Marmosets/GCA_011100555.2_mCalJa1.2.pat.X_chrM_ucsc_chrs.chrsizes"
TMP_DIR="/oak/stanford/groups/smontgom/klawren/nph_gtex/data/tmp/marmoset"

# pull out all exons
awk -v OFS='\t' '$3=="exon" {
    # extract exon_id and gene_id
    match($0, /exon_id "([^"]+)"/, eid)
    match($0, /gene_id "([^"]+)"/, gid)
    print $1, $4-1, $5, eid[1], gid[1]
}' "${ALL_GTF}" | \
    bedtools sort -i - > "${TMP_DIR}/all_exons.bed"

# pull out exons from the non-overlapping gene file used in rnqseqc
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


# find the nearest region that can accomedate each exon 
# notebooks/matched_nongenic_null.ipynb



