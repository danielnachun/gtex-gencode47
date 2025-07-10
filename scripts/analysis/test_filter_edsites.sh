#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit


# create a merged vcf of dbsnp and gtex wgs sites, which we filter out
bcftools merge /home/klawren/oak/gtex/data/edsite_references/All_20180418.vcf.gz \
        /home/klawren/oak/gtex/data/edsite_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz \
        -o /home/klawren/oak/gtex/data/edsite_references/combined_common_vars.vcf.gz
echo "done"
gatk IndexFeatureFile \
     -I /home/klawren/oak/gtex/data/edsite_references/combined_common_vars.vcf.gz

# index the excluded regions
gatk IndexFeatureFile \
     -I /home/klawren/oak/gtex/data/edsite_references/ENCFF356LFX.bed

# sort and index splice sites
# splice sites created from gencode with `notebooks/create_edsite_regions.ipynb`
sort -k1,1 -k2,2n /home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.bed >/home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.sorted.bed
gatk IndexFeatureFile \
     -I /home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.sorted.bed 


# fitler the mutect samples
gatk FilterMutectCalls \
   -R /home/klawren/oak/gtex/data/edsite_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   -V /home/klawren/oak/gtex/output/test_bams/output_kate/mutect_first/GTEX-1C4CL-2126-SM-7IGQC.mutect2.vcf \
   -O /home/klawren/oak/gtex/output/test_bams/output_kate/mutect_first/GTEX-1C4CL-2126-SM-7IGQC.mutect2_filtered.vcf \
   --min-median-base-quality 20 \
   --max-alt-allele-count 1 \
   --min-median-read-position 6 \
   --unique-alt-read-count 3 \
   --read-filter NotSupplementaryAlignmentReadFilter \
   --exclude-intervals /home/klawren/oak/gtex/data/edsite_references/ENCFF356LFX.bed \
   --exclude-intervals /home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.sorted.bed \
   --exclude-intervals /home/klawren/oak/gtex/data/edsite_references/combined_common_vars.vcf.gz
#   --ob-priors /home/klawren/oak/gtex/output/test_bams/output_kate/mutect_first/GTEX-1C4CL-2126-SM-7IGQC.f1r2.tar.gz  # this gives an error



# filter the haplotypecaller samples
bcftools view /home/klawren/oak/gtex/output/test_bams/output_kate/haplotype_caller_first/GTEX-1C4CL-2126-SM-7IGQC.hc.vcf.gz \
    --include "INFO/DP>=10 & MQ>=40 & MQRankSum>=-12.5 & QD>=2 & ReadPosRankSum>=-8 & ReadPosRankSum<8 & AD > 3" \
    --targets-file "^/home/klawren/oak/gtex/data/edsite_references/ENCFF356LFX.bed.gz" \
    --targets-file "^/home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.bed.gz" \
    --min-alleles 2 \
    --max-alleles 2 \
    --targets-file "^/home/klawren/oak/gtex/data/edsite_references/combined_common_vars.vcf.gz" \
    --output /home/klawren/oak/gtex/output/test_bams/output_kate/haplotype_caller/GTEX-1C4CL-2126-SM-7IGQC.hc_filtered.vcf.gz



# alternate approach (not working)
# filter the haplotypecaller samples
gatk VariantFiltration \
   -R /home/klawren/oak/gtex/data/edsite_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
   -V /home/klawren/oak/gtex/output/test_bams/output_kate/haplotype_caller_first/GTEX-1C4CL-2126-SM-7IGQC.hc.vcf.gz \
   -O /home/klawren/oak/gtex/output/test_bams/output_kate/haplotype_caller_first/GTEX-1C4CL-2126-SM-7IGQC.hc_filtered.vcf.gz \
   --read-filter NotSupplementaryAlignmentReadFilter \
   --exclude-intervals /home/klawren/oak/gtex/data/edsite_references/ENCFF356LFX.bed \
   --exclude-intervals /home/klawren/oak/gtex/data/edsite_references/gencode.v47.splice_site.sorted.bed \
    --filter-name "LowCoverage" --filter-expression "DP < 10" \
    --filter-name "LowMappingQuality" --filter-expression "MQ < 40" \
    --filter-name "LowMappingQualityRankSum" --filter-expression "MQRankSum < -12.5" \
    --filter-name "LowQualityDepth" --filter-expression "QD < 2" \
    --filter-name "ExtremeReadPositionRankSum" --filter-expression "ReadPosRankSum > -8 && ReadPosRankSum < 8" \
    --filter-name "MultipleAlternates" --filter-expression "N_ALT > 1" \
    --filter-name "LowAlternateAlleleReads" --filter-expression "AD[1] < 3" 