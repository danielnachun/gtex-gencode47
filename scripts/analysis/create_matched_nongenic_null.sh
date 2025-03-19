#!/bin/bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

INPUT_GTF="/home/klawren/oak/gtex/data/references/gencode.v47.annotation.gtf"
OUTPUT_GTF="gencode.v47.background.gtf"
CHROM_SIZES="/home/klawren/oak/gtex/data/references/GRCh38.chrsizes"

# 1. Extract exons and genes, and sort them properly
grep -w "exon" $INPUT_GTF | sort -k1,1V -k4,4n > exons.tmp.gtf
grep -w "gene" $INPUT_GTF | sort -k1,1V -k4,4n > genes.tmp.gtf

# 2. Create BED format for exons
awk -v OFS="\t" '{print $1,$4-1,$5,$9,".",$7}' exons.tmp.gtf | sort -k1,1V -k2,2n > exons.tmp.bed

# 3. Create gene BED format
awk -v OFS="\t" '{print $1,$4-1,$5,$9,".",$7}' genes.tmp.gtf | sort -k1,1V -k2,2n > genes.tmp.bed

# 4. Find intronic regions (regions within genes but not in exons)
bedtools subtract -a genes.tmp.bed -b exons.tmp.bed | \
  awk -v OFS="\t" '{
    # Trim 100bp from both sides of intronic regions
    start=$2+100;
    end=$3-100;
    if(end > start) print $1,start,end,$4,".",$6,"intronic"
  }' | sort -k1,1V -k2,2n > intronic.tmp.bed

# 5. Find intergenic regions using the provided chromosome sizes
bedtools complement -i genes.tmp.bed -g $CHROM_SIZES | \
  awk -v OFS="\t" '{
    # Trim 1000bp from both sides of intergenic regions
    start=$2+1000;
    end=$3-1000;
    if(end > start) print $1,start,end,"intergenic",".",".","intergenic"
  }' | sort -k1,1V -k2,2n > intergenic.tmp.bed

# 6. Combine intronic and intergenic regions
cat intronic.tmp.bed intergenic.tmp.bed | sort -k1,1V -k2,2n > background_regions.tmp.bed

# 7. Find nearest background region for each exon and calculate shift
bedtools closest -a exons.tmp.bed -b background_regions.tmp.bed -D ref | \
  awk -v OFS="\t" '{
    shift=$8-$2;  # Calculate shift distance
    type=$13;     # Get the type (intronic or intergenic)
    print $1,$2+shift,$3+shift,$4,".",$6,type
  }' > shifted_exons.tmp.bed

# 8. Convert back to GTF format and add _background to IDs
awk -v OFS="\t" '{
    split($4,attrs,";");
    new_attrs="";
    for(i in attrs) {
        if(attrs[i] ~ /gene_id|transcript_id/) {
            gsub(/"[^"]*"/, "&_background", attrs[i]);
        }
        new_attrs = new_attrs attrs[i] ";";
    }
    source = ($7 == "intronic") ? "intronic_background" : "intergenic_background";
    print $1,source,"exon",$2+1,$3,".",$6,".",new_attrs;
}' shifted_exons.tmp.bed > $OUTPUT_GTF

# 9. Clean up temporary files
rm *.tmp.*

echo "Created background GTF file: $OUTPUT_GTF"