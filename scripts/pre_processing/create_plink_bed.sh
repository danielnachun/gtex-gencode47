#!/bin/bash

# load plink
module load plink

# paths
VCF_FILE="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/realign_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"
OUTPUT_BASENAME="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/pecotmr_references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.MAF01"

# create plink bed file
plink --vcf ${VCF_FILE} --maf 0.01 --make-bed --out ${OUTPUT_BASENAME}
