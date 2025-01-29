#!/usr/bin/env bash

# sample args
participant_id=GTEX-1A3MV
full_vcf_file=/oak/stanford/groups/smontgom/dnachun/data/gtex/data/references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs


export _JAVA_OPTIONS="-Xmx64g"
bash run_participant_vcfs.sh \
    --participant_id ${participant_id} \
    --full_vcf_file ${full_vcf_file} \
    --output_dir ${vcf_dir} \