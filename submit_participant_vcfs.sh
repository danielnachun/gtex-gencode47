#!/usr/bin/env bash

# sample args
full_vcf_file=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs

parallel -j16 "bash -c 'bash ./run_participant_vcfs.sh \
    --participant_id {} \
    --full_vcf_file ${full_vcf_file} \
    --output_dir ${vcf_dir}'" < ./test_participants
