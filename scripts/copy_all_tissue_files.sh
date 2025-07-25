#!/bin/bash

module load google-cloud-cli

echo "Copying files..."
{
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/rnaseq_qc gs://fa-exchange2/all_tissues/rnaseq_qc
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/coverage gs://fa-exchange2/all_tissues/coverage
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/flagstat gs://fa-exchange2/all_tissues/flagstat
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/gatk gs://fa-exchange2/all_tissues/gatk
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/leafcutter gs://fa-exchange2/all_tissues/leafcutter
    gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues/rsem gs://fa-exchange2/all_tissues/rsem
} > copy_all_tissues.log 2>&1

echo "Done."