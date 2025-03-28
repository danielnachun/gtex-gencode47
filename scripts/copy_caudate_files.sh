#!/bin/bash
set -e

gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/coverage gs://fa-exchange2/caudate_analysis/coverage || exit 1
gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/flagstat gs://fa-exchange2/caudate_analysis/flagstat || exit 1
gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/gatk gs://fa-exchange2/caudate_analysis/gatk || exit 1
gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/leafcutter gs://fa-exchange2/caudate_analysis/leafcutter || exit 1
gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/rnaseq_qc gs://fa-exchange2/caudate_analysis/rnaseq_qc || exit 1
gsutil -m cp -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/rsem gs://fa-exchange2/caudate_analysis/rsem || exit 1
