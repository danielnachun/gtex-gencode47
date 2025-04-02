#!/bin/bash

gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/coverage gs://fa-exchange2/caudate/coverage 
gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/flagstat gs://fa-exchange2/caudate/flagstat 
gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/gatk gs://fa-exchange2/caudate/gatk 
gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/leafcutter gs://fa-exchange2/caudate/leafcutter 
gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/rnaseq_qc gs://fa-exchange2/caudate/rnaseq_qc
gsutil -m cp -c -r /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/caudate_analysis/rsem gs://fa-exchange2/caudate/rsem 
