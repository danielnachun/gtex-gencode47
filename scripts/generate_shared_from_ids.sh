#!/usr/bin/env bash

sample_ids="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/caudate/caudate_samples.txt"
participant_ids="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/caudate/caudate_participants_unique.txt"
output_file="/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/other_references/caudate/caudate_shared_samples.txt"

# clear output file
> "$output_file"

awk '{
    if (NR==FNR) {
        participants[$1]=1
    } else {
        split($1,a,"-")
        if ((a[1]"-"a[2]) in participants) {
            print $0
        }
    }
}' "$participant_ids" "$sample_ids" > "$output_file"