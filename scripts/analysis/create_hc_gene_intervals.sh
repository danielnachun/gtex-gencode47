#!/usr/bin/env bash
set -o xtrace -o nounset -o errexit

gencode_bed=/home/klawren/oak/gtex/data/realign_references/gencode.v47.annotation.bed
# sort the bed
sort -k1,1 -k2,2n ${gencode_bed}> /home/klawren/oak/gtex/data/other_references/gencode.v47.annotation.sorted.bed
# merge the bed
bedtools merge -i /home/klawren/oak/gtex/data/other_references/gencode.v47.annotation.sorted.bed > /home/klawren/oak/gtex/data/edsite_references/gencode.v47.merged.bed