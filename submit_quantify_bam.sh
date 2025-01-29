
# sample args
dir_prefix=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow
sample_id=GTEX-1A3MV-0005-SM-7PC1O
reference_fasta=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
chr_sizes=/oak/stanford/groups/smontgom/dnachun/data/gtex/data/references/GRCh38.chrsizes
genes_gtf=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/references/gencode.v47.annotation.gtf
# intervals_bed = TODO download correct interval_bed
vcf_dir=/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/processed/vcfs

export _JAVA_OPTIONS="-Xmx64g"
bash quantify_bam.sh \
    --dir_prefix ${dir_prefix} \
    --sample_id ${sample_id} \
    --reference_fasta ${reference_fasta} \
    --chr_sizes ${chr_sizes} \
    --genes_gtf ${genes_gtf} \
    --intervals_bed ${intervals_bed} \
    --vcf_dir ${vcf_dir}