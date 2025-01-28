
dir_prefix = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/test_workflow
sample_id = GTEX-1A3MV-0005-SM-7PC1O
reference_fasta = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta
rsem_ref_dir = /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/data/references/rsem_reference_GRCh38_gencode47
star_index = /oak/stanford/groups/smontgom/dnachun/data/gtex/data/references/STAR_genome_GRCh38_noALT_noHLA_noDecoy_v47_oh75

./realign_bam.sh \
    ${dir_prefix} \
    ${sample_id} \
    ${reference_fasta} \
    ${rsem_ref_dir} \
    ${star_index}