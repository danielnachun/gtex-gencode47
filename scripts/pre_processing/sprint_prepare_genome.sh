# see https://github.com/jumphone/SPRINT/tree/master
source <(pixi shell-hook --environment calledsites --manifest-path /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing/pixi.toml)

sprint prepare \
    -t /home/klawren/oak/gtex/data/realign_references/gencode.v47.annotation.gtf \
    /home/klawren/oak/gtex/data/edsite_references/Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta \
    /oak/stanford/groups/smontgom/dnachun/data/gtex/v10/scripts/processing/.pixi/envs/calledsites/bin/bwa
