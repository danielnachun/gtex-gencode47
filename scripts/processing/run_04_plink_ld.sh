#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

options_array=(
    genotype_prefix
    sample_ids
    output_dir
    chr_id
    from_bp
    to_bp
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'fastq_to_star' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --genotype_prefix )
            genotype_prefix="${2}"; shift 2 ;;
        --sample_ids )
            sample_ids="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_dir )
            output_dir="${2}"; shift 2 ;;
        --chr_id )
            chr_id="${2}"; shift 2 ;;
        --from_bp )
            from_bp="${2}"; shift 2 ;;
        --to_bp )
            to_bp="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

THREADS="${THREADS:-8}"

mkdir -p ${output_dir}
echo $(date +"[%b %d %H:%M:%S] Running PLINK LD for chr${chr_id}:${from_bp}-${to_bp}")

# Create the variant list with quality control filters
plink \
  --bfile "$genotype_prefix" \
  --keep "$sample_ids" \
  --chr "$chr_id" \
  --from-bp "$from_bp" \
  --to-bp "$to_bp" \
  --maf 0.009 \
  --geno 0.05 \
  --hwe 1e-6 \
  --make-bed \
  --threads "$THREADS" \
  --out "$output_dir/LD_chr${chr_id}_${from_bp}_${to_bp}"

# Generate LD matrix from the filtered variant list
plink \
  --bfile "$output_dir/LD_chr${chr_id}_${from_bp}_${to_bp}" \
  --r square gz \
  --threads "$THREADS" \
  --out "$output_dir/LD_chr${chr_id}_${from_bp}_${to_bp}"



# Convert gzip compression to xz compression for pecotmr compatibility
ld_file="$output_dir/LD_chr${chr_id}_${from_bp}_${to_bp}.ld.gz"
if [[ -f ${ld_file} ]]; then
    echo $(date +"[%b %d %H:%M:%S] Converting compression format from gzip to xz for R compatibility")
    
    # Decompress the gzip file and recompress with xz
    gunzip -c ${ld_file} | xz -c > ${ld_file}.tmp
    
    # Replace the original file with the xz-compressed version
    mv ${ld_file}.tmp ${ld_file}
    
    echo $(date +"[%b %d %H:%M:%S] Compression conversion completed")
else
    echo "Warning: LD file ${ld_file} not found for compression conversion"
fi