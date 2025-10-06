#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

options_array=(
    tissue_id_list
    coloc_base_dir
    code_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'format_colocboosts' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --tissue_id_list )
            tissue_id_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --coloc_base_dir )
            coloc_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        -- ) 
            shift; break ;;
        * ) 
            echo "Invalid argument ${1} ${2}" >&2; 
            exit 1 ;;
    esac
done

# Activate environment
source <(pixi shell-hook --environment pecotmr-dev --manifest-path  ${code_dir}/pixi.toml)

# map job id to line number and then to tissue id
line_number="${SLURM_ARRAY_TASK_ID}"
tissue_id="$(sed "${line_number}q; d" "${tissue_id_list}")"

coloc_dir="${coloc_base_dir}/${tissue_id}"
echo "Finding files to process in ${coloc_dir}..."
tmpfile=$(mktemp)
count=0

# Count the number of .colocboost.rds and .txt files
num_rds=$(find "${coloc_dir}" -type f -name '*colocboost.rds' | wc -l)
num_txt=$(find "${coloc_dir}" -type f -name '*.txt' | wc -l)
echo "Number of .colocboost.rds files: $num_rds"
echo "Number of .txt files: $num_txt"

if [[ "$num_rds" -ne "$num_txt" ]]; then
    # Find all .rds files, but only keep those that do not already have a corresponding .txt output.
    find "${coloc_dir}" -type f -name '*.txt' | sed 's/\.txt$//' | sort -u > "${tmpfile}.txtprefixes"
    echo "Found $(wc -l < "${tmpfile}.txtprefixes") .txt prefixes"
    find "${coloc_dir}" -type f -name '*colocboost.rds' | while read -r rdsfile; do
        prefix="${rdsfile%.rds}"
        if ! grep -Fxq "$prefix" "${tmpfile}.txtprefixes"; then
            echo "$rdsfile"
        fi
    done > "$tmpfile"

    # Count how many files to process
    total_to_process=$(wc -l < "$tmpfile")
    echo "Found $total_to_process files to process, with rds but no txt"

    # Set the number of parallel jobs to one less than SLURM_CPUS_PER_TASK, with a minimum of 1
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
        num_parallel=$((SLURM_CPUS_PER_TASK - 1))
    else
        num_parallel=1
    fi

    echo "Processing files in parallel..."
    parallel --jobs "${num_parallel}" "${code_dir}/run_format_colocboost.sh" --file {} --code_dir "${code_dir}" :::: "$tmpfile"

    # Count the number of successfully processed files
    count=$(find "${coloc_dir}" -type f -name '*.txt' | wc -l)
    echo "Done. Processed $count files."

    # Cleanup
    rm -f "$tmpfile" "${tmpfile}.txtprefixes"
else
    echo "Number of .colocboost.rds files ($num_rds) matches number of .txt files ($num_txt). No processing needed."
fi

# Aggregate outputs via wrapper

"${code_dir}/run_combine_colocboost.sh" --tissue_id "${tissue_id}" --coloc_output_dir "${coloc_dir}" --code_dir "${code_dir}" --aggregated_output_dir "${coloc_base_dir}"
echo "Aggregation complete."
