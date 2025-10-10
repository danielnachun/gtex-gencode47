#!/usr/bin/env bash
set  -o nounset -o pipefail -o errexit


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
    tissue_id_list
    coloc_base_dir
    code_dir
    regenerate
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'aggregate_colocboosts' -- "$@")
eval set -- "${arguments}"

regenerate="false"

while true; do
    case "${1}" in
        --tissue_id_list )
            tissue_id_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --coloc_base_dir )
            coloc_base_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --regenerate )
            regenerate="${2}"; shift 2 ;;
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

# Create format_file_list directory if it doesn't exist
mkdir -p "${coloc_base_dir}/format_file_list"

tmpfile="${coloc_base_dir}/format_file_list/${tissue_id}.${SLURM_ARRAY_TASK_ID}.file_list.tmp"
count=0

if [[ "${regenerate}" == "true" ]]; then
    # Regenerate: process all .colocboost.rds files, regardless of .txt files
    find "${coloc_dir}" -type f -name '*colocboost.rds' > "$tmpfile"
    total_to_process=$(wc -l < "$tmpfile")
    echo "Regenerate mode: Found $total_to_process .colocboost.rds files to process."

    # Set the number of parallel jobs to one less than SLURM_CPUS_PER_TASK, with a minimum of 1
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
        num_parallel=$((SLURM_CPUS_PER_TASK - 1))
    else
        num_parallel=1
    fi

    echo "Processing files in parallel (regenerate mode)..."
    parallel --jobs "${num_parallel}" "${code_dir}/run_format_colocboost.sh" --file {} --code_dir "${code_dir}" :::: "$tmpfile"

    count=$(find "${coloc_dir}" -type f -name '*.txt' | wc -l)
    echo "Done. Processed $count files."

    rm -f "$tmpfile"
else
    # Only process .colocboost.rds files that are missing the robust output
    # i.e., rerun if <prefix>.xqtl_coloc.robust.txt does NOT exist
    num_rds=$(find "${coloc_dir}" -type f -name '*colocboost.rds' | wc -l)
    echo "Number of .colocboost.rds files: $num_rds"

    # Fast path: if we already have exactly two .txt per .rds (non-robust + robust), skip formatting
    num_txt=$(find "${coloc_dir}" -type f -name '*colocboost*.txt' | wc -l)
    if [[ "$num_txt" -eq $(( 2 * num_rds )) ]]; then
        echo "Detected $num_txt .txt files for $num_rds .rds files (2x). Skipping format step."
        rm -f "$tmpfile"
        # proceed directly to aggregation step after this block
    else

    # Build list of RDS files lacking robust txt
    find "${coloc_dir}" -type f -name '*colocboost.rds' | while read -r rdsfile; do
        prefix="${rdsfile%.rds}"
        robust_txt="${prefix}.xqtl_coloc.robust.txt"
        if [[ ! -f "$robust_txt" ]]; then
            # robust output missing; schedule for processing (even if non-robust exists)
            echo "$rdsfile"
        fi
    done > "$tmpfile"

    total_to_process=$(wc -l < "$tmpfile")
    echo "Found $total_to_process files to process missing robust outputs"

    if [[ "$total_to_process" -gt 0 ]]; then
        if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
            num_parallel=$((SLURM_CPUS_PER_TASK - 1))
        else
            num_parallel=1
        fi

        echo "Processing files in parallel..."
        parallel --jobs "${num_parallel}" "${code_dir}/run_format_colocboost.sh" --file {} --code_dir "${code_dir}" :::: "$tmpfile"

        # Count only robust outputs for reporting
        count=$(find "${coloc_dir}" -type f -name '*.xqtl_coloc.robust.txt' | wc -l)
        echo "Done. Robust outputs present: $count files."
    else
        echo "All .colocboost.rds files have robust outputs. No processing needed."
    fi

    rm -f "$tmpfile"
    fi
fi

# Aggregate outputs via wrapper

"${code_dir}/run_combine_colocboost.sh" --tissue_id "${tissue_id}" --coloc_output_dir "${coloc_dir}" --code_dir "${code_dir}" --aggregated_output_dir "${coloc_base_dir}"
echo "Aggregation complete."
