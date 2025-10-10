#!/usr/bin/env bash

set -o nounset -o pipefail -o errexit

# Wrapper script for combining formatted coloc outputs

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

options_array=(
    tissue_id
    coloc_output_dir
    code_dir
    aggregated_output_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'run_combine_colocboost' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --tissue_id )
            tissue_id="${2}"; shift 2 ;;
        --coloc_output_dir )
            coloc_output_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --aggregated_output_dir )
            aggregated_output_dir="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

aggregated_output_dir="${aggregated_output_dir:-}" 
tissue_id="${tissue_id:-}"

echo "Python version: $(python --version 2>&1)"
echo "Python path: $(which python)"
python "${code_dir}/combine_colocboost.py" \
  --coloc_output_dir "${coloc_output_dir}" \
  ${aggregated_output_dir:+--aggregated_output_dir "${aggregated_output_dir}"} \
  ${tissue_id:+--tissue_id "${tissue_id}"} \
  --quiet
