#!/usr/bin/env bash

# Set bash options to fail immediately on errors or if variables are undefined.
set -o nounset -o pipefail -o errexit

# Wrapper script for processing individual files with format_colocboost.R

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
    file
    code_dir
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'run_format_colocboost' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --file )
            file="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --code_dir )
            code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

out_prefix="${file%.rds}"

# Process the file
if Rscript "${code_dir}/format_colocboost.R" --pecotmr_colocboost "${file}" --output_prefix "${out_prefix}" --quiet; then
  exit 0
else
  echo "Error: Failed to process ${file}"
  exit 1
fi
