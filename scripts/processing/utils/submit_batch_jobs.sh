#!/usr/bin/env bash

set -o nounset -o errexit

# Parse command line arguments with getopt
options_array=(
    config_file
    batch_script
    job_name
    items_list_file
    processing_script
    file_param
    completion_dir
    regenerate_all
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'submit_batch_jobs' -- "$@")
eval set -- "${arguments}"

file_param=${file_param-}
regenerate_all=${regenerate_all-false}

while true; do
    case "${1}" in
        --config_file )
            CONFIG_FILE="${2}"; shift 2 ;;
        --batch_script )
            BATCH_SCRIPT="${2}"; shift 2 ;;
        --job_name )
            JOB_NAME="${2}"; shift 2 ;;
        --items_list_file )
            ITEMS_LIST_FILE="${2}"; shift 2 ;;
        --processing_script )
            PROCESSING_SCRIPT="${2}"; shift 2 ;;
        --file_param )
            file_param="${2}"; shift 2 ;;
        --completion_dir )
            completion_dir="${2}"; shift 2 ;;
        --regenerate_all )
            regenerate_all="${2}"; shift 2 ;;
        --)
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

ADDITIONAL_PARAMS=("$@")

# Source the config file
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE" || { echo "Error: Config file $CONFIG_FILE not found!"; exit 1; }

# Set defaults for variables that might not be in config
num_parallel=${num_parallel:-1}
max_array_size=${max_array_size:-1000}

# Get items to process, filtering out completed items unless regenerate_all is true
all_items=$(cat "${ITEMS_LIST_FILE}")

if [ "${regenerate_all}" = "true" ] || [ "${regenerate_all}" = "TRUE" ]; then
    items_to_process="$all_items"
else
    check_dir="${completion_dir:-${output_dir}/completed/${JOB_NAME}}"
    completed_items=$(find "${check_dir}" -name "*.completed" -exec basename {} \; 2>/dev/null | sed 's|\.completed$||')
    
    if [ -n "$completed_items" ]; then
        items_to_process=$(echo "$all_items" | grep -v -F -f <(echo "$completed_items"))
    else
        items_to_process="$all_items"
    fi
fi

if [ -z "$items_to_process" ]; then
    echo "To be processed: 0"
    echo "All items processed"
    exit 0
fi

to_process_count=$(echo "$items_to_process" | wc -l)
original_count=$(grep -v '^$' "${ITEMS_LIST_FILE}" | wc -l)
completed_count=$((original_count - to_process_count))

num_series=$(( (to_process_count + max_array_size - 1) / max_array_size ))
num_series=$(( (num_series + num_parallel - 1) / num_parallel ))
num_batches=$(( (to_process_count + num_parallel * num_series - 1) / (num_parallel * num_series) ))
final_items_per_batch=$(( (to_process_count + num_batches - 1) / num_batches ))

echo "Completed: ${completed_count}/${original_count}, To process: ${to_process_count}/${original_count}"
echo "Array jobs: ${num_batches} array jobs with ${final_items_per_batch} items per batch, ${num_parallel} parallel jobs * ${num_series} series jobs"

# Create batch files
file_list_folder="$output_dir/file_lists/file_lists_${JOB_NAME}"
rm -rf "${file_list_folder}"
mkdir -p "${file_list_folder}"

echo "$items_to_process" | split -l "${final_items_per_batch}" --additional-suffix=".txt" - "${file_list_folder}/file_list_"

file_list_paths="${output_dir}/file_lists/file_list_paths_${JOB_NAME}.txt"
rm -rf "${file_list_paths}"
printf "%s\n" "${file_list_folder}"/* > "${file_list_paths}"

echo "Batches created: ${num_batches}"

# Configure SLURM job parameters
sbatch_params=(
    --output "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --error "${output_dir}/logs/${JOB_NAME}/%A/%A_%a.log"
    --array "1-${num_batches}%250"
    --cpus-per-task "${num_parallel}"
    --time "${job_time:-12:00:00}"
    --mem "${job_mem:-256G}"
    --job-name "${JOB_NAME}_batch"
    "${BATCH_SCRIPT}" \
    --file_list_paths "${file_list_paths}" \
    --num_parallel "${num_parallel}" \
    --processing_script "${PROCESSING_SCRIPT}" \
    --file_param "${file_param:---bam_file}" \
    -- \
    "${ADDITIONAL_PARAMS[@]}"
)

if [ -n "${job_partition:-}" ]; then
    sbatch_params+=(--partition "${job_partition}")
fi
if [ -n "${job_account:-}" ]; then
    sbatch_params+=(--account "${job_account}")
fi
if [ -n "${job_constraint:-}" ]; then
    sbatch_params+=(--constraint "${job_constraint}")
fi
if [ -n "${job_tmp:-}" ]; then
    sbatch_params+=(--tmp "${job_tmp}")
fi

# Check if job requires long QOS
is_long_job=false
job_time_value="${job_time:-12:00:00}"
if [[ "$job_time_value" =~ ^[0-9]+- ]]; then
    days=$(echo "$job_time_value" | cut -d'-' -f1)
    if [ "$days" -gt 2 ]; then
        is_long_job=true
    fi
elif [[ "$job_time_value" =~ ^[0-9]+:[0-9]+:[0-9]+$ ]]; then
    hours=$(echo "$job_time_value" | cut -d':' -f1)
    if [ "$hours" -gt 48 ]; then
        is_long_job=true
    fi
fi

# Submit job array
if [ "${submit_on}" = 'sherlock' ]; then
    sherlock_params=(--tmp 200G)
    if [ "$is_long_job" = true ]; then
        sherlock_params+=(--partition sfgf,biochem)
    else
        sherlock_params+=(--partition normal,owners)
    fi
    sherlock_params+=(--exclude pritch)
    sbatch "${sherlock_params[@]}" "${sbatch_params[@]}"
elif [ "${submit_on}" = 'scg' ]; then
    sbatch --account smontgom --partition batch --constraint "nvme" "${sbatch_params[@]}"
else
    echo "must submit on either 'sherlock' or 'scg'"
fi
