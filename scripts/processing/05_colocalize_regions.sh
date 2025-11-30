#!/usr/bin/env bash

# Unified co-localization script for all analysis modes
# Usage: 05_colocalize_regions.sh [parameters...]

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o pipefail -o errexit


# Helper functions
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

# Define all possible parameters
options_array=(
    ld_region_list
    ld_region
    tissue_id_list
    genotype_stem
    covariate_dir
    expression_dir
    all_v39_genes_path
    region_padding
    output_dir
    code_dir
    # GWAS-specific parameters (optional)
    gwas_id_list
    gwas_dir
    gwas_meta
    ld_meta
    gwas_column_matching
    strong_only
    p_threshold
    # Analysis mode flags
    multi_tissue
    run_v39
    run_individual
    run_xqtl_only
    run_separate_gwas
    run_joint_gwas
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'coloc' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --ld_region_list ) ld_region_list="${2}"; if [ "${2}" != "none" ]; then check_for_file "${1}" "${2}"; fi; shift 2 ;;
        --ld_region ) ld_region="${2}"; shift 2 ;;
        --tissue_id_list ) tissue_id_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --genotype_stem ) genotype_stem="${2}"; shift 2 ;;
        --covariate_dir ) covariate_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --expression_dir ) expression_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --all_v39_genes_path ) all_v39_genes_path="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --region_padding ) region_padding="${2}"; shift 2 ;;
        --output_dir ) output_dir="${2}"; shift 2 ;;
        --code_dir ) code_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        # GWAS-specific parameters
        --gwas_id_list ) gwas_id_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_dir ) gwas_dir="${2}"; check_for_directory "${1}" "${2}"; shift 2 ;;
        --gwas_meta ) gwas_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --ld_meta ) ld_meta="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --gwas_column_matching ) gwas_column_matching="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --strong_only ) strong_only="${2}"; shift 2 ;;
        --p_threshold ) p_threshold="${2}"; shift 2 ;;
        # Analysis mode flags
        --multi_tissue ) multi_tissue="${2}"; shift 2 ;;
        --run_v39 ) run_v39="${2}"; shift 2 ;;
        --run_individual ) run_individual="${2}"; shift 2 ;;
        --run_xqtl_only ) run_xqtl_only="${2}"; shift 2 ;;
        --run_separate_gwas ) run_separate_gwas="${2}"; shift 2 ;;
        --run_joint_gwas ) run_joint_gwas="${2}"; shift 2 ;;
        --) shift; break;;
        * ) echo "Invalid argument ${1} ${2}" >&2; exit 1
    esac
done

# Set default values for analysis flags
multi_tissue=${multi_tissue:-FALSE}
run_v39=${run_v39:-FALSE}
run_individual=${run_individual:-FALSE}
run_xqtl_only=${run_xqtl_only:-FALSE}
run_separate_gwas=${run_separate_gwas:-FALSE}
run_joint_gwas=${run_joint_gwas:-FALSE}

# Set default values for GWAS filtering parameters
strong_only=${strong_only:-FALSE}
p_threshold=${p_threshold:-5e-8}

# Activate pixi environment to get tabix and other tools
source <(pixi shell-hook --environment pecotmr-dev --manifest-path "${code_dir}/pixi.toml")


# Get LD region - either from parameter or from ld_region_list using SLURM_ARRAY_TASK_ID
if [ -n "${ld_region:-}" ]; then
    # Use provided ld_region parameter
    echo "Using provided ld_region: ${ld_region}"
    # Create a unique identifier for working directory from the region string
    # Use a hash to ensure uniqueness even with parallel execution
    line_number=$(echo "${ld_region}" | md5sum | cut -d' ' -f1 | head -c 16)
else
    # Get from ld_region_list using SLURM array task ID
    if [ -z "${ld_region_list:-}" ] || [ "${ld_region_list}" = "none" ]; then
        echo "Error: Must provide either --ld_region or --ld_region_list"
        exit 1
    fi
    line_number="${SLURM_ARRAY_TASK_ID}"
    ld_region="$(sed "${line_number}q; d" "${ld_region_list}")"
fi

# Create working directory with unique name to avoid conflicts in parallel execution
working_dir="${TMPDIR}/coloc_${line_number}_$$"
mkdir -p "${working_dir}"

# GWAS filtering logic (if running GWAS analysis)
if [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_joint_gwas}" = "TRUE" ] || [ "${run_joint_gwas}" = "true" ]; then
    # Function to check for strong associations in GWAS file
    check_strong_associations() {
        local gwas_file="$1"
        local region_chr="$2"
        local region_start="$3"
        local region_end="$4"
        
        # Convert both values to numbers for proper comparison with scientific notation
        local tabix_output=$(tabix -h "${gwas_file}" "${region_chr}:${region_start}-${region_end}" 2>&1)
        local tabix_exit=$?
        
        if [ "${tabix_exit}" -ne 0 ] && [ -n "${tabix_output}" ]; then
            # Tabix error - return failure
            return 1
        fi
        
        local strong_count=$(echo "${tabix_output}" | \
           awk -v p_thresh="${p_threshold}" 'BEGIN {p_thresh_num = p_thresh + 0} ($11+0) < p_thresh_num {count++} END {print count+0}')
        
        # Check if any variants have strong associations
        if [ "${strong_count}" -gt 0 ]; then
            return 0  # Has strong associations
        else
            return 1  # No strong associations
        fi
    }
    
    # Create filtered GWAS ID list if strong_only is enabled
    if [ "${strong_only}" = "TRUE" ] || [ "${strong_only}" = "true" ]; then
        echo "Filtering GWAS files for strong associations (p < ${p_threshold}) in region ${ld_region}"

        # Parse region coordinates for filtering
        region_chr=$(echo "${ld_region}" | cut -d: -f1)
        region_start=$(echo "${ld_region}" | cut -d: -f2 | cut -d- -f1)
        region_end=$(echo "${ld_region}" | cut -d: -f2 | cut -d- -f2)
        region_chr_tabix=$(echo "${region_chr}" | sed 's/^chr//')
        
        filtered_gwas_id_list="${working_dir}/filtered_gwas_id_list.txt"
        rm -f "${filtered_gwas_id_list}"
        touch "${filtered_gwas_id_list}"
        
        original_count=$(wc -l < "${gwas_id_list}")
        filtered_count=0
        
        while read -r gwas_id; do
            gwas_file="${gwas_dir}/imputed_${gwas_id}.txt.gz"
            if [ -f "${gwas_file}" ]; then
                if check_strong_associations "${gwas_file}" "${region_chr_tabix}" "${region_start}" "${region_end}"; then
                    echo "${gwas_id}" >> "${filtered_gwas_id_list}"
                    echo "  ${gwas_id}: PASSED (has strong associations)"
                    filtered_count=$((filtered_count + 1))
                else
                    echo "  ${gwas_id}: FILTERED (no strong associations)"
                fi
            else
                echo "  ${gwas_id}: SKIPPED (file not found: ${gwas_file})"
            fi
        done < "${gwas_id_list}"
        
        # Update gwas_id_list to use filtered list
        gwas_id_list="${filtered_gwas_id_list}"
        echo "Filtered from ${original_count} to ${filtered_count} GWAS files"
        
        # Exit if no GWAS files passed the filter
        if [ "${filtered_count}" -eq 0 ]; then
            echo "No GWAS files with strong associations found in region ${ld_region}"
            echo "Skipping colocalization analysis for this region"
            exit 0
        fi
    fi
fi

echo $(date +"[%b %d %H:%M:%S] Running colocboost with and without v39 genes on region ${ld_region}")

# Read tissue IDs (preserving spaces in tissue names)
mapfile -t tissue_ids < "${tissue_id_list}"

# Unified colocboost processing function
process_colocboost() {
    if [ "${multi_tissue}" = "TRUE" ] || [ "${multi_tissue}" = "true" ]; then
        echo "Processing multi tissue for region ${ld_region}"
        
        # Use original paths directly (files read once, network I/O acceptable, avoids rsync overhead)
        # Prepare aggregated files using original paths
        gene_bed_list_all="${working_dir}/all_tissues.gene_bed_list.txt"
        covariate_list_all="${working_dir}/all_tissues.covariate_list.txt"
        
        printf "%s\n" "${tissue_ids[@]}" | sed "s|^|${expression_dir}/|;s|\$|.v11.normalized_expression.bed.gz|" > "${gene_bed_list_all}"
        printf "%s\n" "${tissue_ids[@]}" | sed "s|^|${covariate_dir}/|;s|\$|.v11.covariates.txt|" > "${covariate_list_all}"
        conditions_list_individual_all="${working_dir}/all_tissues.conditions_list_individual.txt"
        printf "%s\n" "${tissue_ids[@]}" > "${conditions_list_individual_all}"
        
        echo "Running colocboost aggregated on LD region ${ld_region} across all tissues"
        
        # Build base command
        local cmd="${code_dir}/run_05_colocboost_compare.sh \
            --tissue_id \"all_tissues\" \
            --ld_region \"${ld_region}\" \
            --association_region \"${ld_region}\" \
            --gene_region \"${ld_region}\" \
            --gene_bed_list \"${gene_bed_list_all}\" \
            --covariate_list \"${covariate_list_all}\" \
            --covariate_path none \
            --conditions_list_individual \"${conditions_list_individual_all}\" \
            --genotype_stem \"${genotype_stem}\" \
            --v39_gene_id_path \"${all_v39_genes_path}\" \
            --run_v39 \"${run_v39}\" \
            --run_individual \"${run_individual}\" \
            --run_xqtl_only \"${run_xqtl_only}\" \
            --run_separate_gwas \"${run_separate_gwas}\" \
            --run_joint_gwas \"${run_joint_gwas}\" \
            --output_dir \"${output_dir}/all_tissues\" \
            --code_dir \"${code_dir}\""
        
        # Add GWAS-specific parameters if running GWAS analysis
        if [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_joint_gwas}" = "TRUE" ] || [ "${run_joint_gwas}" = "true" ]; then
            cmd="${cmd} \
                --gwas_id_list \"${gwas_id_list}\" \
                --gwas_phenotype_list \"${gwas_id_list}\" \
                --gwas_meta \"${gwas_meta}\" \
                --ld_meta \"${ld_meta}\" \
                --gwas_column_matching \"${gwas_column_matching}\""
        fi
        
        # Execute the command
        eval "${cmd}"
    else
        echo "Processing single tissue for region ${ld_region}"
        
        # Set up parallel processing
        if [ -n "${OMP_NUM_THREADS:-}" ] && [ "${OMP_NUM_THREADS}" -gt 1 ]; then
            parallel_jobs=$((OMP_NUM_THREADS - 1))
        else
            parallel_jobs=1
        fi
        
        # Define unified process_tissue function
        process_tissue() {
            local tissue_id="${1}"
            
            echo "Running colocboost on LD region ${ld_region} in tissue ${tissue_id}"
            
            # Build base command
            # Create conditions_list_individual file for single tissue
            conditions_list_individual_file="${working_dir}/${tissue_id}.conditions_list_individual.txt"
            echo "${tissue_id}" > "${conditions_list_individual_file}"
            
            local cmd="${code_dir}/run_05_colocboost_compare.sh \
                --tissue_id \"${tissue_id}\" \
                --ld_region \"${ld_region}\" \
                --association_region \"${ld_region}\" \
                --gene_region \"${ld_region}\" \
                --gene_bed_list \"${expression_dir}/${tissue_id}.v11.normalized_expression.bed.gz\" \
                --covariate_list \"${covariate_dir}/${tissue_id}.v11.covariates.txt\" \
                --covariate_path none \
                --conditions_list_individual \"${conditions_list_individual_file}\" \
                --genotype_stem \"${genotype_stem}\" \
                --v39_gene_id_path \"${all_v39_genes_path}\" \
                --run_v39 \"${run_v39}\" \
                --run_individual \"${run_individual}\" \
                --run_xqtl_only \"${run_xqtl_only}\" \
                --run_separate_gwas \"${run_separate_gwas}\" \
                --run_joint_gwas \"${run_joint_gwas}\" \
                --output_dir \"${output_dir}/${tissue_id}\" \
                --code_dir \"${code_dir}\""
            
            # Add GWAS-specific parameters if running GWAS analysis
            if [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_joint_gwas}" = "TRUE" ] || [ "${run_joint_gwas}" = "true" ]; then
                cmd="${cmd} \
                    --gwas_id_list \"${gwas_id_list}\" \
                    --gwas_phenotype_list \"${gwas_id_list}\" \
                    --gwas_meta \"${gwas_meta}\" \
                    --ld_meta \"${ld_meta}\" \
                    --gwas_column_matching \"${gwas_column_matching}\""
            fi
            
            # Execute the command
            eval "${cmd}"
        }
        
        export -f process_tissue
        
        echo "Running process_tissue in parallel for ${parallel_jobs} jobs, across ${#tissue_ids[@]} tissues"
        
        # Run process_tissue in parallel for each tissue_id
        printf "%s\n" "${tissue_ids[@]}" | parallel --tag --line-buffer -j "${parallel_jobs}" process_tissue {} 2>&1
    fi
}

# Validate analysis mode and run processing
if [ "${run_separate_gwas}" = "TRUE" ] || [ "${run_separate_gwas}" = "true" ] || [ "${run_xqtl_only}" = "TRUE" ] || [ "${run_xqtl_only}" = "true" ]; then
    process_colocboost
else
    echo "Error: Must specify either run_separate_gwas=TRUE or run_xqtl_only=TRUE"
    exit 1
fi

# Create completion marker using normalized region string for filesystem safety
# Normalize the region string to handle colons and other special characters
completion_dir="${output_dir}/completed"
mkdir -p "${completion_dir}"
normalized_region=$(echo "$ld_region" | tr '\t' '_' | tr '\n' ' ' | sed 's/[:/]/_/g' | tr -d '\r' | sed 's/  */_/g' | sed 's/^_\|_$//g')
completion_file="${completion_dir}/${normalized_region}.completed"
echo "Processing completed successfully for ld region ${ld_region}" > "${completion_file}"
echo "Completion marker created: ${completion_file}"
