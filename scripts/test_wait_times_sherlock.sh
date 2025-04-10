#!/bin/bash
#SBATCH --job-name=cpu_wait_time_test
#SBATCH --output=cpu_wait_time_test_%j.out
#SBATCH --error=cpu_wait_time_test_%j.err
#SBATCH --time=00:10:00  # Set a time limit for the job
#SBATCH --partition=normal,owners  # Replace with your partition name

# Function to submit a job and measure wait time
submit_job() {
    local cpus=$1
    local job_name="test_job_${cpus}_cpus"
    
    # Record the start time
    start_time=$(date +%s)

    # Submit the job
    job_id=$(sbatch --job-name="$job_name" --cpus-per-task="$cpus" --wrap="sleep 60" | awk '{print $4}')

    # Wait for the job to complete
    squeue -j "$job_id" -h | while read -r line; do
        sleep 1
        if [ -z "$(squeue -j "$job_id" -h)" ]; then
            break
        fi
    done

    # Record the end time
    end_time=$(date +%s)

    # Calculate wait time
    wait_time=$((end_time - start_time))
    echo "Wait time for job requesting $cpus CPUs: $wait_time seconds" >> wait_times.log
}

# Clear previous log
echo "Wait Times Log" > wait_times.log

# Submit jobs for 1, 8, and 16 CPUs
submit_job 1
submit_job 8
submit_job 16