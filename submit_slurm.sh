#!/bin/bash
#SBATCH --job-name=postAna_mgold
#SBATCH -N 1
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --output=logs/anaDir_%j.log
#SBATCH --error=logs/anacDir_%j.err
#SBATCH --array=0-31%8
#SBATCH --mail-user=mgold@unm.edu
#SBATCH --mail-type=ALL
#SBATCH -A m2676


# Usage: sbatch --array=0-<numfiles-1>%<parallel_jobs> submit_slurm.sh <date_tag>
# Example: sbatch --array=0-9%4 submit_slurm.sh 09_10_2024

set -e

DATE_TAG="${1:?Error: date tag not provided. Usage: sbatch submit_slurm.sh <date_tag>}"
WORK_DIR="/global/homes/m/mgold/mgold/bacon2Data"

cd "$WORK_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

# Get list of files matching the date tag
files=(rootData/*"$DATE_TAG"*)

# Check if any files were found
if [[ ${#files[@]} -eq 0 ]]; then
    echo "Error: No files found matching tag '$DATE_TAG'" >&2
    exit 1
fi

# Get the file for this array task
file="${files[$SLURM_ARRAY_TASK_ID]}"

# Check if file exists
if [[ ! -f "$file" ]]; then
    echo "Error: File not found: $file" >&2
    exit 1
fi

echo "Processing file: $file"
echo "Task ID: $SLURM_ARRAY_TASK_ID / ${#files[@]}"
echo "Starting at $(date)"

# Set library path and run analysis
export LD_LIBRARY_PATH="$WORK_DIR:$LD_LIBRARY_PATH"
"$WORK_DIR/compiled/anacg" "$file"

echo "Completed at $(date)"
