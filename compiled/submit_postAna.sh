#!/bin/bash
#SBATCH --qos=debug
#SBATCH --constraint=cpu
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -A m2676


# Usage: sbatch submit_postAna.sh <theStartTag> <theEndTag> <maxEntry>
# Example: sbatch submit_postAna.sh 09_10_2024 09_15_2024 100000

set -e

# Parse command-line arguments
THE_START_TAG="${1:?Error: theStartTag not provided. Usage: sbatch submit_postAna.sh <theStartTag> <theEndTag> <maxEntry>}"
THE_END_TAG="${2:?Error: theEndTag not provided. Usage: sbatch submit_postAna.sh <theStartTag> <theEndTag> <maxEntry>}"
MAX_ENTRY="${3:?Error: maxEntry not provided. Usage: sbatch submit_postAna.sh <theStartTag> <theEndTag> <maxEntry>}"

WORK_DIR="/global/homes/m/mgold/mgold/bacon2Data/compiled"
cd "$WORK_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=========================================="
echo "postAna SLURM Job Started"
echo "=========================================="
echo "Start Date Tag:  $THE_START_TAG"
echo "End Date Tag:    $THE_END_TAG"
echo "Max Entries:     $MAX_ENTRY"
echo "Working Dir:     $WORK_DIR"
echo "Job ID:          $SLURM_JOB_ID"
echo "Started at:      $(date)"
echo "=========================================="

# Run postAna executable with arguments
$BACON2DIR/compiled/postAna "$THE_START_TAG" "$THE_END_TAG" "$MAX_ENTRY"

echo "=========================================="
echo "postAna completed successfully at $(date)"
echo "=========================================="
