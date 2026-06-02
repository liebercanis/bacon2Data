#!/bin/bash
# Note: To add date/time to job name, submit with: sbatch --job-name=postAna_mgold_$(date +%Y%m%d_%H%M%S) submit_post.sh <tag>
## #SBATCH --array=0-31%8
#SBATCH --job-name=postAna_mgold
#SBATCH -N 1
#SBATCH -q shared
#SBATCH -C cpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --output=logs/postAna_%j.log
#SBATCH --error=logs/postAna_%j.err
#SBATCH --mail-user=mgold@unm.edu
#SBATCH --mail-type=ALL
#SBATCH -A m2676

if [ "$#" -lt 1 ]; then
    echo " usage: sbatch submit_post.sh <tag>"
    echo " To add date/time to job name, submit with: sbatch --job-name=postAna_mgold_$(date +%Y%m%d_%H%M%S) submit_post.sh <tag>"
    exit 1
fi


# Load necessary modules (e.g., GCC, OpenMPI)
# module load cpe/23.12
# module load <other-needed-modules>
module load PrgEnv-gnu

WORK_DIR="/global/homes/m/mgold/mgold/bacon2Data/compiled"

cd "$WORK_DIR"

# Run your executable
echo "Job submitted: $(date '+%Y-%m-%d %H:%M:%S') | Starting at $(date) | tag: $1 | id: ${SLURM_JOB_ID}"
srun -n 1 -c 1 postAna $1


