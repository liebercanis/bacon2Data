#!/usr/bin/env python3
"""
SLURM wrapper for anaDir - submits analysis jobs to SLURM cluster
Usage: python3 submit_slurm.py <date_tag> [max_files] [--parallel=N] [--time=HH:MM:SS]
"""

import sys
import os
import subprocess
import argparse
from pathlib import Path


def get_matching_files(date_tag, rootdata_dir='rootData'):
    """Get list of files matching the date tag"""
    files = []
    try:
        for filename in os.listdir(rootdata_dir):
            if date_tag in filename:
                files.append(filename)
    except FileNotFoundError:
        print(f"Error: Directory '{rootdata_dir}' not found", file=sys.stderr)
        return []
    
    return sorted(files)


def submit_slurm_job(date_tag, num_files=None, mem="32G", parallel_jobs=8, time_limit="01:00:00", account="m2676", queue="shared"):
    """Submit SLURM job array"""
   
    rootData = os.getenv("ROOTDATA")
    files = get_matching_files(date_tag,rootData)

    for(i, f) in enumerate(files):
        print(f" file {i}  file {f} ")
    
    if not files:
        print(f"Error: No files found matching tag '{date_tag}'", file=sys.stderr)
        return 1
    
    # Limit number of files if specified
    if num_files is not None:
        files = files[:num_files]
    
    num_tasks = len(files)
    print(f"Found {num_tasks} files matching '{date_tag}'")
    print(f"Submitting {num_tasks} tasks with max {parallel_jobs} parallel jobs")
    
    # Create logs directory
    os.makedirs('logs', exist_ok=True)
    
    # Build sbatch command
    # The special Python variable __file__ contains the pathname of the file from which the module was loaded. 
    script_path = os.path.join(os.path.dirname(__file__), 'submit_slurm.sh')

    print('script_path= ',script_path) 

    sbatch_cmd = [
        'sbatch',
        f'--job-name=postAna_{date_tag}',
        f'--array=0-{num_tasks-1}%{parallel_jobs}',
        f'--mem={mem}',
        f'--time={time_limit}',
        f'-A {account}',
        f'-q {queue}',
        f'-C', 'cpu',
        script_path,
        date_tag
    ]
    
    print('sbatch_cmd = ',sbatch_cmd) 

    print(f"Running: {' '.join(sbatch_cmd)}")
    
    try:
        result = subprocess.run(sbatch_cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        print(f"✓ Job submitted successfully")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error submitting job: {e.stderr}", file=sys.stderr)
        return 1




def main():
    parser = argparse.ArgumentParser(
        description='Submit analysis jobs to SLURM cluster',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 submit_slurm.py 09_10_2024                    # Array job mode
  python3 submit_slurm.py 09_10_2024 10 --parallel=4   # Array job with max 10 files
  python3 submit_slurm.py 04_16_2026 --single           # Single job: srun -n 1 postAna 04_16_2026
  python3 submit_slurm.py 04_16_2026 --single --time=02:00:00 --cpus-per-task=8 --mem=16G
        """
    )
    
    parser.add_argument('date_tag', help='Date tag to search for in filenames')
    parser.add_argument('max_files', nargs='?', type=int, help='Maximum number of files to process')
    parser.add_argument('--parallel', type=int, default=8, help='Number of parallel jobs (default: 8)')
    parser.add_argument('--time', dest='time_limit', default='01:00:00', 
                       help='Time limit in HH:MM:SS format (default: 03:00:00)')
    parser.add_argument('--single', action='store_true', 
                       help='Submit a single job (srun -n 1 postAna <date_tag>) instead of array job')
    parser.add_argument('--postAna-path', default='./compiled/postAna',
                       help='Path to postAna executable (default: ./compiled/postAna)')
    parser.add_argument('--job-name', help='Custom SLURM job name')
    parser.add_argument('--nodes', type=int, default=1, help='Number of nodes (default: 1)')
    parser.add_argument('--ntasks', type=int, default=1, help='Number of tasks (default: 1)')
    parser.add_argument('--cpus-per-task', type=int, default=4, help='CPUs per task (default: 4)')
    parser.add_argument('--mem', default='32G', help='Memory allocation (default: 32G)')
    parser.add_argument('-A', '--account', default='m2676', help='SLURM account (default: m2676)')
    parser.add_argument('-q', '--queue', default='shared', help='SLURM queue (default: shared)')  
    parser.add_argument('--library-path', default='/global/homes/m/mgold/mgold/bacon2Data/:$LD_LIBRARY_PATH',help='Library path to set LD_LIBRARY_PATH')
    parser.add_argument('--work-dir', default='/global/homes/m/mgold/mgold/bacon2Data', help='Working directory for the job')
    parser.add_argument('--output-file', help='Path to save the sbatch script (default: logs/postAna_<date_tag>.sh)')
    args = parser.parse_args()
    
    # Array job mode (original behavior)
    files = get_matching_files(args.date_tag)
    n = len(files)
    ntot = n

         #print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 2):
        n = int(args.max_files)

        print(" number of files to run  %i of %i  " % (n, ntot))
        

        #for i in range(0, n):
        #    print(" file ", i, " file ", files[i]) 
        
        
        return submit_slurm_job(args.date_tag, args.max_files, args.mem, args.parallel, args.time_limit, args.account, args.queue)


if __name__ == '__main__':
    sys.exit(main())
