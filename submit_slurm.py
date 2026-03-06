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
                files.append(os.path.join(rootdata_dir, filename))
    except FileNotFoundError:
        print(f"Error: Directory '{rootdata_dir}' not found", file=sys.stderr)
        return []
    
    return sorted(files)


def submit_slurm_job(date_tag, num_files=None, parallel_jobs=8, time_limit="01:00:00"):
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
        f'--array=0-{num_tasks-1}%{parallel_jobs}',
        f'--time={time_limit}',
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
  python3 submit_slurm.py 09_10_2024
  python3 submit_slurm.py 09_10_2024 10 --parallel=4
  python3 submit_slurm.py 09_10_2024 --time=02:00:00
        """
    )
    
    parser.add_argument('date_tag', help='Date tag to search for in filenames')
    parser.add_argument('max_files', nargs='?', type=int, help='Maximum number of files to process')
    parser.add_argument('--parallel', type=int, default=8, help='Number of parallel jobs (default: 8)')
    parser.add_argument('--time', dest='time_limit', default='01:00:00', 
                       help='Time limit in HH:MM:SS format (default: 01:00:00)')
    
    args = parser.parse_args()
    
    files = get_matching_files(args.date_tag)
    n = len(files)
    ntot = n

     #print(" files %i ", len(p), " files %i ", len(files))
    if (len(sys.argv) > 2):
        n = int(args.max_files)

    print(" number of files to run  %i of %i  " % (n, ntot))
    

    #for i in range(0, n):
    #    print(" file ", i, " file ", files[i]) 
    
    
    return submit_slurm_job(args.date_tag, args.max_files, args.parallel, args.time_limit)


if __name__ == '__main__':
    sys.exit(main())
