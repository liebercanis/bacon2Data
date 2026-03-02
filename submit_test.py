


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