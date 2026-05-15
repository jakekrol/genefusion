import subprocess
import shutil
import os


def sort_bed(input_path, output_path, tmp_dir='/tmp', threads=4, bgzip=False):
    """
    Sort a BED/BEDPE file using external sorting with parallel merge and temporary directory.
    
    Args:
        input_path: Path to input BED/BEDPE file
        output_path: Path to output file
        tmp_dir: Temporary directory for sort scratch space (should be on fast storage)
        threads: Number of parallel sort threads
        bgzip: If True, compress output with bgzip (output path should end with .gz)
    
    Returns:
        output_path on success
    """
    
    # Validate input
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file {input_path} does not exist")
    
    if not os.path.exists(tmp_dir):
        raise FileNotFoundError(f"Temporary directory {tmp_dir} does not exist")
    
    if not os.path.isdir(tmp_dir):
        raise NotADirectoryError(f"{tmp_dir} is not a directory")
    
    # Decompress if needed (use bgzip for speed)
    if input_path.endswith('.gz'):
        bgzip_cmd = shutil.which('bgzip')
        if not bgzip_cmd:
            raise FileNotFoundError("bgzip command not found in PATH")
        decomp_cmd = f"{bgzip_cmd} -d -c '{input_path}'"
    else:
        decomp_cmd = f"cat '{input_path}'"
    
    # Sort command with external merge using temp dir and parallelization
    # LC_ALL=C for faster sorting, -S for buffer size, -T for temp dir, -P for parallel threads
    sort_cmd = f"LC_ALL=C sort --buffer-size=3G -k1,1 -k2,2n -k3,3n -T '{tmp_dir}' --parallel={threads}"
    
    # Optionally compress
    if bgzip:
        bgzip_cmd = shutil.which('bgzip')
        if not bgzip_cmd:
            raise FileNotFoundError("bgzip command not found in PATH")
        if not output_path.endswith('.gz'):
            raise ValueError(f"Output file {output_path} must end with .gz when bgzip=True")
        compress_cmd = f"{bgzip_cmd} -c"
        full_cmd = f"{decomp_cmd} | {sort_cmd} | {compress_cmd} > '{output_path}'"
    else:
        full_cmd = f"{decomp_cmd} | {sort_cmd} > '{output_path}'"
    
    try:
        subprocess.run(full_cmd, shell=True, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"sort_bed failed: {e}")
        raise e
    
    return output_path
