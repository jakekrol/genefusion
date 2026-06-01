import subprocess
import shutil
import os
import pandas as pd

def fusiontsv2bedpe(path_fusion_tsv, path_bedpe_out, bgzip=False, pad=True, usecols=[1,2,7,9]):

    def breakpoint2interval(breakpoint):
        '''
        breakpoint: str
        star-fusion breakpoint strings are formatted as "chr:pos:strand"
        this function extracts the chrom, pos, and strand, converts to intervals and encodes strand as -1/1
        '''
        chrom, pos, strand = breakpoint.split(':')
        if chrom.startswith('chr'):
            chrom = chrom[3:]  # Remove 'chr' prefix if present
        chrom = str(chrom)
        start = int(pos) - 1
        end = int(pos)
        if strand == '-':
            strand = -1
        elif strand == '+':
            strand = 1
        else:
            raise ValueError(f"Invalid strand {strand} in breakpoint {breakpoint}")
        return chrom, start, end, strand
    
    df = pd.read_csv(path_fusion_tsv, sep='\t', usecols=usecols)
    df.columns = ['num_junction', 'num_spanning', 'left_breakpoint', 'right_breakpoint']
    df['total'] = df['num_junction'] + df['num_spanning']
    
    # Generate all output lines
    lines = []
    for i, row in df.iterrows():
        chromA, startA, endA, strandA = breakpoint2interval(row['left_breakpoint'])
        chromB, startB, endB, strandB = breakpoint2interval(row['right_breakpoint'])
        # write the interval to outfile as many times as the total support for that fusion
        for _ in range(row['total']):
            if pad:
                line = f"{chromA}\t{startA}\t{endA}\t{strandA}\t{chromB}\t{startB}\t{endB}\t{strandB}\t0\n"
            else:
                line = f"{chromA}\t{startA}\t{endA}\t{strandA}\t{chromB}\t{startB}\t{endB}\t{strandB}\n"
            lines.append(line)
    
    # Write output with optional bgzip compression
    if bgzip:
        bgzip_cmd = shutil.which('bgzip')
        if not bgzip_cmd:
            raise FileNotFoundError("bgzip command not found in PATH")
        if not path_bedpe_out.endswith('.gz'):
            raise ValueError(f"Output file {path_bedpe_out} must end with .gz when bgzip=True")
        
        # Pipe through bgzip
        process = subprocess.Popen([bgzip_cmd, '-c'], stdin=subprocess.PIPE, stdout=open(path_bedpe_out, 'wb'), stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(input=''.join(lines).encode())
        if process.returncode != 0:
            raise RuntimeError(f"bgzip failed: {stderr.decode()}")
    else:
        with open(path_bedpe_out, 'w') as f:
            f.writelines(lines)

    

def chimeric2bedpe(path_star_chimeric, path_bedpe_out, has_header=False, bgzip=False, pad=True):
    """
    Convert STAR Chimeric.out.junction format to sorted BEDPE format.
    
    Args:
        path_star_chimeric: Path to STAR Chimeric.out.junction file
        path_bedpe_out: Output BEDPE file path
        has_header: Whether the STAR chimeric file has a header line (default False)
    
    Input columns:
    - c1: chrA
    - c2: bpA (end)
    - c3: strand_donorA
    - c4: chrB
    - c5: bpB (end)
    - c6: strand_acceptorB
    - c11: start_alnA (start)
    - c13: start_alnB (start)
    
    Output: BEDPE format with intervals sorted lexicographically by chr, start, end
    """
    
    if not os.path.exists(path_star_chimeric):
        raise FileNotFoundError(f"Input file {path_star_chimeric} does not exist")
    
    awk_script = r"""
    BEGIN { OFS="\t" }
    NR==1 && has_header { next }
    /^#/ { next }
    {
        chrA = $1;    startA = $11;  endA = $2;   strandA = $3
        chrB = $4;    startB = $13;  endB = $5;   strandB = $6
        
        # Determine which interval goes "left" using cascading comparisons
        # Left = smaller chr lexicographically, or same chr with smaller start, or same chr/start with smaller end
        if (chrA < chrB || (chrA == chrB && startA < startB) || (chrA == chrB && startA == startB && endA < endB)) {
            # A is left, B is right
            if (pad) {
                print chrA, startA, endA, strandA, chrB, startB, endB, strandB, 0
            } else {
                print chrA, startA, endA, strandA, chrB, startB, endB, strandB
            }
        } else {
            # B is left, A is right
            if (pad) {
                print chrB, startB, endB, strandB, chrA, startA, endA, strandA, 0
            } else {
                print chrB, startB, endB, strandB, chrA, startA, endA, strandA
            }
        }
    }
    """
    
    if bgzip:
        bgzip_cmd = shutil.which('bgzip')
        if not bgzip_cmd:
            raise FileNotFoundError("bgzip command not found in PATH")
        assert path_bedpe_out.endswith('.gz'), f"Output file {path_bedpe_out} must end with .gz when bgzip=True"
        cmd = f"awk -v has_header={1 if has_header else 0} -v pad={1 if pad else 0} '{awk_script}' '{path_star_chimeric}' | sed 's|\t-\t|\t-1\t|g' | sed 's|\t+\t|\t1\t|g' | {bgzip_cmd} -c > '{path_bedpe_out}'"
    else:
        cmd = f"awk -v has_header={1 if has_header else 0} -v pad={1 if pad else 0} '{awk_script}' '{path_star_chimeric}' | sed 's|\t-\t|\t-1\t|g' | sed 's|\t+\t|\t1\t|g' > '{path_bedpe_out}'"
    
    try:
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"awk chimeric2bedpe failed: {e}")
        raise e
    
    return path_bedpe_out
