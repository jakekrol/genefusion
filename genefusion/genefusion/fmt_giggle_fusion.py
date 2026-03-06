#!/usr/bin/env python3
"""
Format Giggle query output for gene fusion analysis.
Adds header with metadata, column names, and extracts sample basenames.
Outputs TSV format (optionally gzipped).
"""

import sys,os
import gzip
import argparse
from datetime import datetime
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Format Giggle output for gene fusion analysis'
    )
    parser.add_argument(
        '-g', '--gene',
        required=True,
        help='Gene name'
    )
    parser.add_argument(
        '-r', '--region',
        required=True,
        help='Query/gene region (e.g., chr1:13386750-13389469)'
    )
    parser.add_argument(
        '-m', '--metadata',
        default='',
        help='Additional metadata (e.g., tissue type like "kidney")'
    )
    parser.add_argument(
        '-i', '--input',
        default=None,
        help='Input file path (default: read from stdin)'
    )
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output file path (default: write to stdout)'
    )
    parser.add_argument(
        '-z', '--gzip',
        action='store_true',
        help='Compress output with gzip'
    )
    return parser.parse_args()


def extract_basename(sample_path):
    """Extract sample basename from full path.
    
    Example: beds/FI12751.excord.bed.gz -> FI12751.excord.bed.gz
    """
    base = os.path.basename(sample_path)
    return base


def main():
    args = parse_args()
    
    # Read from input file or stdin
    if args.input:
        with open(args.input, 'r') as f:
            lines = f.readlines()
    else:
        lines = sys.stdin.readlines()
    
    # Determine output destination and mode
    if args.output:
        # Write to file
        output_file = args.output
        if args.gzip and not output_file.endswith('.gz'):
            output_file += '.gz'
        opener = gzip.open if args.gzip else open
        mode = 'wt' if args.gzip else 'w'
        out = opener(output_file, mode)
        should_close = True
    else:
        # Write to stdout
        out = sys.stdout
        should_close = False
    
    try:
        # Write header
        out.write(f"# Gene={args.gene}\n")
        out.write(f"# Region={args.region}\n")
        out.write(f"# Date={datetime.now().strftime('%Y-%m-%d')}\n")
        if args.metadata:
            out.write(f"# Metadata={args.metadata}\n")
        
        # Write column names
        columns = [
            'chrm_left', 'start_left', 'end_left', 'strand_left',
            'chrm_right', 'start_right', 'end_right', 'strand_right', 'evidence_type',
            'sample'
        ]
        out.write('\t'.join(columns) + '\n')
        
        # Process and write data lines
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) >= 10:
                # Extract sample basename (10th field, index 9)
                fields[9] = extract_basename(fields[9])
                out.write('\t'.join(fields) + '\n')
    except BrokenPipeError:
        # Handle pipe being closed (e.g., when piping to head)
        # Silence the error and exit cleanly
        sys.stderr.close()
    finally:
        if should_close:
            out.close()
            print(f"Output written to {output_file}", file=sys.stderr)


if __name__ == '__main__':
    main()