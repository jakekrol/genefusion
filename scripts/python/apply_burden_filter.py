#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Filter table by gene blacklist')
    parser.add_argument('--table', required=True, help='Input table file')
    parser.add_argument('--genes', required=True, help='Gene blacklist file')
    parser.add_argument('--cx', type=int, required=True, help='1-indexed column for gene x')
    parser.add_argument('--cy', type=int, required=True, help='1-indexed column for gene y') 
    parser.add_argument('--output', required=True, help='Output file')
    
    args = parser.parse_args()
    
    # Read gene blacklist
    with open(args.genes) as f:
        # skip the first header line
        next(f)
        blacklist = set()
        for line in f:
            if line.strip():
                gene = line.strip().split('\t')[0]  # Take first column (gene name)
                blacklist.add(gene)
    
    # Convert to 0-indexed
    cx, cy = args.cx - 1, args.cy - 1
    
    # Filter table
    i=1
    with open(args.table) as infile, open(args.output, 'w') as outfile:
        header = next(infile)
        outfile.write(header)
        
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) > max(cx, cy):
                if (fields[cx] not in blacklist) and (fields[cy] not in blacklist):
                    outfile.write(line)
            print(f'Processed {i} lines')
            i += 1

if __name__ == '__main__':
    main()
