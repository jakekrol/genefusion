#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Aggregate gene burden across multiple tissue files')
    parser.add_argument('--files', required=True, help='Comma-separated list of input files')
    parser.add_argument('--output', help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    # Parse file list
    file_list = [f.strip() for f in args.files.split(',')]
    
    # Read and combine all files
    all_data = []
    for file_path in file_list:
        df = pd.read_csv(file_path, sep='\t', header=0)
        
        # Handle different column names (breast has 'burden_total_tumor', others have 'burden_total_dna_tumor')
        if 'burden_total_tumor' in df.columns:
            df = df.rename(columns={'burden_total_tumor': 'burden_total_dna_tumor'})
        
        all_data.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)

    
    # Group by gene and sum burden values
    result = combined_df.groupby('gene')['burden_total_dna_tumor'].sum().reset_index()
    
    # sort by burden descending
    result = result.sort_values(by='burden_total_dna_tumor', ascending=False)
    
    # Output
    if args.output:
        result.to_csv(args.output, sep='\t', index=False, header=True)
    else:
        result.to_csv(sys.stdout, sep='\t', index=False, header=False)
    

if __name__ == '__main__':
    main()
