#!/usr/bin/env python3

import argparse
import sys, os, shutil
import tempfile
import duckdb as ddb

parser = argparse.ArgumentParser(description='Convert fusion tsv to bedpe format')
parser.add_argument('-i', '--input', type=str, required=True, help='Input fusion tsv file')
parser.add_argument('-b', '--bed', type=str, required=True, help='Input bed file for gene coordinates')
parser.add_argument('-o', '--output', type=str, required=True, help='Output bedpe file')
parser.add_argument('--temp-dir', default='/data/jake/tmp', type=str, help='Temporary directory for intermediate files')
args = parser.parse_args()

def main():
    # Use temporary directory for any intermediate files if needed
    temp_dir = args.temp_dir or tempfile.gettempdir()
    temp_file = os.path.join(temp_dir, f"temp_joined_{os.getpid()}.tsv")
    
    try:
        # Step 1: Do the joins and write to temporary file
        with ddb.connect() as con:
            print('Step 1: Performing joins...')
            
            # Simplified query - just do the joins, no complex column reordering
            join_query = f"""
            COPY (
                WITH joined_left AS (
                    SELECT 
                        a.*,
                        b.chrom as chrom_left, 
                        b.start as start_left, 
                        b.end as end_left, 
                        b.strand as strand_left,
                        b.name as name_left
                    FROM read_csv_auto('{args.input}', delim='\t', header=true) a
                    LEFT JOIN read_csv_auto('{args.bed}', delim='\t', header=false, 
                                           names=['chrom','start','end','name','strand'],
                                           types={{'chrom': 'VARCHAR', 'start': 'INTEGER', 'end': 'INTEGER', 'name': 'VARCHAR', 'strand': 'VARCHAR'}}) b
                    ON a."left" = b.name
                ),
                final_joined AS (
                    SELECT 
                        a.*,
                        b.chrom as chrom_right, 
                        b.start as start_right, 
                        b.end as end_right, 
                        b.strand as strand_right,
                        b.name as name_right
                    FROM joined_left a
                    LEFT JOIN read_csv_auto('{args.bed}', delim='\t', header=false,
                                           names=['chrom','start','end','name','strand'],
                                           types={{'chrom': 'VARCHAR', 'start': 'INTEGER', 'end': 'INTEGER', 'name': 'VARCHAR', 'strand': 'VARCHAR'}}) b
                    ON a."right" = b.name
                )
                SELECT * FROM final_joined
            ) TO '{temp_file}' (DELIMITER '\t', HEADER true)
            """
            
            con.execute(join_query)
        
        # Move temp file to final output - column reordering can be done externally if needed
        print('Moving result to final output...')
        shutil.move(temp_file, args.output)
        temp_file = None  # Prevent cleanup since we moved it
    
    finally:
        # Clean up temporary file only if it still exists
        if temp_file and os.path.exists(temp_file):
            os.remove(temp_file)
    
    print(f'Successfully created joined file: {args.output}')
    print('Note: Use cut or similar tools to reorder columns to BEDPE format if needed.')

if __name__ == '__main__':
    main()

