'''
combine_susie-coloc_output_files.py

This script combines the locus-level result files from SuSiE-Coloc into a
single file. 

***This needs to be run in the default conda environment ***
'''

import os
import glob
import pandas as pd
import argparse

def combine_files(directory, naming_convention, output_file, set_rows):
    #Convert set_rows to boolean
    if set_rows == 'False' or set_rows == 'F':
        set_rows = False
    else:
        set_row = True
    # Create the full path pattern for glob
    pattern = os.path.join(directory, naming_convention)
    # Find all files matching the pattern
    files = glob.glob(pattern)
    
    if not files:
        print("No files found matching the naming convention.")
        return
    
    # Initialize an empty list to hold dataframes
    dfs = []
    
    for file in files:
        try:
            df = pd.read_csv(file, delimiter="\t")
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue
    
    if not dfs:
        print("No dataframes were loaded. Exiting.")
        return
    
    # Concatenate all dataframes, keeping the header only once
    combined_df = pd.concat(dfs, ignore_index=(not set_rows))
    
    # Save to a new file
    combined_df.to_csv(output_file, index=set_rows, sep="\t")
    print(f"Combined file saved as {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Combine all files in a directory matching a naming convention into a single CSV file.',
        usage='python script.py <--directory> <--naming_convention> [--output OUTPUT_FILE]'
    )
    parser.add_argument('--directory', help='Directory containing files to combine')
    parser.add_argument('--naming_convention', help='File naming convention to match (e.g., *.csv)')
    parser.add_argument('--output', '-o', default='combined_output.csv', help='Name of the output file (default: combined_output.csv)')
    parser.add_argument('--set_rows', default='False', help='Whether to add row names to output file (default: False')
    args = parser.parse_args()
    
    combine_files(args.directory, args.naming_convention, args.output, args.set_rows)

if __name__ == '__main__':
    main()

