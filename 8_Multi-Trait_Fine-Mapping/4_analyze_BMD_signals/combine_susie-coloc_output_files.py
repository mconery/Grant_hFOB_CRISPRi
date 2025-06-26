'''
combine_susie_output_files.py

This script combines the locus-level result files from SuSiE-Coloc into a
single file. 

***This needs to be run in the rsparsepro conda environment ***
'''

import os
import sys
import glob
import pandas as pd

def combine_files(directory, naming_convention, output_file='combined_output.csv'):
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
        df = pd.read_csv(file)
        dfs.append(df)
    
    # Concatenate all dataframes, keeping the header only once
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save to a new file
    combined_df.to_csv(output_file, index=False)
    print(f"Combined file saved as {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <directory> <naming_convention>")
        sys.exit(1)
    directory = sys.argv[1]
    naming_convention = sys.argv[2]
    combine_files(directory, naming_convention)
