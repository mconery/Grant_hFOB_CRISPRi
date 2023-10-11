'''
download_roadmap_epigenomics_chip-seq_data.py

The purpose of this script is to take the list of roadmap epigenomics 
datasets, identify those for H3K27ac, H3K9ac, H3K4me1, and H3K4me3 with 
paired input runs, download them, and rename them to sensible names.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import os
import ast
import numpy as np
import pandas as pd

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -r => <str> Roadmap Epigenomics MetaData File REQUIRED
    -o => <str> Output directory for downloaded renamed files REQUIRED
    -a => <str> Annotation types OPTIONAL (default is ['H3K4me1','H3K4me3','H3K27ac','H3K9ac'])
                If used list must be input without spaces
    -f => <str> FTP Directory Address for file downloads OPTIONAL (default is 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak')
    -d => <int> Minimum acceptable number of reads OPTIONAL (default is 30M)
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "r:o:a:f:d:nh")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    ## Check for help flag and return empty string to break function
    help_ind = options_dict.get('-h', False)
    if help_ind != False:
        help()
    
    ## Required arguments
    try:
        roadmap_file = options_dict['-r']
        out_dir = options_dict['-o']        
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    data_types_str = options_dict.get('-a', "['H3K4me1','H3K4me3','H3K27ac','H3K9ac']")
    data_types = ast.literal_eval(data_types_str)
    ftp_address = options_dict.get('-f', "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak")
    read_depth = options_dict.get('-d', 3e7)
    
    #Check that directories and files exist
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
            print("WARNING: Output directory for peak files did not exist. It has been created.")
        except PermissionError:
            print("ERROR: Output directory for peak files does not exist and location is unwritable.")
            sys.exit(1)
    if exists(os.path.dirname(roadmap_file)) == False:
        print("ERROR: Invalid location specified for roadmap metadata file.")
        sys.exit(1)
    
    #Recast significance/purity thresholds and cutting height
    try:
        read_depth = int(read_depth)
    except ValueError:
        print("ERROR: Input minimum read count should be an integer.")
        sys.exit(1)
    
    print("Acceptable Inputs Given")
    #Call driver function
    driver(roadmap_file, out_dir, data_types, ftp_address, read_depth)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(roadmap_file, out_dir, data_types, ftp_address, read_depth):

    #Add slash to out_dir and ftp address if need be
    out_dir = add_slash(out_dir)
    ftp_address = add_slash(ftp_address)
    #Create a dummy scripts folder in the out_dir which will contain the scripts for downloading
    if not os.path.exists(out_dir + 'scripts'):
        os.makedirs(out_dir + 'scripts')
    
    #Read in file as a pandas dataframe
    roadmap_raw = pd.read_table(roadmap_file, sep=",")
    #Filter the data frame for the selected file types
    roadmap_filt = roadmap_raw[roadmap_raw['MARK'].isin(data_types)]
    #Filter for the minimum read count
    roadmap_filt = roadmap_filt[roadmap_filt['NREADS'] >= read_depth]
    
    #Count number of cell types (For testing)
    #len(roadmap_filt['E-Mnemonic'].unique())
    
    #Prepare the commands to download the resulting files
    chip_combos = roadmap_filt['E-Mnemonic'] + '.' + roadmap_filt['MARK'] 
    download_commands = 'wget -O ' + out_dir + chip_combos + '.narrowPeak.gz ' + ftp_address + roadmap_filt['FNAME'].str.replace(r'\*$', "narrowPeak.gz", regex=True)
    #Prepare names of script files
    script_names = out_dir + 'scripts/' + chip_combos + '.download.sh'
    #Iterate over the commands and send each to the console
    for i in range(len(chip_combos)):
        #Prepare variable names for the chip_combo
        command = download_commands.iloc[i]
        script_name = script_names.iloc[i]
        script_content = f"""#!/bin/bash\n{command}"""
        #Write to file
        with open(script_name, "w") as script_file:
            script_file.write(script_content)
        #Submit job
        os.system('sbatch ' + script_name)
    



###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
