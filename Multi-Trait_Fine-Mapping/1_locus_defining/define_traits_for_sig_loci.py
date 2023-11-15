'''
define_traits_for_sig_loci.py

The purpose of this script is to read over the list of significant loci defined
for BMD multi-trait fine-mapping and identify the traits that have a minimal 
level of significance at each of the loci and hence should be included in the
fine-mapping. 

***This needs to be run on the default Respulbica Python/3.9 as there's  ***
***an incompatability with other versions of python.                     ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
import os
import pandas as pd
import numpy as np
from os.path import exists
from scipy.stats import norm
import json

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -l => <txt> locus file directory REQUIRED
    -g => <txt> munged summary statistics file directory  REQUIRED
    -p => <flt> p-value significance threshold (default is 1e-6) OPTIONAL
    -n => <txt> variant file naming convention (default is CHR.START.END.variants.csv) OPTIONAL
    -a => <txt> gwas summary statistics naming convention (default is TRAIT.sumstats.gz) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:g:p:n:a:nh")
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
        locus_dir = options_dict['-l']
        munge_dir = options_dict['-g']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    p_thresh = options_dict.get('-p', 1e-6)
    variant_name_conv = options_dict.get('-n', 'CHR.START.END.variants.csv')
    gwas_name_conv = options_dict.get('-a', 'TRAIT.sumstats.gz')

    #Confirm that float was given for optional p-value threshold parameter
    try:
        p_thresh = float(p_thresh)
    except ValueError:
        print("ERROR: P-value threshold must be a float.")
        sys.exit(1)
    #Verify that the p-value threshold is between 0 and 1
    if p_thresh <= 0 or p_thresh > 1:
        print("ERROR: p-value threshold must be in (0,1].")
        sys.exit(1)
    #Confirm that the directories exist
    if exists(locus_dir) == False:
        print("ERROR: Input locus directory does not exist.")
        sys.exit(1)
    if exists(munge_dir) == False:
        print("ERROR: Input directory of munged GWAS summary statistics does not exist.")
        sys.exit(1)
    #Verify that TRAIT is present in gwas naming convention
    if 'TRAIT' not in gwas_name_conv:
        print("ERROR: Unable to identify position of TRAIT in gwas file name convention.")
        sys.exit(1)
    #Also check the variant naming convention
    if 'CHR' not in variant_name_conv:
        print("ERROR: Unable to identify position of CHRomosome in variant file name convention.")
        sys.exit(1)
    elif 'START' not in variant_name_conv or 'END' not in variant_name_conv:
        print("ERROR: Unable to identify position of locus START or END in variant file name convention.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(locus_dir, munge_dir, p_thresh, variant_name_conv, gwas_name_conv)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

#Define a function that checks a name conv
def check_file_name(file_name, name_conv):
    conv_temp = name_conv.split(".")
    file_temp = file_name.split(".")
    if len(conv_temp) != len(file_temp):
        return(False)
    temp = [True for x in range(len(conv_temp)) if conv_temp[x] == "TRAIT" or conv_temp[x] == file_temp[x]]
    return all(temp)

###############################################################################
################################  TEST LOCATIONS ##############################
###############################################################################

locus_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files"
munge_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats"
p_thresh=1e-6
variant_name_conv="CHR.START.END.variants.csv"
gwas_name_conv="TRAIT.sumstats.gz"

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(locus_dir, munge_dir, p_thresh, variant_name_conv, gwas_name_conv):
    
    #Add slash to directories if needed
    munge_dir = add_slash(munge_dir)
    locus_dir = add_slash(locus_dir)
    #List all variant files in the gwas dir and get list of traits
    trait_files = os.listdir(munge_dir)
    trait_files = [x for x in trait_files if check_file_name(x, gwas_name_conv)]
    traits = [x.split(".")[gwas_name_conv.split(".").index("TRAIT")] for x in trait_files]
    #List all variant files in the locus dir and get list of loci
    variant_files = os.listdir(locus_dir)
    variant_files = [x for x in variant_files if check_file_name(x, variant_name_conv)]
    temp = variant_name_conv.split(".")
    loci = [str(x.split(".")[temp.index("CHR")]) + "." + str(x.split(".")[temp.index("START")]) + "." + str(x.split(".")[temp.index("END")]) for x in variant_files]
    
    #Calculate p-value threshold in Negative log terms
    neg_log_p_thresh = -np.log10(p_thresh)
    #Make a dictionary to hold the traits for each locus
    out_dict = dict(zip(loci, [[] for x in loci]))
    
    #Loop over the traits and then over the loci to fill the dictionary
    for trait_index in range(len(trait_files)):
        #Read in the full gwas for the current trait
        trait_gwas = pd.read_csv(munge_dir + trait_files[trait_index], sep="\t") #Fixing index column to format of munged file
        trait_gwas.set_index('SNP', inplace=True)
        #Use z-scores to calculate -logP values
        trait_gwas['-logP'] = (-np.log(2) - norm.logcdf(-np.abs(trait_gwas['Z'])))/np.log(10)
        #Filter for rows that are significant at the threshold
        filtered_trait_gwas = trait_gwas.loc[trait_gwas['-logP'] > neg_log_p_thresh]
        #Now loop over the loci and find any that contain one of the remaining filtered variants
        for locus_index in range(len(variant_files)):
            #Read in locus file
            locus_variants = pd.read_csv(locus_dir + variant_files[locus_index], index_col=0, sep=",")
            #Attempt to find significant variants in the locus by merging pandas dfs
            merged_df = pd.merge(filtered_trait_gwas, locus_variants, left_index=True, right_index=True, how='inner')
            #Check if any rows remained and if so then append the trait to the locus entry in the dictionary
            if len(merged_df) > 0:
                out_dict[loci[locus_index]].append(traits[trait_index])
        #Print update
        print("NOTE: Finished identifying loci for " + traits[trait_index])
        
    
    #Create trait_file_name
    output_file = locus_dir + "traits_per_loci.json"
    #Write to file
    with open(output_file, 'w') as json_file:
        json.dump(out_dict, json_file)
        

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
