'''
finemap_SuSiE.py

This script will prepare to fine map a given locus using SuSiE. It will then 
call an R script that will complete the mapping. The 
output will be in a RDS standard format so that it can be used for running the
full SuSiE-Coloc pipeline.

***This needs to be run with the base python 3.9 install.                  ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
from os.path import exists
import os
from scipy.stats import t
import copy
import json
from functools import reduce
import subprocess

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> file prefix for files should consist of chr.start.stop REQUIRED
    -o => <txt> output folder location REQUIRED
    -m => <txt> munged summary statistics directory (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats) OPTIONAL
    -l => <txt> Storage location for map and matrix files (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/ld_matrices) OPTIONAL
    -n => <txt> Location of file of trait sample sizes (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes.tsv) OPTIONAL
    -k => <int> number of signals to map (default is 10) OPTIONAL
    -u => <num> purity threshold for determining real signals (default is 0.1) OPTIONAL
    -f => <num> confidence threshold (default is 0.95) OPTIONAL
    -j => <txt> location of json specifying traits to map at each locus (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json) OPTIONAL
    -t => <txt> directory for temp mapping files (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_mapping_files) OPTIONAL
    -i => <txt> directory for temp mapping scripts (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_script_files) OPTIONAL
    -r => <num> random seed (default is 5) OPTIONAL
    -c => <str> summary stats file naming convention showing where  TRAIT is found relative to periods 
                (default is TRAIT.sumstats.gz) OPTIONAL
    -e => <num> memory allocated to mapping jobs in GB (default is 40GB) OPTIONAL
""")
    sys.exit(exit_num)


###############################################################################
#################################  TEST VARIABLES  ############################
###############################################################################

file_prefix="chr19.44000001.44750000"
out_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/susie_results"
munge_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats"
ld_loc="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/ld_matrices"
size_file="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes_cc.tsv"
num_signals = 10
purity_thresh = 0.1
confidence = 0.95
trait_json = "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json"
temp_dir = "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_mapping_files"
script_dir = "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_script_files"
random_seed = 5
file_format = 'TRAIT.sumstats.gz'
map_script = "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi/8_Multi-Trait_Fine-Mapping/3_fine-mapping/finemap_SuSiE.R"
susie_mem = 40

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:m:l:n:k:u:f:j:t:i:r:c:e:nh")
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
        file_prefix = options_dict['-p']
        out_dir = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Optional arguments
    munge_dir = options_dict.get('-m', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats")
    ld_loc = options_dict.get('-l', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/ld_matrices")
    size_file = options_dict.get('-n', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes.tsv")
    num_signals = options_dict.get('-k', 10)
    purity_thresh = options_dict.get('-u', 0.1)
    confidence = options_dict.get('-f', 0.95)
    trait_json = options_dict.get('-j', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json")
    temp_dir = options_dict.get('-t', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_mapping_files")
    script_dir = options_dict.get('-i', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/temp_script_files")
    random_seed = options_dict.get('-r', 5)
    file_format = options_dict.get('-c', 'TRAIT.sumstats.gz')
    susie_mem = options_dict.get('-e', 40)
    #Confirm that the file/folder locs exist
    if exists(munge_dir) == False:
        print("ERROR: Directory for munged summary statistics does not exist.")
        sys.exit(1)
    if exists(ld_loc) == False:
        print("ERROR: Storage location for LD matrices and map files does not exist.")
        sys.exit(1)
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(temp_dir) == False:
        print("ERROR: Temp directory for mapping files does not exist.")
        sys.exit(1)
    if exists(script_dir) == False:
        print("ERROR: Temp directory for script files does not exist.")
        sys.exit(1)
    if exists(size_file) == False:
        print("ERROR: File of sample sizes does not exist at specified location.")
        sys.exit(1)
    if exists(trait_json) == False:
        print("ERROR: Trait json does not exist at specified location.")
        sys.exit(1)
    #Recast random_seed 
    try:
        random_seed = int(random_seed)
    except ValueError:
        print("ERROR: Random seed must be an integer.")
        sys.exit(1)
    #Recast signal count
    try:
        num_signals = int(num_signals)
        if num_signals < 1:
            print("ERROR: The number of mapped signals must be >= 1.")
    except ValueError:
        print("ERROR: The number of mapped signals must be an integer >= 1.")
        sys.exit(1)
    #Recast susie mapping memory (GB)
    try:
        susie_mem = int(susie_mem)
        if susie_mem < 1:
            print("ERROR: The mapping memory for each SuSiE run must be >= 1GB.")
    except ValueError:
        print("ERROR: The mapping memory for each SuSiE run must be an integer value >=1 GB.")
        sys.exit(1)

    #Verify that TRAIT is present in file format
    if 'TRAIT' not in file_format:
        print("ERROR: Unable to identify position of TRAIT in file name format.")
        sys.exit(1)
    #Recast purity value 
    try:
        purity_thresh = float(purity_thresh)
        if purity_thresh >= 1 or purity_thresh < 0: 
            print("ERROR: Purity level should be in [0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Purity level is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Recast confidence value 
    try:
        confidence = float(confidence)
        if confidence >= 1 or confidence <= 0: 
            print("ERROR: Confidence level should be in (0,1).")
            sys.exit(1)
    except ValueError:
        print("ERROR: Confidence level is not coercible to a float. It should be a proportion in (0,1).")
        sys.exit(1)
        
    print("Acceptable Inputs Given for " + file_prefix)
    
    #Call driver function
    try:
        driver(file_prefix, out_dir, munge_dir, ld_loc, size_file, num_signals, purity_thresh, confidence, trait_json, temp_dir, script_dir, random_seed, file_format, susie_mem)
    except MemoryError:
        print("ERROR: Insufficient memory for fine-mapping job. Try rerunning with additional memory")
        sys.exit(1)

###############################################################################
#############################  HELPFUL FUNCTIONS ##############################
###############################################################################

#Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

#Define a function to write shell scripts
def write_shell(shell_script, job_mem, job_commands, dependencies = []):
    with open(shell_script, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -J ' + os.path.basename(shell_script) + '\n')
        f.write('#SBATCH -o ' + shell_script.replace('.sh', '.log') + '\n')
        f.write('#SBATCH --mem ' + str(job_mem + 1) + 'G\n')
        if dependencies != []:
            f.write('#SBATCH -d afterok:' + ':'.join(dependencies) + '\n')
        f.write('#SBATCH --time=2-00:00:00\n')
        for command in job_commands:
            f.write(command + '\n')
    #Set file permissions
    os.chmod(shell_script, 0o755)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(file_prefix, out_dir, munge_dir, ld_loc, size_file, num_signals, purity_thresh, confidence, trait_json, temp_dir, script_dir, random_seed, file_format, susie_mem): 
    
    #Set fine-mapping script paths
    temp = os.path.realpath(__file__)
    map_script = temp[:temp.find('finemap_SuSiE.py')] + 'finemap_SuSiE.R'
    
    #Add final slash to out_dir and ld_loc if it's missing
    out_dir = add_slash(out_dir)
    ld_loc = add_slash(ld_loc)
    munge_dir = add_slash(munge_dir)
    temp_dir = add_slash(temp_dir)
    script_dir = add_slash(script_dir)
    
    #Extract info from the file_prefix
    temp = file_prefix.split(".")
    chromo = int(temp[0][3:])
    bp_start =  int(temp[1])
    bp_end = int(temp[2])
        
    #Get needed files from file prefix
    try:
        int_matrix_name = ld_loc + [x for x in os.listdir(ld_loc) if x.startswith(file_prefix) and x.endswith('.ld.gz')][0]
        int_map_name = ld_loc + [x for x in os.listdir(ld_loc) if x.startswith(file_prefix) and x.endswith('.map')][0]
    except IndexError:
        print("ERROR: Either or both the map and ld matrix file are mssing for " + file_prefix + '.')
        sys.exit(1)
        
    #Load in trait json
    with open(trait_json, 'r') as json_file:
        trait_dict = json.load(json_file)
    #Get traits for the locus
    try:
        locus_traits=trait_dict[file_prefix]
    except KeyError:
        print("ERROR: Given locus missing from trait json. Check json location or file prefix: " + file_prefix + '.')
        sys.exit(1)
    
    #Read in lists of traits and filter for SNPs in locus
    trait_stats_dict = {}
    for trait in locus_traits:
        trait_stats_loc = munge_dir + file_format.replace('TRAIT', trait)
        try:
            trait_stats_dict[trait] = pd.read_table(trait_stats_loc, sep="\t")
            trait_stats_dict[trait] = trait_stats_dict[trait][(trait_stats_dict[trait]['CHR'] == chromo) & (trait_stats_dict[trait]['BP'] >= bp_start) & (trait_stats_dict[trait]['BP'] <= bp_end)]
        except FileNotFoundError:
            print("ERROR: Munged summary stats not found for " + trait + " in " + munge_dir + '.')
            sys.exit(1)
    
    #Read in the map file
    map_raw = pd.read_table(int_map_name, sep="\t", header=None)
    map_raw.columns = pd.Series(['CHR', 'SNP', 'CM', 'BP'])
    #Make a list of all the variant ids and filter for the lowest common denominator
    variant_ids = [trait_stats_dict[trait]['SNP'] for trait in trait_stats_dict.keys()]
    variant_ids.append(map_raw.iloc[:,1])
    intersected_variants = reduce(lambda x, y: pd.merge(x, y, how='inner'), variant_ids)
    intersected_variants = intersected_variants['SNP']
    
    #Filter the map file for the intersected SNPs
    map_filt = map_raw.loc[map_raw['SNP'].isin(intersected_variants)]
    #Read in the ld matrix and filter it for intersected SNPs
    ld_raw = np.loadtxt(int_matrix_name, dtype=float,  delimiter=",")
    ld_filt = ld_raw[map_filt.index,:]
    ld_filt = ld_filt[:,map_filt.index]
    #Reduce pull on memory
    del ld_raw
    #Make the matrix symmetric
    ld_symmetric = ld_filt + ld_filt.T - np.diag(ld_filt.diagonal())
    #Write ld file to file
    ld_filt_df = pd.DataFrame(ld_symmetric, columns = map_filt['SNP'])
    ld_filt_df.index = ld_filt_df.columns
    ld_temp_loc = temp_dir + file_prefix + '.ld.gz'
    ld_filt_df.to_csv(temp_dir + file_prefix + '.ld.gz', sep = "\t", header =True, index = True, na_rep='nan', compression='gzip')
    #Reduce memory pull
    del ld_filt
    
    #Create a list of job ids for the deletion job
    trait_jobids = []
    #Loop over traits at the locus and make temp files of the summary stats then process jobs to console
    for trait in locus_traits:
        #Write data to file
        temp_stats = trait_stats_dict[trait]
        temp_stats = temp_stats.loc[temp_stats['SNP'].isin(intersected_variants)]
        temp_stats_loc = temp_dir + trait + '.' + file_prefix + '.tsv.gz'
        temp_stats.to_csv(temp_stats_loc, sep = "\t", header =True, index = False, na_rep='nan', compression='gzip')
        #Create SuSiE out file name
        out_file_loc = out_dir + trait + "." + file_prefix + ".susie.rds"
        #Create script file
        trait_commands = ["Rscript " + map_script + " " + temp_stats_loc + " " + ld_temp_loc + " " + size_file + " " + out_file_loc + " " + str(random_seed) + " " + str(num_signals) + " " + str(purity_thresh) + " " + str(confidence)]
        trait_commands.extend(["rm " + temp_stats_loc])
        trait_script = script_dir + file_prefix + '.' + trait + '.mapping.sh'
        write_shell(trait_script, susie_mem, trait_commands)
        #Submit the job
        temp = subprocess.run(['sbatch', trait_script], stdout=subprocess.PIPE, text=True)    
        #Capture the job id
        trait_jobids.append(temp.stdout.strip().split(" ")[3])
    
    #Submit deletion script for the ld matrices
    ld_del_command = ["rm " + ld_temp_loc]
    ld_del_script = script_dir + file_prefix + '.ld_del.sh'
    write_shell(ld_del_script, 16, ld_del_command, dependencies = trait_jobids)
    temp = subprocess.run(['sbatch', ld_del_script], stdout=subprocess.PIPE, text=True)
    
        
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
