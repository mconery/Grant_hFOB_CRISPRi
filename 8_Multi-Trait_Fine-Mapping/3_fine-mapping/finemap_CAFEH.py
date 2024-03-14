'''
finemap_CAFEH.py

This script that will fine map a given locus using CAFEH. It will return as its
output a .pkl pickle file contains the default output of fine mapping plus the
negative log P-values, the base pair positions of the SNPs, and the rsids. The 
output will be in a standard format so that it can be used for plotting 
results, testing the effects of different purity/association thresholds, and 
any other downstream analyses that arise.

***This needs to be run in the cafeh conda environment.                     ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
from os.path import exists
import os
from cafeh.cafeh_summary import fit_cafeh_summary
from scipy.stats import t
import copy
import pickle
import json
from functools import reduce

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> file prefix for files should consist of trait.chr.start.stop.ancestry REQUIRED
    -o => <txt> output folder location REQUIRED
    -m => <txt> munged summary statistics directory (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats) OPTIONAL
    -l => <txt> Storage location for map and matrix files (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/ld_matrices) OPTIONAL
    -n => <txt> Location of file of trait sample sizes (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes.tsv) OPTIONAL
    -k => <int> starting number of components to map (default is 10) OPTIONAL
    -x => <int> max number of components to map (default is 30) OPTIONAL
    -u => <num> purity threshold for determining real signals (default is 0.1) OPTIONAL
    -j => <txt> location of json specifying traits to map at each locus (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json) OPTIONAL
    -r => <num> random seed (default is 5) OPTIONAL
    -c => <str> summary stats file naming convention showing where ANCestry and TRAIT are found relative to periods 
                (default is TRAIT.sumstats.gz) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:m:l:n:k:x:u:j:r:c:nh")
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
    starting_signals = options_dict.get('-k', 10)
    max_signals = options_dict.get('-x', 25)
    purity_thresh = options_dict.get('-u', 0.1)
    trait_json = options_dict.get('-j', "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json")
    random_seed = options_dict.get('-r', 5)
    file_format = options_dict.get('-c', 'TRAIT.sumstats.gz')
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
    #Recast starting signal count
    try:
        starting_signals = int(starting_signals)
    except ValueError:
        print("ERROR: The number of starting signals must be an integer.")
        sys.exit(1)
    try:
        max_signals = int(max_signals)
    except ValueError:
        print("ERROR: The maximum number of signals must be an integer.")
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
        
    print("Acceptable Inputs Given for " + file_prefix)
    
    #Call driver function
    try:
        driver(file_prefix, out_dir, munge_dir, ld_loc, size_file, starting_signals, max_signals, purity_thresh, trait_json, random_seed, file_format)
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

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(file_prefix, out_dir, munge_dir, ld_loc, size_file, starting_signals, max_signals, purity_thresh, trait_json, random_seed, file_format): 
    
    #Add final slash to out_dir and ld_loc if it's missing
    out_dir = add_slash(out_dir)
    ld_loc = add_slash(ld_loc)
    munge_dir = add_slash(munge_dir)
    
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
    #Filter the trait files for the intersected variants list
    trait_stats_dict = {trait:trait_stats_dict[trait][trait_stats_dict[trait]['SNP'].isin(intersected_variants)] for trait in trait_stats_dict.keys()}
    #Reset the indices on all the dataframes
    for trait in trait_stats_dict.keys():
        trait_stats_dict[trait].reset_index(drop=True, inplace=True)
    
    #Prepare the beta inputs to CAFEH 
    beta_df = pd.concat([trait_stats_dict[trait]['BETA'] for trait in trait_stats_dict.keys()], axis=1) #This code only works because the files were pre-sorted in ascending order by bp
    beta_df = beta_df.T
    beta_df.index = trait_stats_dict.keys()
    beta_df.columns = intersected_variants
    #Prepare the se inputs to CAFEH
    stderr_df = pd.concat([trait_stats_dict[trait]['SE'] for trait in trait_stats_dict.keys()], axis=1) #This code only works because the files were pre-sorted in ascending order by bp
    stderr_df = stderr_df.T
    stderr_df.index = trait_stats_dict.keys()
    stderr_df.columns = intersected_variants
    
    #Read in file of sample sizes and get sample sizes for each trait
    size_raw = pd.read_table(size_file, sep="\t", index_col=0)
    n_df = pd.DataFrame(np.ones_like(beta_df))
    n_df.index = beta_df.index
    n_df.columns = beta_df.columns
    n_df = n_df.T.mul(size_raw.loc[beta_df.index,'sample size']).T
    
    #Identify the variants in the map file that are good to use
    map_raw['in_munge'] = map_raw.iloc[:,1].isin(intersected_variants)
    #Read in the matrix
    ld_raw = np.loadtxt(int_matrix_name, dtype=float, delimiter=",")
    #Filter down the ld matrix and map file
    good_map_idx = map_raw[map_raw['in_munge'] == True].index
    ld_filt = ld_raw[list(good_map_idx)]
    del ld_raw #Reduce pull on memory
    ld_filt = ld_filt[:,list(good_map_idx)]
    map_filt = map_raw.iloc[good_map_idx, ]
    #Reset the pandas indices
    map_filt.reset_index(inplace = True, drop = True)
    #Make the matrix symmetric
    ld_symmetric = ld_filt + ld_filt.T - np.diag(ld_filt.diagonal())
    #Convert numpy array to pandas dataframe
    LD_df = pd.DataFrame(ld_symmetric, columns = map_filt.iloc[:,1])
    LD_df.index = LD_df.columns
    
    ##Fit CAFEH with betas and standard errors
    np.random.seed(random_seed) # Set random seed to ensure reproducible results
    cafehs = fit_cafeh_summary(LD_df, beta_df, stderr_df, n=n_df, K=starting_signals)
    #Verify that at least one low-purity signal was found
    K = starting_signals
    while np.sum([x < purity_thresh for x in list(cafehs.purity.values())]) == 0 and K < max_signals:
        K = K + 1 
        cafehs = fit_cafeh_summary(LD_df, beta_df, stderr_df, n=n_df, K=K)
    
    #Store extra basic data in the cafehs object
    cafehs.ID = np.array(map_filt['SNP'])
    cafehs.bp = np.array(map_filt['BP'])
    cafehs.chr = np.array(map_filt['CHR'])
    cafehs.rsid = np.array(list(trait_stats_dict.values())[0]['RSID'])
    cafehs.neglogP = np.array(pd.concat([pd.Series((-np.log(2) - t.logcdf(-np.abs(trait_stats_dict[trait]['BETA']/trait_stats_dict[trait]['SE']), df=n_df.loc[trait][0]-1))/np.log(10)) for trait in trait_stats_dict.keys()], axis=1).T)
    cafehs.n = np.array(n_df) #This works since all SNPs have same sample size for each trait
    
    #Hard print purity in the object as well
    cafehs.realpure = copy.deepcopy(cafehs.purity) #Need to prevent recalculation every time!!!
    #Hard print weight_means in as well
    cafehs.realweight_means = copy.deepcopy(cafehs.weight_means)
    
    ##### Calculate residuals and add them to the cafehs object #####
    for i in range(K):
        #Calculate test statistics for preceding residuals
        if i > 0:
            precede_prediction = np.sum([cafehs.compute_first_moment(l) for l in range(cafehs.dims['K']) if l < i], axis=0)
            precede_test = cafehs.B - precede_prediction #Residuals when predicting betas given the preceding components
            precede_test2 = precede_test/cafehs.S #Estimated Z-scores when predicting given the preceding components
        else:
            precede_test2 = cafehs.B/cafehs.S #Actual Z-scores when looking at the first component
        #Calculate test statistics for absolute residuals
        abs_prediction = np.sum([cafehs.compute_first_moment(l) for l in range(cafehs.dims['K']) if l != i], axis=0)
        abs_test = cafehs.B - abs_prediction #Residuals when predicting betas given all other components
        abs_test2 = abs_test/cafehs.S #Estimated Z-scores when predicting given all other components
        #Calculate pvals for the residuals
        abs_neglogp_residual = (-np.log(2) - t.logcdf(-np.abs(abs_test2), df=cafehs.n-1))/np.log(10)
        precede_neglogp_residual = (-np.log(2) - t.logcdf(-np.abs(precede_test2), df=cafehs.n-1))/np.log(10)
        if i == 0:
            cafehs.abs_neglogp_resid = {cafehs.study_ids[i]:abs_neglogp_residual[i] for i in range(len(cafehs.study_ids))}
            cafehs.precede_neglogp_resid = {cafehs.study_ids[i]:precede_neglogp_residual[i] for i in range(len(cafehs.study_ids))}
        else:
            cafehs.abs_neglogp_resid = {cafehs.study_ids[j]:np.array(pd.concat([pd.DataFrame(cafehs.abs_neglogp_resid[cafehs.study_ids[j]].T), pd.Series(abs_neglogp_residual[j])],axis=1)).T for j in range(len(cafehs.study_ids))}
            cafehs.precede_neglogp_resid = {cafehs.study_ids[j]:np.array(pd.concat([pd.DataFrame(cafehs.precede_neglogp_resid[cafehs.study_ids[j]].T), pd.Series(precede_neglogp_residual[j])],axis=1)).T for j in range(len(cafehs.study_ids))}
    
    #Save pickle file
    cafehs.save(f'{out_dir}/{file_prefix}.pkl', save_data=True)
    
###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
