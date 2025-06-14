'''
extract_BMD_signals.py

This script will extract all the relevant information from the CAFEH pickle 
files. It outputs 4 files as follows:
    1) File of generic info on each signal
    2) Bed file where each row is a SNP-signal combo
    3) Matrix of activity scores for each signal
    4) Matrix of weights for each signal

***This needs to be run on Python 3.9 as there an incompatability with other ***
***versions of the .to_list() function  and the CAFEH package for python.    ***
***This is also incompatible with the petagene module.                       ***
'''

###############################################################################
    
#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
import os
from os.path import exists
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import cafeh

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> Pickled CAFEH output directory REQUIRED
    -o => <txt> output folder location REQUIRED
    -t => <txt> file of trait sizes (default is /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes.tsv) OPTIONAL
    -a => <txt> association type for thresholding: gwas, absolute residual, or preceding residual  (default is gwas) OPTIONAL
    -c => <num> activity threshold for calling a trait active (default is 0.5) OPTIONAL
    -u => <num> purity threshold for thresholding (default is 0.5) OPTIONAL
    -m => <num> Necessary min association p-value for thresholding (default is 1) OPTIONAL
    -r => <boo> Filter for high-quality signals only when using GWAS Association for significance filtering.
                High quality signals carry peak of own residual association (default is False) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:t:a:c:u:m:r:nh")
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
        pickle_dir = options_dict['-p']
        out_dir = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Get optional inputs
    trait_sizes_file = options_dict.get('-t', '/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes.tsv')
    assoc_type = options_dict.get('-a', 'gwas')
    active_thresh = options_dict.get('-c', 0.5)
    purity = options_dict.get('-u', 0.5)
    min_assoc = options_dict.get('-m', 1)
    filt_high_flag = options_dict.get('-r', False)
    #Confirm that the file/folder locs exist
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(pickle_dir) == False:
        print("ERROR: Directory of CAFEH pickle files doesn't exist.")
        sys.exit(1)
    if exists(trait_sizes_file) == False:
        print("ERROR: File of trait sizes needed for standardizing matrix columns doesn't exist.")
        sys.exit(1)
    #Check that the assoc_type makes sense
    if type(assoc_type) != str:
        print("ERROR: Association type is invalid. Select gwas, absolute_residual, or preceding_residual.")
        sys.exit(1)
    else:
        if assoc_type.lower()[0] == 'g':
            assoc_type = 'gwas'
            print("NOTE: GWAS association type selected.")
        elif assoc_type.lower()[0] == 'a':
            assoc_type = 'absolute_residual'
            print("NOTE: Absolute residual association type selected.")
        elif assoc_type.lower()[0] == 'p':
            assoc_type = 'preceding_residual'
            print("NOTE: Preceding residual association type selected.")
        else: 
            print("ERROR: Association type is invalid. Select gwas, absolute_residual, or preceding_residual.")
            sys.exit(1)
    #Recast purity value 
    try:
        purity = float(purity)
        if purity > 1 or purity < 0: 
            print("ERROR: Purity level should be in [0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Purity level is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Recast active value 
    try:
        active_thresh = float(active_thresh)
        if active_thresh > 1 or active_thresh <= 0: 
            print("ERROR: Activity threshold should be in (0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Activity threshold is not coercible to a float. It should be a proportion in [0,1).")
        sys.exit(1)
    #Recast minimum association value needed for retaining
    try:
        min_assoc = float(min_assoc)
        if min_assoc > 1 or min_assoc <= 0: 
            print("ERROR: Minimum association level should be in (0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Minimum association level is not coercible to a float. It should be a p-value in (0,1].")
        sys.exit(1)
    #Recast High-quality signal filtering flag
    if(type(filt_high_flag)) == str:
        if(filt_high_flag.lower() == 'false' or filt_high_flag.lower() == 'f' or filt_high_flag.lower() == 'no' or filt_high_flag.lower() == 'n'):
            filt_high_flag = False
        elif(filt_high_flag.lower() == 'true' or filt_high_flag.lower() == 't' or filt_high_flag.lower() == 'yes' or filt_high_flag.lower() == 'y'):
            filt_high_flag = True
        else:
            print("ERROR: Unrecognized option selected for whether to filter for signals where credible set captures peak of residual association. Select True/False.")
            sys.exit(1)
    elif(type(filt_high_flag) == bool):
        #In this case the indicator is already a boolean, so there's no need to recast.
        #Verify that there isn't a conflict betweent the association type and the high-quality signal filtering flag
        if filt_high_flag == True and assoc_type != 'gwas':
            print("ERROR: Filtering for signals that contain peak of own absolute residual is only compatible with initial filtering by gwas p-value.")
            sys.exit(1)
    else:
        print("ERROR: Unrecognized option selected for whether to filter for signals where credible set captures peak of residual association. Select True/False.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(pickle_dir, out_dir, trait_sizes_file, assoc_type, active_thresh, purity, min_assoc, filt_high_flag)

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

def driver(pickle_dir, out_dir, trait_sizes_file, assoc_type, active_thresh, purity, min_assoc, filt_high_flag): 

    #Add slashes to directories if needed    
    out_dir = add_slash(out_dir)
    pickle_dir = add_slash(pickle_dir)
    #Get list of traits from sample sizes file
    trait_sizes = pd.read_table(trait_sizes_file, sep="\t")
    traits = trait_sizes.loc[:,'trait'].to_list()
    traits = [trait for trait in traits if trait != 'BMD']
    #Get list of pickled outputs
    pickle_files = os.listdir(pickle_dir)
    pickle_files = [x for x in pickle_files if x[-4:] == '.pkl']
    #Extract loci prefixes from pickle_files
    loci = [x.replace(".pkl", "") for x in pickle_files]
    
    #Loop over the loci and extract all the BMD signals at the given thresholds
    bmd_signals_weights = [] #Empty list to hold the unbinarized weights
    bmd_signals_activity = [] #Empty list to hold the unbinarized activity scores
    bmd_signals_activity_binarized = [] #Empty list to hold the binarized activity scores
    bmd_signals_out = [] #Empty list to hold the writeable dataframe
    bmd_bed = [] #Empty list to hold the bed file formatted data
    for locus in loci:
        print('NOTE: Beginning locus ' + locus)
        #Read in pickle file
        with open(pickle_dir + locus + '.pkl', 'rb') as f:
            cafehs = pickle.load(f)
        #Check if the pickle file is empty
        if cafehs is None:
            print("WARNING: Skipping locus that encountered RAM issue, " + locus + '.')
            continue
        #Check for convergence
        if cafehs.check_convergence() == False:
            print("WARNING: Skipping results for a locus that failed to converge, " + locus + '.')
            continue
        #Identify index of BMD data
        bmd_index = [x for x in range(len(cafehs.study_ids)) if cafehs.study_ids[x] == 'BMD'][0]
        #Identify active BMD signals
        bmd_signals = [x for x in range(len(cafehs.active[bmd_index,:])) if cafehs.active[bmd_index,x] > active_thresh]
        #Identify pure, active BMD signals
        pure_bmd_signals = [x for x in bmd_signals if cafehs.realpure[x] > purity]
        #Check significance thresholds (first check the association type used for thresholding)
        if filt_high_flag == True:
            bmd_assoc = cafehs.neglogP[bmd_index,:]
            bmd_residual_assoc = cafehs.abs_neglogp_resid['BMD']
            sig_pure_bmd_signals = [x for x in pure_bmd_signals if np.sum(bmd_assoc[cafehs.credible_sets[x]] > -np.log10(min_assoc)) > 0 and np.argmax(bmd_residual_assoc[x,]) in cafehs.credible_sets[x]]
        elif assoc_type == 'gwas':
            bmd_assoc = cafehs.neglogP[bmd_index,:]
            sig_pure_bmd_signals = [x for x in pure_bmd_signals if np.sum(bmd_assoc[cafehs.credible_sets[x]] > -np.log10(min_assoc)) > 0]
        elif assoc_type == 'absolute_residual':
            assoc_array = cafehs.abs_neglogp_resid
            bmd_assoc = assoc_array['BMD']
            sig_pure_bmd_signals = [x for x in pure_bmd_signals if np.sum(bmd_assoc[x,cafehs.credible_sets[x]] > -np.log10(min_assoc)) > 0]
        else:
            assoc_array = cafehs.abs_neglogp_resid
            bmd_assoc = assoc_array['BMD']
            sig_pure_bmd_signals = [x for x in pure_bmd_signals if np.sum(bmd_assoc[x,cafehs.credible_sets[x]] > -np.log10(min_assoc)) > 0]
        if sig_pure_bmd_signals == []:
            continue
        #Get the trait indices
        trait_indices = [int(np.where(cafehs.study_ids == x)[0]) if x in cafehs.study_ids else -1 for x in traits]
        #Get the trait indices for the significant, active traits
        sig_trait_index_map = {}
        if filt_high_flag == True:
            for signal in sig_pure_bmd_signals:
                temp = cafehs.neglogP[:,cafehs.credible_sets[signal]]
                sig_trait_index_map[signal] = [trait_indices[x] if (trait_indices[x] == -1) else trait_indices[x] if np.any(temp[trait_indices[x],:] > -np.log10(min_assoc)) and cafehs.active[trait_indices[x],signal] > active_thresh and np.argmax(cafehs.abs_neglogp_resid[traits[x]][signal,:]) in cafehs.credible_sets[signal] else -1 for x in range(len(trait_indices)) ]
        elif assoc_type == 'gwas':
            for signal in sig_pure_bmd_signals:
                temp = cafehs.neglogP[:,cafehs.credible_sets[signal]]
                sig_trait_index_map[signal] = [x if (x == -1) else x if np.any(temp[x,:] > -np.log10(min_assoc)) and cafehs.active[x,signal] > active_thresh else -1 for x in trait_indices ]
        else:
            for signal in sig_pure_bmd_signals:
                temp = {x : assoc_array[x][signal,cafehs.credible_sets[signal]] for x in list(assoc_array.keys())}
                sig_trait_index_map[signal] = [x if (x == -1) else x if np.any(temp[cafehs.study_ids[x]] > -np.log10(min_assoc)) and cafehs.active[x,signal] > active_thresh else -1 for x in trait_indices]
        #Reset the weights using the hard-coded stored valued independent of LD matrices
        cafehs.weight_means = cafehs.realweight_means
        #Reorient the weights so that the BMD signal always has a positive weight for all signals
        bmd_normalized_weights = ([1 if x > 0 else -1 for x in cafehs.get_expected_weights()[bmd_index,:]] * cafehs.get_expected_weights())
        #Get the active values and weight values
        for signal in sig_pure_bmd_signals:
            bmd_signals_activity.append(np.vstack([cafehs.active[x,signal] if x != -1 else 0 for x in trait_indices]).T)
            bmd_signals_weights.append(np.vstack([bmd_normalized_weights[x,signal] if x != -1 else 0 for x in sig_trait_index_map[signal]]).T)
            bmd_signals_activity_binarized.append(np.vstack([1 if x != -1 else 0 for x in sig_trait_index_map[signal]]).T)
            #Extract the signal info for a data frame
            temp = []
            sig_trait_indices = [x for x in sig_trait_index_map[signal] if x != -1]
            if sig_trait_indices == []:
                temp = pd.Series([locus + '.' + str(signal), locus, signal, 
                                   len(cafehs.credible_sets[signal]), 
                                   ','.join(cafehs.ID[cafehs.credible_sets[signal]].tolist()), 
                                   ','.join(cafehs.rsid[cafehs.credible_sets[signal]].tolist()),
                                   cafehs.realpure[signal],
                                   cafehs.active[bmd_index,signal],
                                   np.max(cafehs.neglogP[bmd_index,cafehs.credible_sets[signal]]),
                                   str(cafehs.abs_neglogp_resid['BMD'][signal,:].argmax() in cafehs.credible_sets[signal]),
                                   np.max(cafehs.get_pip()[cafehs.credible_sets[signal]]), 
                                   0,
                                   'None', 
                                   'NA',
                                   'NA',
                                   'NA',
                                   'NA'
                                   ])
            else:
                temp = pd.Series([locus + '.' + str(signal), locus, signal, 
                                   len(cafehs.credible_sets[signal]), 
                                   ','.join(cafehs.ID[cafehs.credible_sets[signal]].tolist()), 
                                   ','.join(cafehs.rsid[cafehs.credible_sets[signal]].tolist()),
                                   cafehs.realpure[signal],
                                   cafehs.active[bmd_index,signal],
                                   np.max(cafehs.neglogP[bmd_index,cafehs.credible_sets[signal]]),
                                   str(cafehs.abs_neglogp_resid['BMD'][signal,:].argmax() in cafehs.credible_sets[signal]),
                                   np.max(cafehs.get_pip()[cafehs.credible_sets[signal]]), 
                                   len(sig_trait_indices),
                                   ','.join(cafehs.study_ids[sig_trait_indices].tolist()), 
                                   ','.join([str(x) for x in cafehs.active[sig_trait_indices,signal]]),
                                   ','.join([str(x) for x in np.max(cafehs.neglogP[sig_trait_indices,:][:,cafehs.credible_sets[signal]], axis=1).tolist()]),
                                   ','.join([str(u) for u in np.max(np.vstack([z[signal,cafehs.credible_sets[signal]] for z in [cafehs.abs_neglogp_resid[y] for y in cafehs.study_ids[sig_trait_indices].tolist()]]), axis = 1)]),
                                   ','.join([str(True) if z[signal,:].argmax() in cafehs.credible_sets[signal] else str(False) for z in [cafehs.abs_neglogp_resid[y] for y in cafehs.study_ids[sig_trait_indices].tolist()]])
                                   ])
            bmd_signals_out.append(temp)
            #Extract the bed file info
            temp = np.vstack([cafehs.chr[cafehs.credible_sets[signal]],
                              cafehs.bp[cafehs.credible_sets[signal]],
                              cafehs.bp[cafehs.credible_sets[signal]] + 1,
                              locus + '.' + str(signal) + '.' + cafehs.ID[cafehs.credible_sets[signal]], 
                              cafehs.rsid[cafehs.credible_sets[signal]],
                              cafehs.get_pip()[cafehs.credible_sets[signal]],
                              np.tile(','.join(list(cafehs.study_ids[sig_trait_indices])),cafehs.credible_sets[signal].shape),
                              np.array([','.join(list(cafehs.neglogP[sig_trait_indices,:][:,cafehs.credible_sets[signal]][:,y].astype(str))) for y in range(len(cafehs.credible_sets[signal]))])]).T
            bmd_bed.append(temp)
        
    #Make matrices 
    bmd_signals_matrix = pd.DataFrame(np.vstack(bmd_signals_out), 
                                      columns=["signal_id", "locus", "signal", "num_snps", "snp_ids", "rsids", "purity", "bmd_activity", "max_bmd_neglog_gwas","bmd_capture_max_residual", "max_pip", "num_traits", "traits", "trait_activities","max_other_neglog_gwas", "max_other_neglog_residual", "other_capture_max_residual"])
    bmd_bed_matrix = pd.DataFrame(np.vstack(bmd_bed),
                                  columns = ["chr","start","end","signal.snp_id","rsid","total_pip","active_traits","active_traits_gwas_neglogp"])
    bmd_activity_matrix = pd.DataFrame(np.vstack(bmd_signals_activity), columns=traits)
    bmd_activity_matrix = bmd_activity_matrix.set_index(bmd_signals_matrix['signal_id'])
    bmd_weights_matrix = pd.DataFrame(np.vstack(bmd_signals_weights), columns=traits)
    bmd_weights_matrix = bmd_weights_matrix.set_index(bmd_signals_matrix['signal_id'])
    bmd_activity_matrix_binarized = pd.DataFrame(np.vstack(bmd_signals_activity_binarized), columns=traits)
    bmd_activity_matrix_binarized = bmd_activity_matrix_binarized.set_index(bmd_signals_matrix['signal_id'])
    
    #Write matrices to file after checking for high-quality filtering for naming conventions
    out_file_suffix = '.purity_' + str(purity) + '.activity_' + str(active_thresh) + '.' + assoc_type + '_' + str(min_assoc) + '.'
    if filt_high_flag == True:
        out_file_suffix = out_file_suffix + 'highest_residual_filtered.' 
    bmd_signals_matrix.to_csv(out_dir + 'bmd_signals' + out_file_suffix + 'tsv', index=True, sep='\t')
    bmd_bed_matrix.to_csv(out_dir + 'bmd_variants' + out_file_suffix + 'bed', index=True, sep='\t')
    bmd_activity_matrix.to_csv(out_dir + 'bmd_signal_activity' + out_file_suffix + 'tsv', index=True, sep='\t')
    bmd_weights_matrix.to_csv(out_dir + 'bmd_signal_weights' + out_file_suffix + 'tsv', index=True, sep='\t')    
    bmd_activity_matrix_binarized.to_csv(out_dir + 'bmd_signal_activity' + out_file_suffix + 'binarized.tsv', index=True, sep='\t')
    
###############################################################################

#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

