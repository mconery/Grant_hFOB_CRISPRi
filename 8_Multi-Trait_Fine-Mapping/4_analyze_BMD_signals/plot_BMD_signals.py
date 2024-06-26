'''
plot_BMD_signals.py

This script will plot the results of fine mapping a particular locus with 
CAFEH with a focus on the BMD signals that will be extracted under the matched
purity and significance thresholds.

***This needs to be run on Python 3.9 as there an incompatability with other ***
***versions of the .to_list() function  and the CAFEH package for python.    ***
'''

###############################################################################

#Load libraries
import getopt
import numpy as np
import pandas as pd
import sys
from os.path import exists
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from cafeh.model_queries import *
from cafeh.cafeh_summary import CAFEHSummary
import pickle

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <txt> Pickled CAFEH output REQUIRED
    -o => <txt> output folder location REQUIRED
    -u => <num> purity threshold for plotting non-residual plots (default is 0.5) OPTIONAL
    -c => <num> activity threshold for calling a trait active (default is 0.5) OPTIONAL
    -m => <num> Necessary min association p-value for plotting credible sets on non-residual plots (default is 1) OPTIONAL
    -r => <boo> Filter for high-quality signals only when using GWAS Association for significance filtering.
                High quality signals carry peak of own residual association (default is False) Optional
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:u:c:m:r:nh")
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
        pickle_file = options_dict['-p']
        out_dir = options_dict['-o']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    
    #Get optional inputs
    purity = options_dict.get('-u', 0.5)
    active_thresh = options_dict.get('-c', 0.5)
    max_assoc = options_dict.get('-m', 1)
    filt_high_flag = options_dict.get('-r', False)
    #Confirm that the file/folder locs exist
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(pickle_file) == False:
        print("ERROR: Pickle file of CAFEH output doesn't exist.")
        sys.exit(1)
    #Recast purity value 
    try:
        purity = float(purity)
        if purity >= 1 or purity < 0: 
            print("ERROR: Purity level should be in [0,1).")
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
    #Recast maximum association value 
    try:
        max_assoc = float(max_assoc)
        if max_assoc > 1 or purity <= 0: 
            print("ERROR: Maximum association level should be in (0,1].")
            sys.exit(1)
    except ValueError:
        print("ERROR: Maximum association level is not coercible to a float. It should be a p-value in (0,1].")
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
        pass
    print("Acceptable Inputs Given")
    #Call driver function
    driver(pickle_file, out_dir, purity, active_thresh, max_assoc, filt_high_flag)

###############################################################################
################################  TEST LOCATIONS ##############################
###############################################################################

pickle_file="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/cafeh_results/chr1.21500001.24250000.pkl"
out_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/residual-filtered_signal_plots"
purity=0.5
active_thresh=0.5
max_assoc=5e-8
filt_high_flag=True

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

def driver(pickle_file, out_dir, purity, active_thresh, max_assoc, filt_high_flag): 

    #Extract file prefix from pickle_file
    file_prefix = pickle_file.split("/")[-1]
    file_prefix = file_prefix[:file_prefix.rfind(".")]
    
    #Add slash to out_dir if needed
    out_dir = add_slash(out_dir)
    
    #Read in pickle file
    with open(pickle_file, 'rb') as f:
        cafehs = pickle.load(f)
    
    #Check if the pickle file is empty
    if cafehs is None:
        print("WARNING: Attempted to plot unmapped locus that encountered RAM issue, " + file_prefix + '.')
        sys.exit(0)
    
    #Check for convergence
    if cafehs.check_convergence() == False:
        print("WARNING: Plotting results for a locus that failed to converge, " + file_prefix + '.')

    #Get credible sets
    credible_sets = cafehs.credible_sets

    # Use Agg backend for non-interactive plotting
    matplotlib.use('Agg')  
    #Set subplot height
    subplot_height = 8     
    # Define a custom color map with enough distinct colors for 'k'
    num_k_colors = cafehs.active.shape[1]  # Adjust this based on the expected number of 'k' values
    color_map = plt.cm.get_cmap('tab10', num_k_colors)
    
    ## Identify signals active in BMD that pass activity and purity thresholds ##
    #Get position of BMD in the list of traits
    BMD_trait_id = [x for x in list(range(len(cafehs.study_ids))) if list(cafehs.study_ids == 'BMD')[x]][0]
    BMD_active_signals = [x for x in list(cafehs.credible_sets.keys()) if (cafehs.active[BMD_trait_id,] > active_thresh).tolist()[x]]
    BMD_pure_active_signals = [x for x in BMD_active_signals if cafehs.realpure[x] > purity]
    
    #Make plot
    #Identify BMD signals passing signficance thresholds after checking whether the residual filtering flag has been thrown
    if filt_high_flag == False:
        BMD_sig_pure_active_signals = [x for x in BMD_pure_active_signals if np.any(cafehs.neglogP[BMD_trait_id, cafehs.credible_sets[x]] >= -(np.log10(max_assoc)))]
    else:
        BMD_sig_pure_active_signals = [x for x in BMD_pure_active_signals if np.any(cafehs.neglogP[BMD_trait_id, cafehs.credible_sets[x]] >= -(np.log10(max_assoc))) and np.argmax(cafehs.abs_neglogp_resid['BMD'][x,:]) in cafehs.credible_sets[x]]
    #Quit if there are no mappable signals
    if BMD_sig_pure_active_signals == []:
        print("NOTE: No mappable signals for " + file_prefix + '.')
        sys.exit(0)
    #Identify traits ids and trait names that are active and significant for at least one of the signals
    #Get the trait indices for the significant traits
    sig_trait_index_map = {}
    for signal in BMD_sig_pure_active_signals:
        temp = cafehs.neglogP[:,cafehs.credible_sets[signal]]
        if filt_high_flag == False:
            sig_trait_index_map[signal] = [list(cafehs.study_ids)[x] for x in range(len(cafehs.study_ids)) if np.any(temp[x,:] > -np.log10(max_assoc)) and cafehs.active[x,signal] > active_thresh]
        else:
            sig_trait_index_map[signal] = [list(cafehs.study_ids)[x] for x in range(len(cafehs.study_ids)) if np.any(temp[x,:] > -np.log10(max_assoc)) and cafehs.active[x,signal] > active_thresh and np.argmax(cafehs.abs_neglogp_resid[cafehs.study_ids[x]][signal,:]) in cafehs.credible_sets[signal]]
    plottable_trait_ids = [x for x in range(len(cafehs.study_ids)) if np.any([True if cafehs.study_ids[x] in z else False for z in sig_trait_index_map.values()])]
    plottable_traits = cafehs.study_ids[plottable_trait_ids]
    plottable_traits_sorted = np.sort(plottable_traits).tolist()
    plottable_traits_sorted.remove('BMD') #Removing to reorder for plot ordering purposes
    temp = ['BMD']
    temp.extend(plottable_traits_sorted)
    plottable_traits_sorted = temp[:]
    del temp
    #Make map of plottable trait ids to their position in the stacked figure putting BMD first and alphabetical after that
    position_map = {x : plottable_traits_sorted.index(x) for x in plottable_traits}            
    plt.clf() #Clear plot
    #Establish plot for number of traits
    fig, axs = plt.subplots(len(plottable_trait_ids), figsize=(10, len(plottable_trait_ids) * subplot_height))
    #Iterate over trait to make the subplots
    for trait_id in plottable_trait_ids:
        #Get trait
        trait=cafehs.study_ids[trait_id]
        #Get plotting position
        plot_position = position_map[trait]
        #Set an array of SNPs that aren't in any credible sets for the trait
        remain_snps = cafehs.snp_ids
        #Get association array for the trait
        assoc_array = cafehs.neglogP[trait_id,]
        for k in BMD_sig_pure_active_signals:
            #Skip signal if not relevant to the trait
            if trait not in sig_trait_index_map[k]:
                continue
            #Get signal color
            signal_color = color_map(k % num_k_colors)
            #Remove credible set SNPs from remaining snps
            remain_snps = np.delete(remain_snps, np.argwhere(np.isin(remain_snps,cafehs.snp_ids[credible_sets[k]])))
            #Plot the cred set after checking how to handle the index
            if len(plottable_traits) > 1:
                axs[plot_position].scatter(
                    cafehs.bp[credible_sets[k]],
                    assoc_array[credible_sets[k]],
                    label='k{}'.format(k),
                    zorder = num_k_colors+1-k,
                    color = signal_color)
                axs[plot_position].set_ylabel('')
            else:
                axs.scatter(
                    cafehs.bp[credible_sets[k]],
                    assoc_array[credible_sets[k]],
                    label='k{}'.format(k),
                    zorder = num_k_colors+1-k,
                    color = signal_color)
                axs.set_ylabel('')
        if len(plottable_traits) > 1:
            axs[plot_position].scatter(
                cafehs.bp[np.isin(cafehs.snp_ids, remain_snps)],
                assoc_array[np.isin(cafehs.snp_ids, remain_snps)],
                color="#DEDEDE",
                zorder = 1)
            axs[plot_position].set_title(trait.replace("_"," ").replace("whole body","whole-body").title())
            if(plot_position == len(plottable_traits)-1): #Add x-axis title for bottom subplot
                axs[plot_position].set_xlabel('chr' + str(np.max(cafehs.chr)))
        else:
            axs.scatter(
                cafehs.bp[np.isin(cafehs.snp_ids, remain_snps)],
                assoc_array[np.isin(cafehs.snp_ids, remain_snps)],
                color="#DEDEDE",
                zorder = 1)
            axs.set_title(trait.replace("_"," ").replace("whole body","whole-body").title())
            if(plot_position == len(plottable_traits)-1): #Add x-axis title for bottom subplot
                axs.set_xlabel('chr' + str(np.max(cafehs.chr)))
    #Add shared y-label to figure
    fig.text(0, 0.5, '-log(P-value)', va='center', rotation='vertical')
    fig.tight_layout()
    #Print info that will only go once on the figure
    if filt_high_flag == False:
        plt.savefig(f'{out_dir}{file_prefix}.purity-{purity}.{max_assoc}.pvalue.png', bbox_inches='tight')
    else:
        plt.savefig(f'{out_dir}{file_prefix}.purity-{purity}.{max_assoc}.residual-filtered.pvalue.png', bbox_inches='tight')
                
###############################################################################

#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

