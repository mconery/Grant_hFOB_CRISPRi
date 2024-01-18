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
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:u:c:m:nh")
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
    print("Acceptable Inputs Given")
    #Call driver function
    driver(pickle_file, out_dir, purity, active_thresh, max_assoc)

###############################################################################
################################  TEST LOCATIONS ##############################
###############################################################################

pickle_file="C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/pickle_files/chr1.10000001.11750000.pkl"
out_dir="C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping"
purity=0.5
active_thresh=0.5
max_assoc=5e-8

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

def driver(pickle_file, out_dir, purity, active_thresh, max_assoc): 

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
    #Identify BMD signals passing signficance thresholds
    BMD_sig_pure_active_signals = [x for x in BMD_pure_active_signals if np.any(cafehs.neglogP[BMD_trait_id, cafehs.credible_sets[x]] >= -(np.log10(max_assoc)))]
    #Identify traits ids and trait names that are active for at least one of the signals
    plottable_trait_ids = [x for x in list(range(len(cafehs.study_ids))) if np.any(cafehs.active[x,BMD_sig_pure_active_signals] > active_thresh)]
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
        #Get active components for the given trait
        active_components = np.array([x for x in BMD_pure_active_signals if cafehs.active[trait_id,x] > active_thresh]) #Components active in any trait
        assoc_array = cafehs.neglogP[trait_id,]
        for k in np.arange(cafehs.dims['K'])[active_components]:
            #Get signal color
            signal_color = color_map(k % num_k_colors)
            #Remove credible set SNPs from remaining snps
            remain_snps = np.delete(remain_snps, np.argwhere(np.isin(remain_snps,cafehs.snp_ids[credible_sets[k]])))
            #Plot the cred set
            axs[plot_position].scatter(
                cafehs.bp[credible_sets[k]],
                assoc_array[credible_sets[k]],
                label='k{}'.format(k),
                zorder = num_k_colors+1-k,
                color = signal_color)
        axs[plot_position].scatter(
            cafehs.bp[np.isin(cafehs.snp_ids, remain_snps)],
            assoc_array[np.isin(cafehs.snp_ids, remain_snps)],
            color="#DEDEDE",
            zorder = 1)
        axs[plot_position].set_title(trait.title())
        if(plot_position == len(plottable_traits)-1): #Add x-axis title for bottom subplot
            axs[plot_position].set_xlabel('chr' + str(np.max(cafehs.chr)))
    #Add shared y-label to figure
    fig.text(0.04, 0.5, '-log(P-value)', va='center', rotation='vertical')
    #Add shared  and adjust spacing
    fig.suptitle(file_prefix, fontsize=40)
    fig.tight_layout()
    fig.subplots_adjust(top=0.968)
    #Print info that will only go once on the figure
    plt.savefig(f'{out_dir}{file_prefix}.purity-{purity}.{max_assoc}.pvalue.png', bbox_inches='tight')
        
                
###############################################################################

#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

