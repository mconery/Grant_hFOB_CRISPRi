'''
plot_CAFEH.py

This script will plot the results of fine mapping a particular locus with 
CAFEH. It will take a pickle file output from CAFEH and make a p-value plot.
If the option is selected, it will also make absolute or preceding residual
association plots.

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
    -t => <txt> plot type: gwas, absolute residual, or preceding residual  (default is gwas) OPTIONAL
    -u => <num> purity threshold for plotting non-residual plots (default is 0.5) OPTIONAL
    -m => <num> Necessary max association p-value for plotting credible sets on non-residual plots (default is 1) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:t:u:m:nh")
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
    plot_type = options_dict.get('-t', 'gwas')
    purity = options_dict.get('-u', 0.5)
    max_assoc = options_dict.get('-m', 1)
    #Confirm that the file/folder locs exist
    if exists(out_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    if exists(pickle_file) == False:
        print("ERROR: Pickle file of CAFEH output doesn't exist.")
        sys.exit(1)
    #Check that the plot_type makes sense
    if type(plot_type) != str:
        print("ERROR: Plot type is invalid. Select gwas, absolute_residual, or preceding_residual.")
        sys.exit(1)
    else:
        if plot_type.lower()[0] == 'g':
            plot_type = 'gwas'
            print("NOTE: GWAS plot type selected.")
        elif plot_type.lower()[0] == 'a':
            plot_type = 'absolute_residual'
            print("NOTE: Absolute residual plot type selected.")
        elif plot_type.lower()[0] == 'p':
            plot_type = 'preceding_residual'
            print("NOTE: Preceding residual plot type selected.")
        else: 
            print("ERROR: Plot type is invalid. Select gwas, absolute_residual, or preceding_residual.")
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
    driver(pickle_file, out_dir, plot_type, purity, max_assoc)

###############################################################################
################################  TEST LOCATIONS ##############################
###############################################################################

pickle_file="C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/pickle_files/chr1.10000001.11750000.pkl"
out_dir="C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping"
plot_type="absolute-residual"
purity=0.1
max_assoc=1

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

def driver(pickle_file, out_dir, plot_type, purity, max_assoc): 

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
    
    #Check what plot_type we have an print the corresponding plots
    #Non-residual case
    if plot_type == 'gwas':
        plt.clf() #Clear plot
        #Establish plot for number of traits
        fig, axs = plt.subplots(len(cafehs.study_ids), figsize=(10, len(cafehs.study_ids) * subplot_height))
        #Iterate over trait to make the subplots
        for trait_id in range(len(cafehs.study_ids)):
            #Get trait
            trait=cafehs.study_ids[trait_id]
            #Set an array of SNPs that aren't in any credible sets for the trait
            remain_snps = cafehs.snp_ids
            #Get components active for the given trait
            active_components = cafehs.active[trait_id,] > 0.5 #Components active in any trait
            assoc_array = cafehs.neglogP[trait_id,]
            for k in np.arange(cafehs.dims['K'])[active_components]:
                #Get signal color
                signal_color = color_map(k % num_k_colors)
                #Check that the credible set passes the purity and association thresholds
                if np.max(assoc_array[credible_sets[k]]) >= -np.log10(max_assoc) and cafehs.realpure[k] >= purity: 
                    #Remove credible set SNPs from remaining snps
                    remain_snps = np.delete(remain_snps, np.argwhere(np.isin(remain_snps,cafehs.snp_ids[credible_sets[k]])))
                    #Plot the cred set
                    axs[trait_id].scatter(
                        cafehs.bp[credible_sets[k]],
                        assoc_array[credible_sets[k]],
                        label='k{}'.format(k),
                        zorder = num_k_colors+1-k,
                        color = signal_color)
            axs[trait_id].scatter(
                cafehs.bp[np.isin(cafehs.snp_ids, remain_snps)],
                assoc_array[np.isin(cafehs.snp_ids, remain_snps)],
                color="#DEDEDE",
                zorder = 1)
            axs[trait_id].set_title(trait.title())
            if(trait_id == len(cafehs.study_ids)-1): #Add x-axis title for bottom subplot
                axs[trait_id].set_xlabel('chr' + str(np.max(cafehs.chr)))
        #Add shared y-label to figure
        fig.text(0.04, 0.5, '-log(P-value)', va='center', rotation='vertical')
        #Add shared  and adjust spacing
        fig.suptitle(file_prefix, fontsize=40)
        fig.tight_layout()
        fig.subplots_adjust(top=0.968)
        #Print info that will only go once on the figure
        plt.savefig(f'{out_dir}{file_prefix}.purity-{purity}.{max_assoc}.pvalue.png', bbox_inches='tight')
        
    #Residual case
    else:
        #Get the correct set of residuals
        if plot_type == 'absolute_residual':
            assoc_array = cafehs.abs_neglogp_resid
        else:
            assoc_array = cafehs.precede_neglogp_resid
        plt.clf() #Clear plot
        #Set subplot height
        subplot_width = 8
        #Establish plot for number of traits and residuals
        fig, axs = plt.subplots(num_k_colors, len(cafehs.study_ids), figsize=(num_k_colors * subplot_width, len(cafehs.study_ids) * subplot_height))
        #Iterate over trait to make the subplots
        for trait_id in range(len(cafehs.study_ids)):
            #Get trait
            trait=cafehs.study_ids[trait_id]
            #Set an array of SNPs that aren't in any credible sets for the trait
            remain_snps = cafehs.snp_ids
            for k in range(num_k_colors):
                #Get trait_assoc_array
                trait_assoc_array=assoc_array[trait]
                #Get signal color
                signal_color = color_map(k % num_k_colors)
                #Calculate SNP sets and residuals for plot
                logp_test = trait_assoc_array[k,:]
                cred_set = cafehs.snp_ids[credible_sets[k]]
                rest_set = np.array([i for i in cafehs.snp_ids if i not in cred_set])
                logp_test_cred = logp_test[np.sort(np.where(np.isin(cafehs.snp_ids, cred_set)))]
                logp_test_rest = logp_test[np.sort(np.where(np.isin(cafehs.snp_ids, rest_set)))]
                axs[k, trait_id].scatter(
                    cafehs.bp[np.sort(np.where(np.isin(cafehs.snp_ids, cred_set)))],
                    logp_test_cred,
                    zorder = 2,
                    color = signal_color)
                axs[k, trait_id].scatter(
                    cafehs.bp[np.sort(np.where(np.isin(cafehs.snp_ids, rest_set)))],
                    logp_test_rest,
                    zorder = 1,
                    color = "#DEDEDE")
                if k == 0:
                    axs[k, trait_id].set_title(trait.title())
                elif k == num_k_colors: #Add x-axis title for bottom subplot
                    axs[k, trait_id].set_xlabel('chr' + str(np.max(cafehs.chr)))
        #Add shared y-label to figure
        fig.text(0.04, 0.5, '-log(residual association)', va='center', rotation='vertical')
        #Add shared title and adjust spacing
        fig.suptitle(file_prefix, fontsize=40)
        fig.tight_layout()
        fig.subplots_adjust(top=0.968)
        #Save the figure
        plt.savefig(f'{out_dir}{file_prefix}.{plot_type}.png', bbox_inches='tight')
                
###############################################################################

#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

