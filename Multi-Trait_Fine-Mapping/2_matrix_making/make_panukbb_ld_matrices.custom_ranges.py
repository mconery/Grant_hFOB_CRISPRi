'''
make_panukbb_ld_matrices.custom_ranges.py

The purpose of this script is to make LD matrices and variant index files that 
can later be filtered for specific variants and ultimately used in fine-mapping 
applications. The ranges will be specified in an input file and then the script
will check the output directory to verify that the files don't already exist 
before it goes about making them.  
'''

#import modules that will be needed
import pandas as pd
import sys #Taking command line parameters as inputs
import getopt
import time
import hail as hl
from ukbb_pan_ancestry import *
from os.path import exists
from shutil import rmtree
import numpy as np
import os


#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -p => <str> Population to be used for the analysis REQUIRED
    -o => <str> folder for output matrix and index files REQUIRED
    -w => <str> comma-separated file of loci ranges REQUIRED
    -r => <int> number of digits to which to round ld matrices (default is no rounding/15bp) OPTIONAL
""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "p:o:w:r:nh")
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
        ukbb_pop = options_dict['-p']
        output_folder = options_dict['-o']
        loci_file = options_dict['-w']

    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()

    rounding = options_dict.get('-r', 15)
    #Confirm that the file locs exist
    if exists(loci_file) == False:
        print("ERROR: Locis ranges file does not exist.")
        sys.exit(1)
    #Confirm that rounding is an integer
    try:
        rounding = int(rounding)
    except ValueError:
        print("ERROR: Number of digits to round to must be positive integers.")
        sys.exit(1)
    print("Acceptable Inputs Given")
    #Call driver function
    driver(ukbb_pop, output_folder, loci_file, rounding)

###############################################################################
################################  TEST LOCATIONS ##############################
###############################################################################

ukbb_pop="EUR"
loci_file="/Volumes/endor/LD_matrices/bmd.sig_loci.csv"
rounding=5
output_folder="/Volumes/endor/LD_matrices"

###############################################################################
###########################  Helpful Functions  ###############################
###############################################################################

#Make a function that reads in a file loc and rounds the file there to a certain number of digits
def round_ld(matrix_loc, round_loc, digits):
    #Read file into a numpy array
    ld_raw = np.loadtxt(matrix_loc, dtype=float, delimiter=",")
    #Round the values
    ld_round = np.round(ld_raw, decimals=digits)
    del ld_raw
    #Write to file
    np.savetxt(round_loc, ld_round, delimiter=',')
    
#Make a function similar to round_ld that also deletes the original file
def round_ld_and_del(matrix_loc, round_loc, digits):
    #Read file into a numpy array
    ld_raw = np.loadtxt(matrix_loc, dtype=float, delimiter=",")
    #Round the values
    ld_round = np.round(ld_raw, decimals=digits)
    del ld_raw
    #Write to file
    np.savetxt(round_loc, ld_round, delimiter=',')
    #Delete original file
    os.remove(matrix_loc)

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
def driver(ukbb_pop, output_folder, loci_file, rounding):
    
    #Check if output_folder ends in a "/" and add if not
    output_folder = add_slash(output_folder)
    
    #Read in loci ranges file
    loci_ranges = []
    with open(loci_file) as loci_read:
        for line in loci_read.readlines():
            loci_ranges.append(line.strip().split(","))    

    #Make a list to hold the loci we actually want to make matrices for
    num_tiles = len(loci_ranges)
    #Get starting time
    start_time = time.time()
    #Get matrix
    bm = hl.linalg.BlockMatrix.read(get_ld_matrix_path(pop=ukbb_pop))
    #Get file path to variant ids
    ht_idx = hl.read_table(get_ld_variant_index_path(pop=ukbb_pop))
    #Loop over list of loci and make matrices
    for i in range(len(loci_ranges)):
        #Define locus
        tile = loci_ranges[i]
        #Make prefixes
        locus_prefix=tile[0] + ':' + tile[1] + '-' + tile[2]
        unrounded_loc = output_folder + 'panukbb.' + ukbb_pop + '.chr' + str(tile[0]) + '.' + str(tile[1]) + '.' + str(tile[2])
        if rounding == 15:
            export_loc = output_folder + 'panukbb.' + ukbb_pop + '.chr' + str(tile[0]) + '.' + str(tile[1]) + '.' + str(tile[2]) + '.unrounded'
        else:
            export_loc = output_folder + 'panukbb.' + ukbb_pop + '.chr' + str(tile[0]) + '.' + str(tile[1]) + '.' + str(tile[2]) + '.' + str(rounding) + '-digit'
        #Confirm that map file has not already been created yet. If it has, then assume both the matrix and map are needed.
        if not exists(unrounded_loc + '.map'):
    
            # filter by interval
            interval = hl.parse_locus_interval(locus_prefix)
            ht_idx_cur_tile = ht_idx.filter(interval.contains(ht_idx.locus))
            
            #Filter the block matrix
            cur_idx = ht_idx_cur_tile.idx.collect()
            if len(cur_idx) == 0: #Verify that there are at least some variants in the current window
                end_time = time.time()
                duration = round(end_time - start_time,2)
                print(export_loc + ' lacked variants. Completed in ' + str(duration) +' seconds (' + str(i+1) + ' of ' + str(num_tiles) + ')')
                #Reset start time
                start_time = end_time
                continue
            cur_bm = bm.filter(cur_idx, cur_idx)
    
            #Store the variant indices as a pandas df
            var_table = ht_idx_cur_tile.to_pandas()
            int_table = var_table.drop('idx', axis=1)
            #Prep the table for export
            temp_table = pd.DataFrame(int_table['alleles'].to_list())
            temp2_table = int_table['locus'].astype(str).str.split(":", expand=True)
            temp2_table['zero'] = 0
            unique_ids = int_table['locus'].astype(str) + "_" + temp_table[0] + "_" + temp_table[1]
            export_table = pd.concat([temp2_table[0], unique_ids, temp2_table['zero'], temp2_table[1]], axis=1)
    
            #Write pandas df to file (We're using the same map file name regardless of rounding)
            export_table.to_csv(unrounded_loc + '.map', sep='\t', header=False, index=False)
            #Check for whether rounding occurs and then export the filtered matrix
            if rounding == 15:
                # Note: when you apply any operation on BlockMatrix,
                # you need to write it first to storage before export
                cur_bm.write(unrounded_loc, force_row_major=True)
                hl.linalg.BlockMatrix.export(unrounded_loc,export_loc + '.ld.gz', delimiter=',')
                #Delete the directory
                rmtree(unrounded_loc)
            else: 
                # Note: when you apply any operation on BlockMatrix,
                # you need to write it first to storage before export
                cur_bm.write(unrounded_loc, force_row_major=True)
                hl.linalg.BlockMatrix.export(unrounded_loc, unrounded_loc + '.ld.gz', delimiter=',')
                #Delete the directory
                rmtree(unrounded_loc)
                #Round the new matrix
                round_ld_and_del(unrounded_loc + '.ld.gz', export_loc + '.ld.gz', rounding)
        
        #Check to see if only the matrix exists now. If it doesn't, then the map file must have been created for a different rounded value
        elif not exists(export_loc + '.ld.gz'): 
            #Check for whether rounding occurs and then export the filtered matrix
            if rounding == 15:
                #### In this case we need to make the matrix from scratch ####
                #Check for whether the unrounded_loc directory exists (if it does, the job was interrupted the last time, and it will be necessary to remove the half-written directory)
                if exists(unrounded_loc):
                    rmtree(unrounded_loc)
                # filter by interval
                interval = hl.parse_locus_interval(locus_prefix)
                ht_idx_cur_tile = ht_idx.filter(interval.contains(ht_idx.locus))
                
                #Filter the block matrix
                cur_idx = ht_idx_cur_tile.idx.collect()
                if len(cur_idx) == 0: #Verify that there are at least some variants in the current window
                    end_time = time.time()
                    duration = round(end_time - start_time,2)
                    print(export_loc + ' lacked variants. Completed in ' + str(duration) +' seconds (' + str(i+1) + ' of ' + str(num_tiles) + ')')
                    continue
                cur_bm = bm.filter(cur_idx, cur_idx)
                
                # Note: when you apply any operation on BlockMatrix,
                # you need to write it first to storage before export
                cur_bm.write(unrounded_loc, force_row_major=True)
                hl.linalg.BlockMatrix.export(unrounded_loc,export_loc + '.ld.gz', delimiter=',')
                #Delete the directory
                rmtree(unrounded_loc)
            else:
                #### Here we can check to see if we can get away without recreating the matrix, which we can do if the unrounded matrix already exists ####
                if exists(unrounded_loc + '.unrounded.ld.gz'):
                    #Check for whether the unrounded_loc directory exists (if it does, the job was interrupted the last time, and the unrounded file is unusable)
                    if exists(unrounded_loc):
                        os.remove(unrounded_loc + '.unrounded.ld.gz')
                        hl.linalg.BlockMatrix.export(unrounded_loc, unrounded_loc + '.ld.gz', delimiter=',')
                        round_ld_and_del(unrounded_loc + '.ld.gz', export_loc + '.ld.gz', rounding)
                    else:
                        #Round the existing file
                        round_ld(unrounded_loc + '.unrounded.ld.gz', export_loc + '.ld.gz', rounding)
                else: #The unrounded matrix must then be made from scratch
                    #Check for whether the unrounded_loc directory exists (if it does, the job was interrupted the last time, and the unrounded directory must be deleted else it will create a file write conflict)
                    if exists(unrounded_loc):
                        rmtree(unrounded_loc)    
                    # filter by interval
                    interval = hl.parse_locus_interval(locus_prefix)
                    ht_idx_cur_tile = ht_idx.filter(interval.contains(ht_idx.locus))
                    
                    #Filter the block matrix
                    cur_idx = ht_idx_cur_tile.idx.collect()
                    if len(cur_idx) == 0: #Verify that there are at least some variants in the current window
                        end_time = time.time()
                        duration = round(end_time - start_time,2)
                        print(export_loc + ' lacked variants. Completed in ' + str(duration) +' seconds (' + str(i+1) + ' of ' + str(num_tiles) + ')')
                        continue
                    cur_bm = bm.filter(cur_idx, cur_idx)
                    # Note: when you apply any operation on BlockMatrix,
                    # you need to write it first to storage before export
                    cur_bm.write(unrounded_loc, force_row_major=True)
                    hl.linalg.BlockMatrix.export(unrounded_loc, unrounded_loc + '.ld.gz', delimiter=',')
                    #Delete the directory
                    rmtree(unrounded_loc)
                    #Round the new matrix
                    round_ld_and_del(unrounded_loc + '.ld.gz', export_loc + '.ld.gz', rounding)
        
        #Case where the map and ld file both existed
        else:
            #Check the rounding because the two rounding types have different handling
            #In the first case a 15-digit matrix is written directly from the directory
            if rounding == 15 and exists(unrounded_loc):
                #Remove the final unrounded ld file in this case
                os.remove(export_loc + '.ld.gz')
                #Repurpose the directory to write the matrix
                hl.linalg.BlockMatrix.export(unrounded_loc, export_loc + '.ld.gz', delimiter=',')
                #Delete the directory
                rmtree(unrounded_loc)
            #In the next case a rounded ld file is being written from the unrounded file
            elif exists(unrounded_loc + '.ld.gz'): 
                #remove the rounded ld file
                os.remove(export_loc + '.ld.gz')
                #Round the new matrix
                round_ld_and_del(unrounded_loc + '.ld.gz', export_loc + '.ld.gz', rounding)
            #In this final case, everything was done already.
            else:
                #Make custom message for when the map and ld files both already existed
                end_time = time.time()
                duration = round(end_time - start_time,2)
                print(export_loc + ' previously existing. Completed in ' + str(duration) +' seconds (' + str(i+1) + ' of ' + str(num_tiles) + ')')
                #Reset start time
                start_time = end_time
                continue
        
        #Get new end time
        end_time = time.time()
        duration = round(end_time - start_time,2)
        #Print update
        print(export_loc + ' completed in ' + str(duration) +' seconds (' + str(i+1) + ' of ' + str(num_tiles) + ')')
        #Reset start time
        start_time = end_time



###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])
