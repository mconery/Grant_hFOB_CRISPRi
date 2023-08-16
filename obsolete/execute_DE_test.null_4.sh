#!/bin/bash

#SBATCH -o execute_DE_test.null_4.log
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task 60
#SBATCH --mem=180G

################################### Set User Variables ###############################

de_test_script=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/Grant_hFOB_CRISPRi/test_differential_expression_null.R
test_type=ntc
cell_type=all

######################################################################################

#Activate modules
module load R/4.2.3

Rscript $de_test_script $test_type $cell_type

