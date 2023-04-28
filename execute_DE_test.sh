#!/bin/bash

#SBATCH -o execute_DE_tests.log
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task 40
#SBATCH --mem=180G

################################### Set File Locations ##################################

de_test_script=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/Grant_hFOB_CRISPRi/test_differential_expression.R

######################################################################################

#Activate modules
module load R/4.2.2
module load hdf5/1.10.1
module load Pandoc/2.10

Rscript $de_test_script

