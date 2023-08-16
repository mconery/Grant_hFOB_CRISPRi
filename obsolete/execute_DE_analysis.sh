#!/bin/bash

#SBATCH -o execute_DE_analysis.log
#SBATCH --time=24:00:00
#SBATCH --mem=64G

################################### Set File Locations ##################################

de_analysis_script=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/Grant_hFOB_CRISPRi/analyze_differential_expression.R

######################################################################################

#Activate modules
module load R/4.2.3

Rscript $de_analysis_script

