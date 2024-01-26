#!/bin/bash

#SBATCH -o extract_all_results.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
extract_script=$script_dir/Multi-Trait_Fine-Mapping/4_analyze_BMD_signals/extract_BMD_signals.py

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
cafeh_dir=$finemap_dir/cafeh_results
out_dir=$finemap_dir/summary_files
temp_script_dir=$finemap_dir/temp_script_files

###########################################################################################################################################################

#Set parameters for testing
assoc_threshes=(1 1e-6 5e-8)
assoc_types=(gwas absolute_residual)
purity_threshes=(0.95)

#Change the directory to the temp directory
cd $temp_script_dir

#Loop over the loci and plot each
for assoc_type in ${assoc_types[@]}; do 
	for assoc_thresh in ${assoc_threshes[@]}; do
		for purity_thresh in ${purity_threshes[@]}; do 
			echo "source activate cafeh; python $extract_script -p $cafeh_dir -o $out_dir -a $assoc_type -u $purity_thresh -m $assoc_thresh" > $temp_script_dir/$assoc_type.$purity_thresh.$assoc_thresh.extract.sh
			sed -i '1i#!/bin/bash' $temp_script_dir/$assoc_type.$purity_thresh.$assoc_thresh.extract.sh
			sbatch --mem=20G -t 2:00:00 --job-name $assoc_type.$purity_thresh.$assoc_thresh.extract $temp_script_dir/$assoc_type.$purity_thresh.$assoc_thresh.extract.sh

		done
	done
done




