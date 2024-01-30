#!/bin/bash

#SBATCH -o plot_all_extracted_signals.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
plot_script=$script_dir/Multi-Trait_Fine-Mapping/4_analyze_BMD_signals/plot_BMD_signals.py

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
cafeh_dir=$finemap_dir/cafeh_results
plot_dir=$finemap_dir/bmd_signal_plots
residual_filtered_plot_dir=$finemap_dir/residual-filtered_signal_plots
temp_script_dir=$finemap_dir/temp_script_files

#Set global variables
min_purity=0.5
p_thresh=5e-8
activity_thresh=0.95

###########################################################################################################################################################

#Switch folder for temp file generation
cd $temp_script_dir

#Loop over the loci and plot each
for file in $(ls $cafeh_dir); do 
	file_prefix=${file/%.pkl/}
	#Make signal plot
	echo "source activate cafeh; python $plot_script -p $cafeh_dir/$file -o $plot_dir -u $min_purity -c $activity_thresh -m $p_thresh" > $temp_script_dir/$file_prefix.plot.sh
	sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.plot.sh
	sbatch --mem=4G -t 2:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.plot.sh
	#Make residual-filtered signal plot
	echo "source activate cafeh; python $plot_script -p $cafeh_dir/$file -o $plot_dir -u $min_purity -c $activity_thresh -m $p_thresh -r True" > $temp_script_dir/$file_prefix.residual-filtered.plot.sh
	sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.residual-filtered.plot.sh
	sbatch --mem=4G -t 2:00:00 --job-name $file_prefix.filtered $temp_script_dir/$file_prefix.residual-filtered.plot.sh
done



