#!/bin/bash

#SBATCH -o plot_all_results.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
plot_script=$script_dir/Multi-Trait_Fine-Mapping/3_fine-mapping/plot_CAFEH.py

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
cafeh_dir=$finemap_dir/cafeh_results
plot_dir=$finemap_dir/signal_plots
temp_script_dir=$finemap_dir/temp_script_files

###########################################################################################################################################################

#Switch folder for temp file generation
cd $temp_script_dir

#Loop over the loci and plot each
for file in $(ls $cafeh_dir); do 
file_prefix=${file/%.pkl/}
echo "source activate cafeh; python $plot_script -p $cafeh_dir/$file -o $plot_dir -t absolute_residual" > $temp_script_dir/$file_prefix.plot.sh
sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.plot.sh
sbatch --mem=4G -t 2:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.plot.sh
done



