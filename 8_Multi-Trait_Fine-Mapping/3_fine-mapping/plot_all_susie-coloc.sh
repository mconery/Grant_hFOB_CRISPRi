#!/bin/bash

#SBATCH -o plot_all_susie-coloc.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
plot_script=$script_dir/8_Multi-Trait_Fine-Mapping/3_fine-mapping/plot_susie_coloc_locus.R

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
rds_dir=$finemap_dir/susie_results
susie_coloc_dir=$finemap_dir/susie_coloc_results
plot_dir=$finemap_dir/susie-coloc_plots
temp_script_dir=$finemap_dir/temp_script_files

###########################################################################################################################################################

#Switch folder for temp file generation
cd $temp_script_dir

#Loop over the loci and plot each
for file in $(ls $susie_coloc_dir/chr*.variants.bed); do 
temp=$(basename $file)
file_prefix=${temp%.variants.bed}
echo "Rscript $plot_script --locus_prefix $file_prefix --rds_dir $rds_dir --out_dir $plot_dir" > $temp_script_dir/$file_prefix.susie-coloc_plot.sh
sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.susie-coloc_plot.sh
sbatch --mem=8G -t 2:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.susie-coloc_plot.sh
done



