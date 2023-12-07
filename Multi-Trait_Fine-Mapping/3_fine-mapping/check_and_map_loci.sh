#!/bin/bash

#SBATCH -o check_and_map_loci.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
cafeh_script=$script_dir/Multi-Trait_Fine-Mapping/3_fine-mapping/finemap_CAFEH.py

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
cafeh_dir=$finemap_dir/cafeh_results
munge_dir=$finemap_dir/munged_summary_stats
ld_dir=$finemap_dir/ld_matrices
size_file=$finemap_dir/trait_sample_sizes.tsv
temp_script_dir=$finemap_dir/temp_script_files
locus_dir=$finemap_dir/loci_files
locus_file=$locus_dir/bmd.sig_loci.csv
trait_json=$locus_dir/traits_per_loci.json

###########################################################################################################################################################

#Change into temp script directory so all slurm outputs accumulate there
cd $temp_script_dir

#Loop over the loci and map each that needs it
awk -F "," '{print "chr"$1"."$2"."$3}' $locus_file | while read file_prefix; do 

#Check to see whether the files have been made yet
if [ ! -e $cafeh_dir/$file_prefix.pkl ]; then
	echo "source activate cafeh; python $cafeh_script -p $file_prefix -o $cafeh_dir -m $munge_dir -l $ld_dir -n $size_file -j $trait_json" > $temp_script_dir/$file_prefix.sh
      sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.sh
      sbatch --mem=300G -t 7-00:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.sh
     
fi

done



