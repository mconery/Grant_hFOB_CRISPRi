#!/bin/bash

#SBATCH -o check_and_map_loci_SuSiE.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
mapping_script=$script_dir/8_Multi-Trait_Fine-Mapping/3_fine-mapping/finemap_locus_SuSiE.py

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
out_dir=$finemap_dir/susie_results
munge_dir=$finemap_dir/munged_summary_stats
ld_dir=$finemap_dir/ld_matrices
size_file=$finemap_dir/trait_sample_sizes_cc.tsv
temp_mapping_dir=$finemap_dir/temp_mapping_files
temp_script_dir=$finemap_dir/temp_script_files
locus_dir=$finemap_dir/loci_files
locus_file=$locus_dir/bmd.sig_loci.csv
trait_json=$locus_dir/traits_per_loci.json

#Define mapping parameters
purity_thresh=0.01
confidence=0.95
num_signals=10
susie_mem=32

###########################################################################################################################################################

#Change into temp script directory so all slurm outputs accumulate there
cd $temp_script_dir

#Loop over the loci and map each that needs it
awk -F "," '{print "chr"$1"."$2"."$3}' $locus_file | while read file_prefix; do 

	#Call the python script for each locus
	echo "python $mapping_script -p $file_prefix -o $out_dir -m $munge_dir -l $ld_dir -n $size_file -k $num_signals -u $purity_thresh -f $confidence -j $trait_json -t $temp_mapping_dir -i $temp_script_dir -c TRAIT.sumstats.gz -e $susie_mem" > $temp_script_dir/$file_prefix.susie.sh
      sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.susie.sh
      sbatch --mem=64G -t 7-00:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.susie.sh


done



