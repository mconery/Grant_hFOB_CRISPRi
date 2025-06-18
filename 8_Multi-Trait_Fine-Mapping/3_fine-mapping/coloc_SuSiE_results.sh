#!/bin/bash

#SBATCH -o coloc_SuSiE_results.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=1G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
coloc_script=$script_dir/8_Multi-Trait_Fine-Mapping/3_fine-mapping/coloc_BMD_pairwise.R

#Define other file paths
finemap_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related
locus_dir=$finemap_dir/loci_files
locus_file=$locus_dir/bmd.sig_loci.csv
rds_dir=$finemap_dir/susie_results
bmd_trait="BMD"
out_dir=$finemap_dir/susie_coloc_results
temp_script_dir=$finemap_dir/temp_script_dir

###########################################################################################################################################################

#Change into temp script directory so all slurm outputs accumulate there
cd $temp_script_dir

#Loop over the loci and map each that needs it
awk -F "," '{print "chr"$1"."$2"."$3}' $locus_file | while read file_prefix; do 
	#Check if output file exists
	if [ ! -f $out_dir/$file_prefix".susie-coloc.tsv" ]; then
		#Call the Rscript for each locus
		echo "Rscript $coloc_script $file_prefix $rds_dir BMD $out_dir" > $temp_script_dir/$file_prefix.coloc.sh
      		sed -i '1i#!/bin/bash' $temp_script_dir/$file_prefix.coloc.sh
      		sbatch --mem=8G -t 1:00:00 --job-name $file_prefix $temp_script_dir/$file_prefix.coloc.sh
	fi
done



