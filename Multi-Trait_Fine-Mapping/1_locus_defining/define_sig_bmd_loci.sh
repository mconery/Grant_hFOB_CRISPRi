#!/bin/bash

#SBATCH -o define_sig_bmd_loci.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G

################################################################# Set files and locations #################################################################

#Set scripts directory and script locations
script_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi
locus_defining_script=$script_dir/Multi-Trait_Fine-Mapping/1_locus_defining/define_sig_loci_for_bmd.py
trait_identifying_script=$script_dir/Multi-Trait_Fine-Mapping/1_locus_defining/define_sig_loci_for_bmd.py

#Define other file paths
chromo_lengths=$script_dir/Multi-Trait_Fine-Mapping/1_locus_defining/grch37_chr_bounds_mhc_only.bed
munge_dir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsr/snps/bone_Mitch/non_merge
raw_bmd_gwas_file="/mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz"
locus_dir="/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files"

#Define global variables
prime_p_thresh=5e-8
sec_p_thresh=1e-6
tile_size=250000

###########################################################################################################################################################

#Define the significant loci for bmd and identify the SNPs falling into each locus
python $locus_defining_script -g $raw_bmd_gwas_file -o $locus_dir -p $prime_p_thresh -s $sec_p_thresh -c $chromo_lengths -t $tile_size

#Identify the traits that should be mapped at each locus
python $trait_identifying_script -l $locus_dir -g $munge_dir -p $sec_p_thresh
