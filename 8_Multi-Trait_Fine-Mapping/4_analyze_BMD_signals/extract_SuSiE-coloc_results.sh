#!/bin/bash

#SBATCH -o extract_SuSiE-coloc_results.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G

################################################################# Set files and locations #################################################################

#Define directory for SuSiE-Coloc outputs
susie_coloc_dir=/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/susie_coloc_results
#Give location of combination script
combination_script="/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi/8_Multi-Trait_Fine-Mapping/4_analyze_BMD_signals/combine_susie-coloc_output_files.py"

###########################################################################################################################################################

#Process signal files first
python $combination_script --directory $susie_coloc_dir --naming_convention "*.susie-coloc.tsv" --output $susie_coloc_dir/master.susie-coloc.tsv

#Process PP4s files second
python $combination_script --directory $susie_coloc_dir --naming_convention "*.signed_pp4s.tsv" --output $susie_coloc_dir/master.signed_pp4s.tsv
