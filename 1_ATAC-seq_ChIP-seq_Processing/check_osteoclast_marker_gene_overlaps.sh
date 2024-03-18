#!/bin/bash

#SBATCH -o check_osteoclast_marker_gene_overlap.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G

################################### Set Locations ####################################

#Set raw peak files
osteoclast_day4_full_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/GRCh38/osteoclasts_grant_Day4/atac/d95ef553-dba2-4b40-a087-8208a8f7f53c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
osteoclast_bae_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/GRCh38/osteoclasts_bae_2022/atac/228cfbe7-0792-4939-ba5a-74a69714afcd/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
osteoclast_bae_H3K27ac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/ChIPSeq/GRCh38/osteoclasts/Bae_2022/chip/c4621afa-8b76-4787-ae30-69bac519d788/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz

#Set marker gene promoter file
marker_gene_prom=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi/1_ATAC-seq_Processing/osteoclast_marker_gene_promoters.bed

#Set out directory
out_dir=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/marker_gene_intersects

#Set files and types (ADJUST IN FUTURE)
raw_files=($osteoclast_day4_full_atac_peaks $osteoclast_bae_atac_peaks $osteoclast_bae_H3K27ac_peaks) 
file_types=(osteoclast_day4_full osteoclast_bae_atac osteoclast_bae_H3K27ac)

######################################################################################

#Step 0: Make out directory if need be
mkdir -p $out_dir

#Step 1: Run intersection with promoter file
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
	bedtools intersect -a $marker_gene_prom -b ${raw_files[$i]} -c > $out_dir/${file_types[$i]}.promoter_overlaps.bed
done
