#!/bin/bash

#SBATCH -o check_osteoclast_marker_gene_overlap.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G

################################### Set Locations ####################################

#Set raw peak files
osteoclast_day0_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day0/downsample/atac/f9c4b95b-1a2d-4585-a87b-6b64bb186651/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day4_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day4/downsample/atac/8cfab6c7-da4a-4091-98d3-07f05b31de20/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day8_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day8/downsample/atac/7252c593-db46-477a-9c81-b66b333d992d/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day0_full_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day0/atac/d79f8990-f683-4d7f-8b29-31cd3fdd94f6/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
osteoclast_day4_full_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day4/atac/388a6af8-cbec-40eb-bf89-a72ee06d521c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz

#Set marker gene promoter file
marker_gene_prom=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/osteoclast_marker_gene_promoters.bed

#Set out directory
out_dir=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/marker_gene_intersects

#Set files and types (ADJUST IN FUTURE)
raw_files=($osteoclast_day0_down_atac_peaks $osteoclast_day4_down_atac_peaks $osteoclast_day8_down_atac_peaks $osteoclast_day0_full_atac_peaks $osteoclast_day4_full_atac_peaks) 
file_types=(osteoclast_day0_down osteoclast_day4_down osteoclast_day8_down osteoclast_day0_full osteoclast_day4_full)

######################################################################################

#Step 0: Make out directory if need be
mkdir -p $out_dir

#Step 1: Loop over peak files and make temporary merged peak files
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
	bedtools merge ${raw_files[$i]} > $out_dir/${file_types[$i]}.merged_peaks.temp.bed
done

#Step 2: Run intersection with promoter file
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
	bedtools intersect -a $marker_gene_prom -b $out_dir/${file_types[$i]}.merged_peaks.temp.bed -loj > $out_dir/${file_types[$i]}.promoter_overlaps.bed
done

#Step 3:Remove temp files
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
	rm $out_dir/${file_types[$i]}.merged_peaks.temp.bed
done