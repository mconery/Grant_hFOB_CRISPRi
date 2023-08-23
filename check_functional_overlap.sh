#!/bin/bash

#SBATCH -o check_functional_overlap.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G

################################### Set Locations ####################################

finemap_file=/mnt/isilon/sfgi/conerym/data/fine-mapping_results/kanai_ukbb/release1.1/UKBB_94traits_release1.bed.gz
chip_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/ChIPSeq/hFOBs/chip/8591fb64-3c72-489c-a358-c65767083c25/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz
atac_peaks=/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/bone_diff/hFOBsDiff/out/peak/macs2/overlap/optimal_set/hFOB_diff_ppr.naive_overlap.filt.narrowPeak.gz
work_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection

######################################################################################

#Filter the fine-mapped results for the eBMD results for high-quality variants not in LD with SVs
zcat $finemap_file | awk '$12=="eBMD" && $19!=-1 && $22=="FALSE" && $23=="FALSE"' | sort -k1,1 -k2,2n > $work_dir/kanai_fine-mapped_eBMD_variants.bed

#Sort the peak files and add tags so we know if it's a chip or atac intersection later on
zcat $atac_peaks | sort -k1,1 -k2,2n | awk '{print $0"\tatac"}' > $work_dir/hFOB.atac_peaks.bed
zcat $chip_peaks | sort -k1,1 -k2,2n | awk '{print $0"\tchip"}' > $work_dir/hFOB.chip_peaks.bed

#Intersect the fine-mapped variants with the peak files
bedtools intersect -a $work_dir/kanai_fine-mapped_eBMD_variants.bed -b $work_dir/hFOB.atac_peaks.bed $work_dir/hFOB.chip_peaks.bed -wo -sorted > $work_dir/kanai_fine-mapped_eBMD_variants.atac_chip_peaks.bed