#!/bin/bash

#SBATCH -o calculate_hFOB_heritability_enrichment.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G

################################### Set Locations ####################################

diff_hfob_chip_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/ChIPSeq/hFOBs/chip/8591fb64-3c72-489c-a358-c65767083c25/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz
undiff_hfob_atac_peaks=/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/bone_diff/hFOBs/out/peak/macs2/overlap/optimal_set/hFOBs_ppr.naive_overlap.filt.narrowPeak.gz
diff_hfob_atac_peaks=/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/bone_diff/hFOBsDiff/out/peak/macs2/overlap/optimal_set/hFOB_diff_ppr.naive_overlap.filt.narrowPeak.gz
hmsc_bmp2_atac_peaks=/mnt/isilon/sfgi/chesia/analyses/grant/atacSeq/bone_cells/BMP2/ENCODE/all_reps/out/peak/macs2/overlap/optimal_set/BMP2_ppr.naive_overlap.filt.narrowPeak.gz
osteoclast_day0_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day0/downsample/atac/f9c4b95b-1a2d-4585-a87b-6b64bb186651/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day4_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day4/downsample/atac/8cfab6c7-da4a-4091-98d3-07f05b31de20/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day8_down_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day8/downsample/atac/7252c593-db46-477a-9c81-b66b333d992d/call-reproducibility_idr/execution/idr.conservative_peak.narrowPeak.gz
osteoclast_day0_full_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day0/atac/d79f8990-f683-4d7f-8b29-31cd3fdd94f6/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
osteoclast_day4_full_atac_peaks=/mnt/isilon/sfgi/conerym/analyses/grant/atacSeq/osteoclasts/Day4/atac/388a6af8-cbec-40eb-bf89-a72ee06d521c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz
bmd_gwas=/mnt/isilon/sfgi/conerym/data/gwas_summary_stats/pan_ukbb_raw/continuous-3148-both_sexes-irnt.appended.tsv
bmd_gwas_n=242809
ldsc_dir=/mnt/isilon/sfgi/conerym/analyses/grant/ldsc/bone_cells

#Set files and types (ADJUST IN FUTURE)
raw_files=($diff_hfob_chip_peaks $undiff_hfob_atac_peaks $diff_hfob_atac_peaks $hmsc_bmp2_atac_peaks $osteoclast_day0_down_atac_peaks $osteoclast_day4_down_atac_peaks $osteoclast_day8_down_atac_peaks $osteoclast_day0_full_atac_peaks $osteoclast_day4_full_atac_peaks) 
file_types=(diff_hfob_chip undiff_hfob_atac diff_hfob_atac hmsc_bmp2_atac osteoclast_day0_down osteoclast_day4_down osteoclast_day8_down osteoclast_day0_full osteoclast_day4_full)

######################################################################################

#Set file subdirectory names
bim_dir=$ldsc_dir/bim_files
sorted_peaks_dir=$ldsc_dir/sorted_peaks
merged_peaks_dir=$ldsc_dir/merged_peaks
scripts_dir=$ldsc_dir/scripts
directories=($bim_dir $sorted_peaks_dir $merged_peaks_dir $scripts_dir)
#Create directories if need be
for directory in ${directories[@]}; do
  mkdir -p $directory
done

#LD score Step 1: Sort Peak files
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
  zcat ${raw_files[$i]} | sort -k1,1 -k2,2n | gzip -c > $sorted_peaks_dir/"${file_types[$i]}.sorted.bed.gz"
done

#LD score Step 2: Merge Peak files
merged_files=()
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
  bedtools merge -i $sorted_peaks_dir/"${file_types[$i]}.sorted.bed.gz" | gzip -c > $merged_peaks_dir/"${file_types[$i]}.merged.bed.gz"
merged_files[$i]=$merged_peaks_dir/"${file_types[$i]}.merged.bed.gz"
done

#BONUS Step : Intersect Differentiated hFOB files (ADJUST IN FUTURE)
diff_hfob_intersect_peaks=$merged_peaks_dir/diff_hfob_intersect.merged.bed.gz
bedtools intersect -a $merged_peaks_dir/diff_hfob_chip.merged.bed.gz -b $merged_peaks_dir/diff_hfob_atac.merged.bed.gz | gzip -c > $diff_hfob_intersect_peaks
file_types[4]=diff_hfob_intersect
merged_files[4]=$merged_peaks_dir/"${file_types[4]}.merged.bed.gz"

#LD score Step 3: Make annotation files
cd $ldsc_dir
#Generate annotation files for each
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
  file=${merged_files[$i]}
  mkdir -p ${file_types[$i]}
echo 'source activate ldsc
  for j in {1..22}; do
     python /home/conerym/ldsc/make_annot.py \
      --bed-file '${merged_files[$i]}' \
      --bimfile /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/1000G_plinkfiles/1000G.mac5eur.$j.bim \
      --annot-file '$ldsc_dir/${file_types[$i]}/${file_types[$i]}'.$j.annot.gz
      done' > scripts/${file_types[$i]}.sh
      sed -i '1i#!/bin/bash' scripts/${file_types[$i]}.sh
      sbatch --mem=4G -t 6:00:00 --job-name make_anno scripts/${file_types[$i]}.sh
done

#LD score Step 4: Munge BMD Summary Stats
#Prep a munge friendly summary stats file (The Finucane provided one doesn't work with the Pan-UKBB stats)
cd $ldsc_dir
awk -v y=$bmd_gwas_n '{ if (NR == 1) { print "SNP\tA1\tA2\tbeta\tP\tN" } else if ($27 + 0 > 0 && $33 + 0 > 0) { print $1"\t"$5"\t"$6"\t"$27"\t"10^(-1*$39)"\t"y } else { print $1 "\t\t\t\t\t" } }' $bmd_gwas | awk '{ if ($5==0) {print $1"\t"$2"\t"$3"\t"$4"\t1.0e-308\t"$6} else {print $0} }' | gzip -c > $ldsc_dir/bmd.raw.gwas.gz
#Run the munge sumstats script now
conda activate ldsc
python ~/ldsc/munge_sumstats.py --sumstats $ldsc_dir/bmd.raw.gwas.gz --merge-alleles /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/w_hm3.snplist --chunksize=500000 --out bmd
conda deactivate
#Clean up temporary raw file
rm $ldsc_dir/bmd.raw.gwas.gz
  
#LD score Step 5: Calculate ld for each Custom Annotation
cdir=$ldsc_dir
cd $cdir
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
 file=${merged_files[$i]}
 echo 'source activate ldsc
   for j in {1..22}; do 
     python /home/conerym/ldsc/ldsc.py \
       --l2 \
       --bfile /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/1000G_plinkfiles/1000G.mac5eur.$j \
       --ld-wind-cm 1 \
       --annot '$ldsc_dir/${file_types[$i]}/${file_types[$i]}'.$j.annot.gz \
       --thin-annot \
       --out '$ldsc_dir/${file_types[$i]}/${file_types[$i]}'.$j \
       --print-snps /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/hapmap3_snps_fromgtz/hm.$j.snp
   done' > scripts/${file_types[$i]}.calcLD.sh
 sed -i '1i#!/bin/bash' $ldsc_dir/scripts/${file_types[$i]}.calcLD.sh
 sbatch --mem=32G -t 48:00:00 --job-name calc_ld $ldsc_dir/scripts/${file_types[$i]}.calcLD.sh
done 

#LD score Step 6: Run partitioned LD score regression 
#(first line is the output, P value should be adusted for multiple hypothesis testing with other Pvalues you use
#(if comparing multiple celltypes or GWAS)
mkdir -p $ldsc_dir/out
cdir=$ldsc_dir
cd $cdir
gwas_sumstats=($(echo *.sumstats.gz))
for gwas_sumstat in ${gwas_sumstats[@]}; do
 gwas_base=($(basename $gwas_sumstat | cut -d . -f 1))
 cd $cdir
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
 file=${merged_files[$i]}
 echo 'source activate ldsc
 python /mnt/isilon/sfgi/programs/ldsc/ldsc.py \
  --h2 '$gwas_sumstat' \
  --ref-ld-chr '$ldsc_dir/${file_types[$i]}/${file_types[$i]}'.,/mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/baseline/baseline. \
  --w-ld-chr /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/1000G_frq/1000G.mac5eur. \
  --out '$ldsc_dir'/out/'${file_types[$i]}'_'$gwas_base'' > "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 sed -i '1i#!/bin/bash' "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 sbatch --mem=4G -t 2:00:00 --job-name run_pldscore_null "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 done
done

#LD score Step 7: (OPTIONAL) Run partitioned LD score regression without baseline model
#(first line is the output, P value should be adusted for multiple hypothesis testing with other Pvalues you use
#(if comparing multiple celltypes or GWAS)
mkdir -p $ldsc_dir/out_no_baseline
cdir=$ldsc_dir
cd $cdir
gwas_sumstats=($(echo *.sumstats.gz))
for gwas_sumstat in ${gwas_sumstats[@]}; do
 gwas_base=($(basename $gwas_sumstat | cut -d . -f 1))
 cd $cdir
for ((i = 0 ; i < ${#file_types[@]} ; i++)); do
 file=${merged_files[$i]}
 echo 'source activate ldsc
 python /mnt/isilon/sfgi/programs/ldsc/ldsc.py \
  --h2 '$gwas_sumstat' \
  --ref-ld-chr '$ldsc_dir/${file_types[$i]}/${file_types[$i]}'. \
  --w-ld-chr /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr /mnt/isilon/sfgi/pahlm/analyses/grant/disease/partitioned_ldsr_phase1/1000G_frq/1000G.mac5eur. \
  --out '$ldsc_dir'/out_no_baseline/'${file_types[$i]}'_'$gwas_base'' > "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 sed -i '1i#!/bin/bash' "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 sbatch --mem=4G -t 2:00:00 --job-name run_pldscore_null "scripts/run_pldsr_"${file_types[$i]}"_"$gwas_base.sh
 done
done

