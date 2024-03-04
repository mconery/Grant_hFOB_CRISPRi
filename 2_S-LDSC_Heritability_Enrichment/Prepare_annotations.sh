#!/bin/bash

cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

# Prepare BED files #########################
dirin=/mnt/isilon/sfgi/conerym
## for roadmap epigenomics ChIP-seq
for file in $(ls $dirin/data/ChIP-seq/roadmap_epigenomics/*narrowPeak.gz | while read file ; do
        name=$(basename $file .narrowPeak.gz)
        zcat $file | cut -f 1-3 | sort -u | sort-bed - | awk -v c="$name" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > tmp
         ~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/roadmap_epigenomics/$name.bed out
done

## for 12 in-house data sets
dir=$cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak
zcat $dirin/analyses/grant/atacSeq/GRCh38/hFOBs/hFOB_Diff/atac/dc2ad981-9479-4fa7-919a-4c35bdecd9e2/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="hFOB_Diff_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/hFOB_Diff_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/hFOBs/hFOB_Perm/atac/e969a4d6-bc3d-4f0d-a260-c5cb0ac7df77/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="hFOB_Perm_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/hFOB_Perm_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/monocytes/control/atac/5df22173-71e7-446a-9025-a8463db5c9e6/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="monocytes_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/monocytes_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/pediatric_hMSCs/3D_Diff/atac/020fa98d-7ec4-4539-b716-650b1dea66c7/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="pediatric_hMSCs_Osteoblasts_3D_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/pediatric_hMSCs_Osteoblasts_3D_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/pediatric_hMSCs/6D_Diff/atac/ab104833-4251-4bce-bd9b-5bf974f12ad5/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="pediatric_hMSCs_Osteoblasts_6D_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/pediatric_hMSCs_Osteoblasts_6D_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/pediatric_hMSCs/Control/atac/a687fa3c-3d48-46f5-9620-b7226bba99b1/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="pediatric_hMSCs_Control_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/pediatric_hMSCs_Control_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/chondrocytes/atac/662ef7b0-7f8e-4b16-b95a-afcb57488a73/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="chondrocytes_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/chondrocytes_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/hMSC-Osteoblasts/atac/d8b71fb5-a388-4490-b05e-563406e404c8/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="adult_hMSC_Osteoblasts_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/adult_hMSC_Osteoblasts_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/osteoclasts_bae_2022/atac/228cfbe7-0792-4939-ba5a-74a69714afcd/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="osteoclasts_bae_2022_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/osteoclasts_bae_2022_atac.bed
zcat $dirin/analyses/grant/atacSeq/GRCh38/osteoclasts_grant_Day4/atac/d95ef553-dba2-4b40-a087-8208a8f7f53c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="osteoclasts_Grant_D4_atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/osteoclasts_Grant_D4_atac.bed
zcat $dirin/analyses/grant/ChIPSeq/GRCh38/hFOBs/chip/d172c438-f3bb-40e0-af7c-103ad389f5ab/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="hFOBs_chip" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/hFOBs_chip.bed
zcat $dirin/analyses/grant/ChIPSeq/GRCh38/osteoclasts/Bae_2022/chip/c4621afa-8b76-4787-ae30-69bac519d788/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="osteoclasts_Bae_2022_chip" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/osteoclasts_Bae_2022_chip.bed


# Make annotation #########################
files=($(ls $cdir/BEDfiles/bone/roadmap_epigenomics/*.bed $cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak/*.bed))
scriptdir=script_roadmap_epigenomics
logdir=log_roadmap_epigenomics
annodir=annotationfile_roadmap_epigenomics
outdir=out_roadmap_epigenomics
mkdir -p $cdir/$scriptdir
mkdir -p $cdir/$logdir
mkdir -p $cdir/$outdir

## Generate annotation files for each BED files
for file in ${files[@]};do
	base=$(basename $file ".bed" )
	dir="$base""_ldscore"
	mkdir -p $cdir/$annodir/$dir
	
	echo 'source ~/.bashrc
	conda activate ldsc
	for j in {1..22}; do
   		python ~/tools/ldsc/make_annot.py  --bed-file '$file' --bimfile '$plink'.$j.bim  --annot-file '$cdir/$annodir/$dir/$base'.$j.annot.gz
    done' > $cdir/$scriptdir/make_anno_$base.sh
    sed -i '1i#!/bin/bash' $cdir/$scriptdir/make_anno_$base.sh
    sbatch -J make_anno_$base -o $cdir/$logdir/make_anno_$base  $cdir/$scriptdir/make_anno_$base.sh
done

## Calculate ld for each annotation
for file in ${files[@]}; do
	base=$(basename $file ".bed" )
    dir="$base""_ldscore"
	basescr="calcLD_""$base"
	echo '
	source ~/.bashrc
	conda activate ldsc
	for j in {1..22}; do 
		python ~/tools/ldsc/ldsc.py --l2 --bfile '$plink'.$j --ld-wind-cm 1 --annot '$cdir/$annodir/$dir/$base'.$j.annot.gz --thin-annot  --out '$cdir/$annodir/$dir/$base'.$j --print-snps '$cdir/'list.txt
 	done' >  $cdir/$scriptdir/$basescr.sh
 	sed -i '1i#!/bin/bash' $cdir/$scriptdir/$basescr.sh
 	sbatch -c 8 --mem-per-cpu 16G -t 12:00:00 -J calcLD_$base -o $cdir/$logdir/$basescr.out $cdir/$scriptdir/$basescr.sh
done 
