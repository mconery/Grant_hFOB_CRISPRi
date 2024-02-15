#!/bin/bash

cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

# Prepare BED files #########################
dirin=/mnt/isilon/sfgi/conerym
for file in $(ls $dirin/data/ChIP-seq/roadmap_epigenomics/*narrowPeak.gz | while read file ; do
        name=$(basename $file .narrowPeak.gz)
        zcat $file | cut -f 1-3 | sort -u | sort-bed - | awk -v c="$name" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > tmp
         ~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/roadmap_epigenomics/$name.bed out
done

zcat $dirin/analyses/grant/ChIPSeq/hFOBs/chip/8591fb64-3c72-489c-a358-c65767083c25/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="hFOB.H3K27ac.Cottone" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/BEDfiles/bone/tmp
~/tools/liftOver $cdir/BEDfiles/bone/tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/hFOB.H3K27ac.Cottone.bed out

zcat $dirin/analyses/grant/atacSeq/osteoclasts/Bae_2022/atac/c77ca17a-0943-459b-901f-e14e8c3f022c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="OSC.RANKL.Diff.Atac.Bae" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/BEDfiles/bone/tmp
~/tools/liftOver $cdir/BEDfiles/bone/tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/OSC.RANKL.Diff.Atac.Bae.bed out

zcat $dirin/analyses/grant/ChIPSeq/osteoclasts/Bae_2022/chip/e9ccfec0-dc8d-40cd-9c52-8c03e216dc47/call-reproducibility_overlap/execution/overlap.optimal_peak.narrowPeak.gz | cut -f 1-3 | sort -u | sort-bed - | awk -v c="OSC.RANKL.Diff.H3K27ac.Bae" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/BEDfiles/bone/tmp
~/tools/liftOver $cdir/BEDfiles/bone/tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/OSC.RANKL.Diff.H3K27ac.Bae.bed out

zcat $dirin/analyses/grant/atacSeq/osteoclasts/Day4/atac/388a6af8-cbec-40eb-bf89-a72ee06d521c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz | cut -f1-3 | sort -u | sort-bed - | awk -v c="OSC.Day4.Diff.atac" 'BEGIN{OFS="\t"}{name=c"_"NR; print $0,name}' > $dir/BEDfiles/bone/tmp
~/tools/liftOver $cdir/BEDfiles/bone/tmp ~/tools/hg19ToHg38.over.chain.gz $cdir/BEDfiles/bone/OSC.Day4.Diff.Atac.bed out

cut -f1-4 /mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/TOPMED_hg38/chondrocytes/merge_techreps/peaks/macs2/pool_reps/pooled_peaks_reproducible_filtered.bed > $cdir/BEDfiles/bone/chondrocytes.bed


# Make annotation #########################
files=($(ls $cdir/BEDfiles/bone/roadmap_epigenomics/*.bed $cdir/BEDfiles/bone/*.bed))
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
 	sbatch -c 8 --mem-per-cpu 16G -t 18:00:00 -J calcLD_$base -o $cdir/$logdir/$basescr.out $cdir/$scriptdir/$basescr.sh
done 