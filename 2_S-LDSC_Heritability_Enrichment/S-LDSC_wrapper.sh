#!/bin/bash

cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC
files=($(ls $cdir/BEDfiles/bone/roadmap_epigenomics/*.bed $cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak/*.bed))
scriptdir=script_roadmap_epigenomics
logdir=log_roadmap_epigenomics
annodir=annotationfile_roadmap_epigenomics
outdir=out_roadmap_epigenomics

# Run partitioned LD score regression #########################
gwas_sumstats=($(echo $cdir/snps/bone_related_traits/non_merge/*.sumstats.gz ))
for gwas_sumstat in ${gwas_sumstats[@]}; do
    gwas_base=($(basename $gwas_sumstat ".sumstats.gz"))
       for file in ${files[@]};do
        base=$(basename $file ".bed" )
        dir="$base""_ldscore"
        basesrc="pldsr_""$base"
        echo ' 
		source ~/.bashrc  
		conda activate ldsc
		python ~/tools/ldsc/ldsc.py --h2 '$gwas_sumstat' --ref-ld-chr '$cdir/$annodir/$dir/$base'.,'$baseline'.  --w-ld-chr '$weight'.  --overlap-annot  --frqfile-chr '$frq'.  --out '$cdir/$outdir/$base'_'$gwas_base'' > $cdir/$scriptdir/$basesrc.$gwas_base.sh
        sed -i '1i#!/bin/bash'  $cdir/$scriptdir/$basesrc.$gwas_base.sh
        sbatch -c 7  --mem 8G -J pldsr.$base.$gwas_base -o $cdir/$logdir/$basesrc.$gwas_base.out $cdir/$scriptdir/$basesrc.$gwas_base.sh
    done
done

# Summarize results
for gwas_sumstat in ${gwas_sumstats[@]}; do
	gwas_base=($(basename $gwas_sumstat ".sumstats.gz"))
	echo "Celltype	Prop._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p" > $cdir/$outdir/summary_$gwas_base.txt
	ls $cdir/$outdir/*_${gwas_base}.results | while read file ; do
        base=$(basename $file "_${gwas_base}.results" )
        res=$(grep "L2_0" $cdir/$outdir/${base}_${gwas_base}.results | cut -f2- )
        echo $base $res | sed "s/ /\t/g" >> $cdir/$outdir/summary_$gwas_base.txt
	done
done
