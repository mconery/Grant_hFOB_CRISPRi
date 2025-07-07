#!/bin/bash

#Set general file locations
cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

scriptdir=/mnt/isilon/sfgi/conerym/analyses/grant/ldsc/bone_cells/scripts
annodir=annotationfile_roadmap_epigenomics
outdir=/mnt/isilon/sfgi/conerym/analyses/grant/ldsc/bone_cells/joint_osteoblasts

#Get BMD and osteoblast file locations
bmd_sumstat=$cdir/sumstats/bone_related_traits/non_merge/BMD.sumstats.gz
osteoblast_files=($(ls $cdir/BEDfiles/bone/roadmap_epigenomics/*OSTEO*.bed $cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak/*blasts*.bed $cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak/hFOB_Diff*.bed $cdir/BEDfiles/bone/hFOB_CRISPRi_GRCh38_Peak/hFOBs_chip*.bed))

# Prep the LDSC input for all osteoblast annotations #########################
gwas_base=($(basename $bmd_sumstat ".sumstats.gz"))
ld_chr_inp=""
for file in ${osteoblast_files[@]};do
	base=$(basename $file ".bed" )
        dir="$base""_ldscore"
	if [ "$ld_chr_inp" == "" ]; then
		ld_chr_inp=$cdir/$annodir/$dir/$base.
        else
		ld_chr_inp=$ld_chr_inp,$cdir/$annodir/$dir/$base.
	fi
done
ld_chr_inp=$ld_chr_inp,$baseline.

#Create and launch the script
echo ' 
	source ~/.bashrc  
	conda activate ldsc
	python ~/ldsc/ldsc.py --h2 '$gwas_sumstat' --ref-ld-chr '$ld_chr_inp'  --w-ld-chr '$weight'.  --overlap-annot  --frqfile-chr '$frq'.  --out '$outdir'/joint_osteoblasts_BMD' > $scriptdir/joint_osteoblasts_BMD.sh
sed -i '1i#!/bin/bash'  $scriptdir/joint_osteoblasts_BMD.sh
sbatch -c 7  --mem 8G -J joint_osteoblasts_BMD -o $scriptdir/joint_osteoblasts_BMD.out $scriptdir/joint_osteoblasts_BMD.sh

# Summarize results
echo "Celltype	Prop._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p" > $outdir/joint_osteoblasts_BMD.summary.txt
i=2
for file in ${osteoblast_files[@]};do
	base=$(basename $file ".bed" )
	awk -v i=$i -v t=$base 'NR==i{print t"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $outdir/joint_osteoblasts_BMD.results >> $outdir/joint_osteoblasts_BMD.summary.txt
	((i++))
done