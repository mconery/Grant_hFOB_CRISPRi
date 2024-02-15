#!/bin/bash

cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC
gwas_sumstats=($(echo $cdir/snps/bone_related_traits/merge/*.sumstats.gz ))
n_trait=$(ls $cdir/snps/bone_related_traits/merge/*.sumstats.gz | wc -l )

# Run genetic correlation for each traits versus itself and others
for i in $(seq 1 $n_trait) ; do
    e=$((i-1))
    one=${gwas_sumstats[$i-1]}
    self=$(echo $one | awk -F "/" '{print $NF}' | sed "s/.merged.sumstats.gz//g")
    left=($(echo "${gwas_sumstats[@]:$e:$n_trait}" | sed "s/ /,/g"))
    echo '
    source ~/.bashrc
    conda activate ldsc
    python ~/tools/ldsc/ldsc.py --rg '$one','$left' --ref-ld-chr '$baseline'. --w-ld-chr '$weight'. --out '$cdir'/Rg/bone_related_traits/'$self'.merged.overAll ' > $cdir/script/$self.merged.overAll.sh
    sed -i '1i#!/bin/bash'  $cdir/script/$self.merged.overAll.sh
    cd $cdir/script
    sbatch -c 16 --mem-per-cpu 8G -J ldsc_rg.$self.merged.overAll  $cdir/script/$self.merged.overAll.sh
done

# Summarize the results after runs finished
for file in $(ls $cdir/Rg/bone_related_traits/*.log) ; do
        sed -n '/p1/,/^$/p' $file | sed '$d'| sed '1d' >> $cdir/Rg/bone_related_traits/gR_alltraits.txt
        grep " SNPs with valid alleles" $file >> $cdir/Rg/bone_related_traits/gR_overlapSNPS.txt
done
sed -i "s/\/mnt\/isilon\/sfgi\/trangk\/analyses\/grant\/ldsr\/snps\/bone_related_traits\/merge\///g" $cdir/Rg/bone_related_traits/gR_alltraits.txt
sed -i "s/^ //g" $cdir/Rg/bone_related_traits/gR_alltraits.txt
sed -i "s/.merged.sumstats.gz//g" $cdir/Rg/bone_related_traits/gR_alltraits.txt
sed -i "1ip1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" $cdir/Rg/bone_related_traits/gR_alltraits.txt
sed -i 's/[ ][ ]*/ /g' $cdir/Rg/bone_related_traits/gR_alltraits.txt