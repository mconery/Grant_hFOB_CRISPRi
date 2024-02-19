#!/bin/bash

cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
cd $cdir/snps/bone_related_traits

mkdir TSV
ls /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/pan_ukbb_raw/*.tsv.gz | while read file ; do
    name=$(basename $file ".tsv.gz")
    if grep -q "neglog10_pval_EUR" <(zcat $file | head -n 1) ; then
        col=$(zcat $file | head -n 1 | sed "s/\s\+/\n/g" | grep -n "neglog10_pval_EUR" | cut -d':' -f1)
        zcat $file | sed '1d' | awk -v col="$col" '{p= 10^(-$col) ; print p}' > p
        sed -i "1ipval_EUR" p
        paste -d'\t' <(zcat $file) $cdir/snps/bone_related_traits/pan_ukbb_raw/panukbb_rsids.txt p > $cdir/snps/bone_related_traits/TSV/$name.tsv
    else
        paste -d'\t' <(zcat $file ) $cdir/snps/bone_related_traits/pan_ukbb_raw/panukbb_rsids.txt > $cdir/snps/bone_related_traits/TSV/$name.tsv
    fi
    gzip $cdir/snps/bone_related_traits/TSV/$name.tsv

    if grep -q "af_cases_EUR"  <(zcat $file | head -n 1) ; then
        MAF="af_cases_EUR"
    else
        MAF="af_EUR"
    fi
    echo '
    source ~/.bashrc
    conda activate ldsc
    python ~/tools/ldsc/munge_sumstats.py --sumstats '$cdir'/snps/bone_related_traits/TSV/'$name'.tsv.gz --n-min 0 --N-col an_EUR --frq '$MAF' --maf-min 0 --a1 alt --a2 ref --snp rsid  --signed-sumstats beta_EUR,0 --p pval_EUR --chunksize 500000 --out '$cdir'/snps/bone_related_traits/non_merge/'$name' 
    python ~/tools/ldsc/munge_sumstats.py --sumstats '$cdir'/snps/bone_related_traits/TSV/'$name'.tsv.gz --n-min 0 --N-col an_EUR --frq '$MAF' --maf-min 0 --a1 alt --a2 ref --snp rsid  --signed-sumstats beta_EUR,0 --p pval_EUR  --merge-alleles '$cdir'/w_hm3.snplist --chunksize 500000 --out '$cdir'/snps/bone_related_traits/merge/'$name'.merged ' > $cdir/snps/bone_related_traits/mung.$name.sh
    sed -i '1i#!/bin/bash' $cdir/snps/bone_related_traits/mung.$name.sh
    sbatch -c 7 --mem 8G -J munge.$name -o $cdir/snps/bone_related_traits/log/$name.log $cdir/snps/bone_related_traits/mung.$name.sh
done


zcat /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits/PDB_GWAS_results_2022_UKBiobank.txt.gz | sed '1d' | awk '{print "chr" $1 "\t" $2-1 "\t" $2 "\t" $3}' > PDB.bed
~/tools/liftOver PDB.bed ~/tools/hg19ToHg38.over.chain.gz PDB.hg38 PDB.unmap
for i in {1..22} "X" "Y" ; do
    intersectBed -wo -a PDB.hg38 -b /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz  > PDB.hg38.rsid.chr$i
done
cat PDB.hg38.rsid.chr* | awk '{if($2==$6 && $3==$7) print}' | cut -f1-4,8 > PDB.hg38.rsid
join <(cut -f4,5 PDB.hg38.rsid | sort -k1,1 ) <(zcat /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits/PDB_GWAS_results_2022_UKBiobank.txt.gz | cut -f3- | sort -k1,1 ) > PDB.gwas
sed -i "s/ /\t/g" PDB.gwas
sed -i "1iID\tSNP\tREF\tALT\tA1\tOR\tSE\tP" PDB.gwas
conda activate ldsc
python ~/ldsc/munge_sumstats.py --sumstats PDB.gwas --a1 A1 --a2 REF --snp SNP --N-cas 750 --N-con 1002 --no-alleles --out PDB
python ~/ldsc/munge_sumstats.py --sumstats PDB.gwas --a1 A1 --a2 REF --snp SNP --N-cas 750 --N-con 1002 --merge-alleles $cdir/w_hm3.snplist --chunksize 500000 --out PDB.merged
