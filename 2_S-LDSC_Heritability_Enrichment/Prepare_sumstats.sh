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


ln -s /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits ./

# Paget's disease
zcat bone_traits/PDB_GWAS_results_2022_UKBiobank.txt.gz | sed '1d' | awk '{print "chr" $1 "\t" $2-1 "\t" $2 "\t" $3}' > PDB.bed
~/tools/liftOver PDB.bed ~/tools/hg19ToHg38.over.chain.gz PDB.hg38 PDB.unmap
for i in {1..22} "X" "Y" ; do
        sbatch -J chr$i  -o tmp$i -c 8 --mem 20G ~/misc_bash/get_RSID.sh /mnt/isilon/sfgi/trangk/analyses/grant/ldsr/snps/bone_related_traits PDB.hg38 $i
done
cat PDB.hg38.rsid.chr* | awk '{if($2==$6 && $3==$7) print}' | cut -f1-4,8 > PDB.hg38.rsid
join <(cut -f4,5 PDB.hg38.rsid | sort -k1,1 ) <(zcat PDB_GWAS_results_2022_UKBiobank.txt.gz | cut -f3- | sort -k1,1 ) > PDB.gwas
sed -i "s/ /\t/g" PDB.gwas
sed -i "1iID\tSNP\tREF\tALT\tA1\tOR\tSE\tP" PDB.gwas
python ~/tools/ldsc/munge_sumstats.py --sumstats PDB.gwas --a1 A1 --a2 REF --snp SNP --N-cas 750 --N-con 1002 --no-alleles --out PDB
python ~/tools/ldsc/munge_sumstats.py --sumstats PDB.gwas --a1 A1 --a2 REF --snp SNP --N-cas 750 --N-con 1002 --merge-alleles $cdir/w_hm3.snplist --chunksize 500000 --out PDB.merged

# Bone mineral density
zcat bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz | sed '1d' | cut -f2 | grep "rs" |  sed "s/rs//g" > tmp
python ~/tools/LiftRsNumber.py tmp > tmp.lift
paste -d'\t' <( cut -f2 tmp.lift | sed "s/^ /rs/g" ) <(zcat bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz| grep "rs" | cut -f1,3- ) > BMD.sumst
awk '{print "chr" $3 "\t" $4-1 "\t" $4 "\t" $1}' BMD.sumst  > tmp
~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz tmp.lift tmp.non
grep -v "#" tmp.non | cut -f1,4 > tmp
rm tmp.non
for chr in $(cut -f1 tmp | sort -u) ; do
        i=$(echo $chr | sed "s/chr//g")
        grep -w "^$chr" tmp | cut -f2 > snps
        zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz | grep -w -f snps | cut -f1-4 >> tmp.lift
done
chr="chr23"
i=X
grep -w "^$chr" tmp | cut -f2 > snps
zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz | grep -w -f snps | cut -f1-4 >> tmp.lift
for i in {1..22} "X" ; do
        chr="chr"$i"_"
        for pat in $(grep "^$chr" tmp.lift | cut -f1 | sort -u) ; do
                sed  -i "s/$pat/chr$i/g"  tmp.lift
        done
done
echo "CHROM POS SNP EA NEA EAF INFO BETA SE P P.I P.NI N" > BMD
join -1 3 -2 1 -o 1.1,1.2,1.3,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11 <(cut -f1,3,4 tmp.lift | sort -k3,3 ) <(cut -f1,5- BMD.sumst| sort -k1,1 ) >> BMD
zcat bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz | grep -v "rs" | awk '{print "chr" $3 "\t" $4-1 "\t" $4 "\t" $1}' | sed '1d' > tmp
sed -i "s/chr23/chrX/g" tmp
~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz tmp.lift tmp.non
for i in {1..22} "X" ; do
        chr="chr"$i"_"
        for pat in $(grep "^$chr" tmp.lift | cut -f1 | sort -u) ; do
                sed  -i "s/$pat/chr$i/g"  tmp.lift
        done
done
for chr in $(cut -f1 tmp.lift | sort -u) ; do
        i=$(echo $chr | sed "s/chr//g")
        grep -w "^$chr" tmp.lift | cut -f3 > pos
        grep -w -f pos <(zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz ) | cut -f1,3,4 >> tmp2
done
join -1 2 -2 4 -o 1.1,1.3,1.4,2.3 <(awk '{$2=$1"_"$3 ; print $0}' tmp.lift | sort -k2,2) <(awk '{$4=$1"_"$2; print $0}' tmp2 | sort -k4,4) | sort -u > tmp
join -1 3 -2 1 -o 1.1,1.2,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11 <(sort -k3,3 tmp) <(zcat bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz | sed '1d' | grep -v "rs" | cut -f1,5- | sort -k1,1) | sort -u >> BMD
sed -i "s/ /\t/g" BMD
sed 's/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/g;s/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/g' BMD > BMD.capped
gzip BMD
gzip BMD.capped
python ~/tools/ldsc/munge_sumstats.py --sumstats BMD.capped.gz --a1 EA --a2 NEA --snp SNP --n-min 1 --maf-min 0.0001 --info-min 0 --p "P.NI" --ignore "P" --no-alleles --out non_merge/BMD
python ~/tools/ldsc/munge_sumstats.py --sumstats BMD.capped.gz --a1 EA --a2 NEA --snp SNP --n-min 1 --maf-min 0.0001 --info-min 0 --p "P.NI" --ignore "P" --merge-alleles $cdir/w_hm3.snplist --chunksize 500000 --out merge/BMD

# Bone fracture
zcat bone_traits/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz |  sed '1d' | cut -f2 | grep "rs" |  sed "s/rs//g" > tmp
python ~/tools/LiftRsNumber.py tmp > tmp.lift
paste -d'\t' <( cut -f2 tmp.lift | sed "s/^ /rs/g" ) <(zcat bone_traits/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz | grep "rs" | cut -f1,3- ) > BF.sumst
awk '{print "chr" $3 "\t" $4-1 "\t" $4 "\t" $1}' BF.sumst > tmp
~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz tmp.lift tmp.non
grep -v "#" tmp.non | cut -f1,4 > tmp
rm tmp.non
for chr in $(cut -f1 tmp | sort -u) ; do
        i=$(echo $chr | sed "s/chr//g")
        grep -w "^$chr" tmp | cut -f2 > snps
        zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz | grep -w -f snps | cut -f1-4 >> tmp.lift
done
chr="chr23"
i=X
grep -w "^$chr" tmp | cut -f2 > snps
zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz | grep -w -f snps | cut -f1-4 >> tmp.lift
for i in {1..22} "X" ; do
        chr="chr"$i"_"
        for pat in $(grep "^$chr" tmp.lift | cut -f1 | sort -u) ; do
                sed  -i "s/$pat/chr$i/g"  tmp.lift
        done
done
echo "CHROM POS SNP ALLELE1 ALLELE0 A1FREQ INFO logOR logOR.SE OR L95 U95 P P.I P.NI N" > BF
join -1 3 -2 1 -o 1.1,1.2,1.3,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(cut -f1,3,4 tmp.lift | sort -k3,3 ) <(cut -f1,5- BF.sumst| sort -k1,1 ) >> BF
zcat bone_traits/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz |  sed '1d' | grep -v "rs"  | awk '{print "chr" $3 "\t" $4-1 "\t" $4 "\t" $1}' > tmp
sed -i "s/chr23/chrX/g" tmp
~/tools/liftOver tmp ~/tools/hg19ToHg38.over.chain.gz tmp.lift tmp.non
for i in {1..22} "X" ; do
        chr="chr"$i"_"
        for pat in $(grep "^$chr" tmp.lift | cut -f1 | sort -u) ; do
                sed  -i "s/$pat/chr$i/g"  tmp.lift
        done
done
for chr in $(cut -f1 tmp.lift | sort -u) ; do
        i=$(echo $chr | sed "s/chr//g")
        grep -w "^$chr" tmp.lift | cut -f3 > pos
        grep -w -f pos <(zcat /mnt/isilon/sfgi/trangk/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz ) | cut -f1,3,4 >> tmp2
done
join -1 2 -2 4 -o 1.1,1.3,1.4,2.3 <(awk '{$2=$1"_"$3 ; print $0}' tmp.lift | sort -k2,2) <(awk '{$4=$1"_"$2; print $0}' tmp2 | sort -k4,4) | sort -u > tmp
join -1 3 -2 1 -o 1.1,1.2,1.4,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 <(sort -k3,3 tmp) <(zcat bone_traits/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz | sed '1d' | grep -v "rs" | cut -f1,5- | sort -k1,1) | sort -u >> BF
sed 's/ /\t/g;s/[0-9].[0-9]E-[1-9][0-9][0-9][0-9]/1.0E-300/g;s/[0-9].[0-9]E-[3-9][0-9][0-9]/1.0E-300/g' BF > BF.capped
gzip BF.capped
python ~/tools/ldsc/munge_sumstats.py --sumstats BF.capped.gz --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --snp SNP --n-min 1 --maf-min 0.0001 --info-min 0 --p "P.I" --ignore "P" --no-alleles --out non_merge/BF
python ~/tools/ldsc/munge_sumstats.py --sumstats BF.capped.gz --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --snp SNP --n-min 1 --maf-min 0.0001 --info-min 0 --p "P.I" --ignore "P" --merge-alleles $cdir/w_hm3.snplist --chunksize 500000 --out merge/BF
