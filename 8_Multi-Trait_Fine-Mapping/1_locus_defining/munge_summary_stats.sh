#!/bin/bash

#SBATCH -o munge_summary_stats.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=96G

#Munge BMD
zcat /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11}' | sed 's/SNPID/SNP/g' | awk '{ if (NR == 1) {print $0} else {$7 = -1 * $7; print $3":"$4"_"$5"_"$6"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5"\t"$7"\t"$8"\t"$9} }' | sed 's/23:/X:/g' | gzip -c > /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats/BMD.sumstats.gz

#Munge Fracture
zcat /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/bone_traits/Biobank2-British-FracA-As-C-Gwas-SumStats.txt.gz | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9"\t"$10"\t"$14}' | sed 's/SNP/RSID/g' | sed 's/RSIDID/SNP/g' | sed 's/ALLELE1/EA/g' | sed 's/ALLELE0/NEA/g' | sed 's/logOR.SE/SE/g' | sed 's/logOR/BETA/g' | awk '{ if (NR == 1) {print $0} else {$7 = -1 * $7; print $3":"$4"_"$5"_"$6"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5"\t"$7"\t"$8"\t"$9} }' | sed 's/23:/X:/g' | gzip -c > /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats/BF.sumstats.gz

#Munge Pan-UKBB
pop=EUR; for file in $(ls /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/pan_ukbb_raw/*.tsv.gz); do trait=$(basename $file | sed 's/.tsv.gz//g'); declare -a array=($(zcat $file | paste -d '\t' - /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/pan_ukbb_raw/full_variant_qc_metrics.txt | head -n 1)); declare -a array2=($(echo "varid rsid chr pos ref alt  beta_$pop se_$pop neglog10_pval_$pop")); declare -a array3=($(echo "varid rsid chr pos ref alt beta se pval")); for (( i=0; i<${#array[@]}; i++ )); do for (( j=0; j<${#array2[@]}; j++ )); do if [ ${array2[$j]} == ${array[$i]} ]; then temp=$((i + 1)); declare ${array3[$j]}=$temp; fi; done; done; zcat $file | paste -d '\t' - /mnt/isilon/sfgi/conerym/data/gwas_summary_stats/pan_ukbb_raw/full_variant_qc_metrics.txt | awk -F "\t" -v x=$varid -v y=$rsid -v a=$chr -v b=$pos -v c=$alt -v d=$ref -v f=$beta -v g=$se -v h=$pval '{print $x"\t"$y"\t"$a"\t"$b"\t"$c"\t"$d"\t"$f"\t"$g"\t"$h}' | awk '{ if (NR == 1) {print "SNP\tRSID\tCHR\tBP\tEA\tNEA\tBETA\tSE\tP"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"10^(-$9)} }' | awk '{ if ($3 == "X") {print $1"\t"$2"\t23\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} else {print $0} }' | gzip -c > /mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/munged_summary_stats/$trait.sumstats.gz; done