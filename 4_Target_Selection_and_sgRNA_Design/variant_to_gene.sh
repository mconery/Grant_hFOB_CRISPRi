#!/bin/bash

#SBATCH -o variant_to_gene.log
#SBATCH --mem=64G
#SBATCH --time=6-00:00:00

#Load modules
#module load R/4.0.2
#module load perl/5.26.1


#Pull needed locs and variables from input parameters
work_dir=$1
merged_proxies=$2
ATACseq_file=$3
ibed_1frag=$4
ibed_4frag=$5
out_prefix=$6

#Change directory
cd $work_dir

####### prepare proxy.all.bed #######
awk '{print $3"\t"$6"\t"$1"\t"$2"\t"$5"\t"$7}' $merged_proxies | sed 1d | awk 'BEGIN{OFS="\t"; FS="\t"}{chr="chr"$3; name=$1" "$2; end=$5; start=end-1; print chr, start, end, name}' | bedtools sort -i > proxy.bed

gencodeV19_all="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/source/unique_transcript_model.bed"
gencodeV19_coding="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/source/unique_transcript_model.protein_coding_only.bed"

#sort the gencode files
bedtools sort -i $gencodeV19_all > gencodeV19_all.sort.bed
bedtools sort -i $gencodeV19_coding > gencodeV19_coding.sort.bed

closest-features --closest --dist --delim '\t' proxy.bed gencodeV19_coding.sort.bed > tmp1.bed
closest-features --closest --dist --delim '\t' proxy.bed gencodeV19_all.sort.bed > tmp2.bed

perl ~/captureC/scripts/perl/makeFinalSnpBed.pl | sortBed -i - > proxy.all.bed

rm tmp1.bed tmp2.bed proxy.bed gencodeV19_coding.sort.bed gencodeV19_all.sort.bed


####### prepare rsquare.all.csv #######
awk '{print $3"\t"$6"\t"$1"\t"$2"\t"$5"\t"$7}' $merged_proxies | awk 'BEGIN{OFS=","; FS="\t"}{print $1,$2,$6}' > rsquare.all.csv


####### find proxy associated fragment #######
bed_1frag="/mnt/isilon/sfgi/referenceSequences/hg19/hg19_DpnII_frags.bed" # change if you use Arima
intersectBed -a proxy.all.bed -b $bed_1frag -wa -wb > proxy.all.1frag.txt

bed_4frag="/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/hg19_DpnII_4frag.bed" # change if you use Arima
intersectBed -a proxy.all.bed -b $bed_4frag -wa -wb > proxy.all.4frag.txt


####### incorporate ATAC-seq #######
zcat $ATACseq_file | awk '{print $1"\t"$2"\t"$3"\t"$4}' > atac.bed

####### FIND PROXY THAT COINCIDES IN OPEN CHROMATIN #######
intersectBed -a proxy.all.bed -b atac.bed -wa -u > openProxy.bed
intersectBed -a proxy.all.bed -b atac.bed -wa -wb > openProxy.txt


####### FIND OPEN PROXY ASSOCIATED FRAGMENT #######
bed_1frag="/mnt/isilon/sfgi/referenceSequences/hg19/hg19_DpnII_frags.bed" # change if you use Arima
intersectBed -a openProxy.bed -b $bed_1frag -wa -wb > openProxy.1frag.txt

bed_4frag="/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/hg19_DpnII_4frag.bed" # change if you use Arima
intersectBed -a openProxy.bed -b $bed_4frag -wa -wb > openProxy.4frag.txt


####### FIND ALL THE BAIT FRAGMENT IN OPEN CHROMATIN #######
bait_1frag="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/baits_1frag.bed" # change if you use Arima
bait_4frag="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/baits_4frag.bed" # change if you use Arima

intersectBed -a $bait_1frag -b atac.bed -wa -u > baitOpen.1frag.bed
intersectBed -a $bait_4frag -b atac.bed -wa -u > baitOpen.4frag.bed


####### FIND PROMOTER OCR #######
promoter_bed="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.promoter.bed"
bedtools intersect -a atac.bed -b $promoter_bed -wa -wb > openPromoter.txt


####### SNP LOOKUP SnpAtOe #######
## 1frag
baitmap_file="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap" # change if you use Arima
tad_anno_file="/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/hg19_DpnII_no-centromere.annot" # change if you use Arima
openSNPFrag_file=openProxy.1frag.txt
openBait_file=baitOpen.1frag.bed
output=$out_prefix.snpOE_interaction.1frag.txt 

perl ~/captureC/scripts/perl/getCisSnp2BaitInteractions.v2.pl \
$baitmap_file \
$openSNPFrag_file \
$ibed_1frag \
$output \
$tad_anno_file \
$openBait_file

# 4frag
baitmap_file="/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap" # change if you use Arima
tad_anno_file="/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/hg19_DpnII_4frag_no-centromere.annot" # change if you use Arima
openSNPFrag_file=openProxy.4frag.txt
openBait_file=baitOpen.4frag.bed
output=$out_prefix.snpOE_interaction.4frag.txt

perl ~/captureC/scripts/perl/getCisSnp2BaitInteractions.v2.pl \
$baitmap_file \
$openSNPFrag_file \
$ibed_4frag \
$output \
$tad_anno_file \
$openBait_file

# merge snp lookup result (add rsquare info)
rsquare_csv="rsquare.all.csv"
output_1frag=$out_prefix."snpOE_interaction.1frag.txt"
output_4frag=$out_prefix."snpOE_interaction.4frag.txt"
output_bothFrag=$out_prefix."snpOE_interaction.both.txt"
perl ~/captureC/scripts/perl/mergeSnpLookUp.v2.pl $rsquare_csv $output_1frag $output_4frag > $output_bothFrag


####### SNP LOOKUP SnpAtGenePromoter #######
openProxy_file="openProxy.txt"
openPromoter_file="openPromoter.txt"
rsquare_csv="rsquare.all.csv"
output=$out_prefix."snpOpenPromoter.txt"

perl ~/captureC/scripts/perl/openSnpAtPromoter.pl $openProxy_file $openPromoter_file $rsquare_csv > $output
