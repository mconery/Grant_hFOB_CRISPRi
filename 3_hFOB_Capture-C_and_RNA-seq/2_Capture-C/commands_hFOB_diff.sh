#Analysis of capture C for differentiated hFOBs

#Run Hicup

#Set up rawData folder
cd /mnt/isilon/sfgi/rawData/grant/captureC/Promoterome/hFOB_diff
ln -s /mnt/isilon/sfgi/novaseq_raw/*HNNHHDMXX*/Data/Intensities/BaseCalls/hFOB* .
R

#Rename fastq
files = list.files()
files = files[grepl("fastq.gz", files)]
new.files = gsub("_5DD_pCAPc","diff", files)
new.files = gsub("S._","", new.files)
new.files = gsub("_001","", new.files)

file.rename(files, new.files)

#Link files to my folder
cd /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/hFOBdiff_rep1
ln -s /mnt/isilon/sfgi/rawData/grant/captureC/Promoterome/hFOB_diff/*rep1* .

cd /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/hFOBdiff_rep2
ln -s /mnt/isilon/sfgi/rawData/grant/captureC/Promoterome/hFOB_diff/*rep2* .

cd /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/hFOBdiff_rep3
ln -s /mnt/isilon/sfgi/rawData/grant/captureC/Promoterome/hFOB_diff/*rep3* .

#PREP HICUP

#Run below script for each rep in parallel

#Rep1
PROMOTEROME_PARENT_DIR="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome"
REP="rep1" # change
PROMOTEROME_DIR="hFOBdiff_"$REP # change
FASTQ_R1="hFOBdiff_"$REP"_R1.fastq.gz" # change
FASTQ_R2="hFOBdiff_"$REP"_R2.fastq.gz" # change
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

#Rep2
PROMOTEROME_PARENT_DIR="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome"
REP="rep2" # change
PROMOTEROME_DIR="hFOBdiff_"$REP # change
FASTQ_R1="hFOBdiff_"$REP"_R1.fastq.gz" # change
FASTQ_R2="hFOBdiff_"$REP"_R2.fastq.gz" # change
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

#Rep3
PROMOTEROME_PARENT_DIR="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome"
REP="rep3" # change
PROMOTEROME_DIR="hFOBdiff_"$REP # change
FASTQ_R1="hFOBdiff_"$REP"_R1.fastq.gz" # change
FASTQ_R2="hFOBdiff_"$REP"_R2.fastq.gz" # change
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup


#####for each rep

#### CREATE FASTQ CHUNK FOR R1 AND R2
qsub -cwd -j y -o split_1_1.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R1 200000000 1 R1
qsub -cwd -j y -o split_1_2.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R2 200000000 1 R2

#############for each rep
#Zip files and run hicup

total_split=$(ls lane_1_R1* | wc -l)
total_split_0base=$((total_split-1))

### FOR EACH LANE, CREATE CHUNK DIRECTORIES
for i in `seq 1 $total_split`;do mkdir lane_1_${i};done

### GZIP FASTQ
for i in `seq -w 0 $total_split_0base`;do qsub -cwd -j y -o logs -N gzip_1_R1_${i}  ~/captureC/scripts/bash/gzip.sh lane_1_R1_${i};done
for i in `seq -w 0 $total_split_0base`;do qsub -cwd -j y -o logs -N gzip_1_R2_${i}  ~/captureC/scripts/bash/gzip.sh lane_1_R2_${i};done

#############for each rep

total_split=$(ls lane_1_R1* | wc -l)
### MOVE GZIP FASTQ TO DIR
perl ~/captureC/scripts/perl/mvSplits.pl $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup 1 R1 $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

perl ~/captureC/scripts/perl/mvSplits.pl $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup 1 R2 $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

### PREPARE JOB BASH AND CONFIG FILE FOR HICUP
#############for each rep

total_split=$(ls lane_1_R1* | wc -l)
### MOVE GZIP FASTQ TO DIR
perl ~/captureC/scripts/perl/mvSplits.pl $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup 1 R1 $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

perl ~/captureC/scripts/perl/mvSplits.pl $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup 1 R2 $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

### PREPARE JOB BASH AND CONFIG FILE FOR HICUP
perl ~/captureC/scripts/perl/writeHicupScripts.pl $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $total_split 1
perl ~/captureC/scripts/perl/writeHicupScripts.pl $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $total_split 2

### RUN HICUP IN THE QUEUE
perl ~/captureC/scripts/perl/runHicupScripts.pl  $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $total_split 1 8G


### CHECK HICUP RUNNING (IN DIFFERENT CHUNK) ARE SUCCESSFULLY COMPLETED
grep "HiCUP processing complete" lane*/hicup.log | wc -l
# check whether it matches split chunk number
total_split=$(ls -d lane_1_*/ | wc -l)
echo $total_split
total_split_0base=$((total_split-1))

### COMBINE CHUNKS AND RE-DO DEDUPLICATION
perl ~/captureC/scripts/perl/runHicupMergeDedup.pl 1 $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup 4G

### CONVERT TO CHICAGO INPUT (1 FRAGMENT RESOLUTION)
DIR_1FRAG=$REP"_1frag"
qsub -l h_vmem=32G -cwd -o bam2chicago_1frag.o -e bam2chicago_1frag.e ~/captureC/scripts/bash/bam2chicago.sh $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/chicago_1frag.baitmap /mnt/isilon/sfgi/programs/CHiCAGO/auxiliary/hg19_DpnII.rmap $DIR_1FRAG

### CONVERT TO CHICAGO INPUT (4 FRAGMENT RESOLUTION)
DIR_4FRAG=$REP"_4frag"
qsub -l h_vmem=32G -cwd -o bam2chicago_4frag.o -e bam2chicago_4frag.e ~/captureC/scripts/bash/bam2chicago.sh $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/chicago_4frag.baitmap /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/hg19_DpnII_4frag.rmap $DIR_4FRAG


########Run each time

### COLLECT INFO AND OUTPUT AS HICUP SUMMARY
total_split=$(ls -d lane_1_*/ | wc -l)
total_split_0base=$((total_split-1))

deduplicator_FILE=$(ls hicup_deduplicator_summary*)
SUMMARY_FILE=$PROMOTEROME_DIR"_hicupSummary.txt"
DIR_1FRAG=$REP"_1frag"
chinput_FILE=$DIR_1FRAG".chinput"
perl ~/captureC/scripts/perl/collectHicupResults.pl 1 $total_split $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$deduplicator_FILE $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$DIR_1FRAG/$chinput_FILE > $SUMMARY_FILE

### COPY the the first 5 lines of EBV_Bcells_3_hicupSummary.txt to excel template to make figure summaryChicago

# GET 4FRAG CAPTURED
DIR_4FRAG=$REP"_4frag"
chinput_FILE=$DIR_4FRAG".chinput"
TAG_NUM=$(sed '6q;d' $PROMOTEROME_DIR"_hicupSummary.txt" | awk '{print $24}')
OUTPUT_FILE=$PROMOTEROME_DIR"_4frag_captured.txt"
perl ~/captureC/scripts/perl/captureStatsFromChinput.pl $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$DIR_4FRAG/$chinput_FILE $TAG_NUM > $OUTPUT_FILE


# CLEAN-UP
rm -f lane_*_*/*fastq.gz
rm -f lane_*/*fq.gz
rm -f lane_*/*pair.bam

#####################################################
####CHICAGO#####################################################
#####################################################

main_dir="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff"
cd $main_dir
mkdir 1frag
mkdir 4frag

################## PREPARE CHICAGO INPUT FILE ################

cd $main_dir/1frag
files=($(ls /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/hFOBdiff*/hicup/rep*_1frag/*_1frag.chinput | grep -v "lane1"))
for file in ${files[@]}; do
        new_file=$(echo $file | cut -d "/" -f 10 | cut -d "_" -f 2)
        new_file=$new_file"_1frag.chinput"
        # echo -e "$file\t$new_file"
        ln -s $file $new_file
done

cd $main_dir/4frag
files=($(ls /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/hFOBdiff*/hicup/rep*_4frag/*_4frag.chinput | grep -v "lane1"))
for file in ${files[@]}; do
        new_file=$(echo $file | cut -d "/" -f 10 | cut -d "_" -f 2)
        new_file=$new_file"_4frag.chinput"
        # echo -e "$file\t$new_file"
        ln -s $file $new_file
done

##############
cd $main_dir
mkdir $main_dir/features
cd $main_dir/features
ln -s /mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/bone_diff/hFOBsDiff/out/peak/macs2/idr/conservative_set/hFOB_diff_rep3-rep4.IDR0.1.filt.narrowPeak.gz feature_atacseq.list

# RUN CHICAGO 1FRAG AND 4FRAG
featureList="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff/features/feature_atacseq.list"
cd $main_dir/1frag
chinput=$(ls /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff/1frag/rep*_1frag.chinput | tr "\n" "," | sed "s/,$//")
qsub -l m_mem_free=300G -l h_vmem=300G -cwd -o runChicago.o -e runChicago.e ~/captureC/scripts/bash/runChicago.sh /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/ $featureList $chinput chicagoRes

cd $main_dir/4frag
chinput=$(ls /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff/4frag/rep*_4frag.chinput | tr "\n" "," | sed "s/,$//")
qsub -l m_mem_free=200G -l h_vmem=200G -cwd -o runChicago.o -e runChicago.e ~/captureC/scripts/bash/runChicago.sh /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt /mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/ $featureList $chinput chicagoRes

cd $main_dir/1frag
# qsub -cwd -l h_vmem=200G /home/suc1/scripts/bash/exportIbed.sh
module load R/3.3.2
R
vi /home/suc1/scripts/R/exportIbed.R
# library(Chicago)
# cd=readRDS("chicagoRes/data/chicagoRes.Rds")
# exportResults(cd,outfileprefix="chicagoRes/data/chicagoRes",format=c("interBed"))

cd $main_dir/4frag
module load R/3.3.2
R
vi /home/suc1/scripts/R/exportIbed.R
library(Chicago)
cd=readRDS("chicagoRes/data/chicagoRes.Rds")
# library(Chicago)
# cd=readRDS("chicagoRes/data/chicagoRes.Rds")
# exportResults(cd,outfileprefix="chicagoRes/data/chicagoRes",format=c("interBed"))
