#Star commands
#hFOB undiff
source ~/.bashrc
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_1_N_R1.fastq.gz  hFOB_1_N_R2.fastq.gz

source ~/.bashrc
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_2_N_R1.fastq.gz  hFOB_2_N_R2.fastq.gz

source ~/.bashrc
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_3_N_R1.fastq.gz  hFOB_3_N_R2.fastq.gz

#hFOB diff
cd /mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/STAR/hFOB_diff/
source ~/.bashrc 
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_diff_rep1_R1.fastq.gz hFOB_diff_rep1_R2.fastq.gz

cd /mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/STAR/hFOB_diff/
source ~/.bashrc 
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_diff_rep2_R1.fastq.gz hFOB_diff_rep2_R2.fastq.gz

cd /mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/STAR/hFOB_diff/ 
source ~/.bashrc 
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs \
--readFilesCommand zcat --readFilesIn hFOB_diff_rep3/hFOB_diff_rep3_R1.fastq.gz hFOB_diff_rep3_R2.fastq.gz 

#HTseq
dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/HTSeq"
cd $dir
############# HTseq RETRIVE COUNT ##################
# run htseq directly from STAR output (GTF file without gene feature ribosome)
STAR_dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/STAR"
htseq_dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/HTSeq"

mkdir $htseq_dir/htseq_count
mkdir $htseq_dir/scripts

cd $htseq_dir
bam_files=($(ls $STAR_dir/*/*/Aligned.*.bam))
for bam_file in ${bam_files[@]}; do
        condition=$(dirname $bam_file | rev | cut -d "/" -f 1,2|sed 's/\//_/g'| rev)
        # echo $condition
        output=$condition".htseq.txt"
        job_file=$condition"_htseq.sh"
        echo "module load python/2.7
htseq-count -f bam -r name -s reverse -t exon -m intersection-strict $bam_file /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $htseq_dir/htseq_count/$output
" > $htseq_dir/scripts/$job_file
done

## SUBMIT JOBS
cd $htseq_dir/scripts/
files=($(ls *_htseq.sh))
for file in ${files[@]}; do
        qsub -cwd -l h_vmem=8G $file
done

#Gather files and make read count table
perl ~/functionalGenomics/scripts/perl/mergeReadCountFile.pl <(ls htseq_count/*.htseq.txt) > htseq_count.txt

