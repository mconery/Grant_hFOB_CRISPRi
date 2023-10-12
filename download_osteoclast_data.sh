#!/bin/bash

#SBATCH -o download_osteoclast_data.log
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G

################################### Set Files ########################################

chip_dir=/mnt/isilon/sfgi/conerym/data/ChIP-seq/osteoclasts
atac_dir=/mnt/isilon/sfgi/conerym/data/ATAC-seq/osteoclasts

bae_atac_acc_file=$atac_dir/Bae_2022.ATAC.SRR_Acc_List.txt
bae_chip_acc_file=$chip_dir/Bae_2022.H3K27ac.SRR_Acc_List.txt
tesfaye_chip_acc_file=$chip_dir/Tesfaye_2023.H3K27ac.SRR_Acc_List.txt

######################################################################################

for file in $bae_atac_acc_file $bae_chip_acc_file $tesfaye_chip_acc_file; do 
	out_dir=$(dirname $file)
	while read accession; do 
		if [ ! -f "$out_dir/$accession.fastq" ]; then
			fastq-dump $accession -O $out_dir
		fi
	done < $file
done