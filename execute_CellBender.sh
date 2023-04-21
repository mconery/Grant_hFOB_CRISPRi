#!/bin/bash

#SBATCH -o execute_CellBender.log
#SBATCH --time=2-00:00:00
#SBATCH --partition gpuq
#SBATCH --gres=gpu:p100:2
#SBATCH --mem-per-gpu=16G

################################### Set Directories ##################################

work_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs
out_dir=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellbender_outputs/learning_rate_0.00005
pool_ref=/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/Grant_hFOB_CRISPRi/cellbender_counts.txt

######################################################################################

#Activate environment
source activate cellbender

for i in {1..8}; do 
	table_row=$((i + 1))
	expected_cells=$(awk -v a="$table_row" 'NR==a{print $2}' $pool_ref)
	num_droplets=$(awk -v a="$table_row" 'NR==a{print $3}' $pool_ref)
	cellbender remove-background \
        	--input $work_dir/Pool_"$i"/raw_feature_bc_matrix.h5 \
                --output $out_dir/Pool_"$i"/cell_bender.output.h5 \
                --cuda \
                --expected-cells $expected_cells \
                --total-droplets-included $num_droplets \
		--learning-rate 0.00005 \
                --fpr 0.01 \
                --epochs 150
done
