#!/bin/bash

#SBATCH -o run_variant_to_gene.log
#SBATCH --mem=64G
#SBATCH --time=6-00:00:00


#Set files
work_dir="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/bone_cells"
bmd_proxies="/mnt/isilon/sfgi/conerym/data/gwas_summary_stats/clumped_sig/lead_snps_only/bmd_cojo_0.8_ld_proxies.with_pediatric.txt.ld"
ATACseq_file="/mnt/isilon/sfgi/pahlm/analyses/grant/atacSeq/bone_diff/hFOBsDiff/out/peak/macs2/pooled_rep/hFOBsDiff_rep1_R1.PE2SE.nodup.tn5_pooled.pf.pval0.01.300K.filt.narrowPeak.gz"
ibed_1frag="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff/merge/chicagoRes_1frag.ibed"
ibed_4frag="/mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/chicago/hFOBdiff/merge/chicagoRes_4frag.ibed"
out_prefix="ld_proxies.hFOBsDiff.with_pediatric"


#Call the pipeline
sbatch variant_to_gene.sh $work_dir/bone_mineral_density $bmd_proxies $ATACseq_file $ibed_1frag $ibed_4frag $out_prefix

