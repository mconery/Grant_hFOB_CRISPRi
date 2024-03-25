#!/bin/bash

#SBATCH -o run_pygenome.log
#SBATCH --mem=16G
#SBATCH --time=6-00:00:00

#Set files
work_dir="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/osteoblasts_and_adipose"
successful_grnas="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/osteoblasts_and_adipose/hfob_crispri.sig_perturbations.tsv"
hfob_interaction_file="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/bone_cells/bone_mineral_density/ld_proxies.hFOBsDiff.with_pediatric.snpOE_interaction.both.txt"
hmsc_interaction_file="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/bone_cells/bone_mineral_density/ld_proxies.hMSC_Osteo.with_pediatric.snpOE_interaction.both.txt"
adipo_interaction_file="/mnt/isilon/sfgi/conerym/analyses/grant/varianttogene/hMSC_adipocytes/bmd/ld_proxies.hMSC_Adipo.with_pediatric.snpOE_interaction.both.txt"
gene_bed="/mnt/isilon/sfgi/pahlm/annotationFiles/genome_coords/gencodeV19.basic.annotation.bed12"

#Make a folder for slurm outputs
mkdir -p $work_dir/temp_shell_files
cd $work_dir/temp_shell_files

#Iterate over unique list of gRNAs
for grna in $(awk 'NR > 1 {print $1}' $successful_grnas | sort -u); do 

	#Identify genes for the grnas in order to name the directory
	genes=$(awk -v g=$grna '$1 == g {print $3}' $successful_grnas | sort -u | tr '\n' '_' | sed 's/.$//')
	#Make a directory for the locus
	mkdir -p $work_dir/$genes
	
	#Identify SNPs for the grnas to make the SNPs file
	snps=($(awk -v g=$grna '$1 == g {print $2}' $successful_grnas | sort -u | tr '\n' ' '))
	#Make the Targeted SNP file by looping over the SNPs and getting there locations from the proxy file
	>$work_dir/$genes/target_snp.bed
	for snp in ${snps[@]}; do
		awk -v s=$snp '$2 == s {print $4"\t"$5"\t"$5+1"\t"$2"\t1\t."}' $successful_grnas >> $work_dir/$genes/target_snp.bed 
	done
	
	### Make BMP2-hMSC Arcs ###
	#Create temp file with list of intervals for other ends (non-baits) from interaction file (This is for the arcs files)
	>$work_dir/$genes/BMP2.tmp.arcs
	>$work_dir/$genes/BMP2.interaction.arcs
	for snp in ${snps[@]}; do
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $hmsc_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF>7{print $0}' >> $work_dir/$genes/BMP2.tmp.arcs
		for count in {0..1}; do
			awk -v count=$count '{print $(1 + 3*count)"\t"$(2 + 3*count)"\t"$(3 + 3*count)"\t"$(7 + 3*count)"\t"$(8 + 3*count)"\t"$(9 + 3*count)"\t"$(13 + count)}' $work_dir/$genes/BMP2.tmp.arcs >> $work_dir/$genes/BMP2.interaction.arcs
		done
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $hmsc_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF==7{print $0}' >> $work_dir/$genes/BMP2.interaction.arcs
	done
	rm $work_dir/$genes/BMP2.tmp.arcs
	
	### Make hFOB Arcs ###
	#Create temp file with list of intervals for other ends (non-baits) from interaction file (This is for the arcs files)
	>$work_dir/$genes/hFOB.tmp.arcs
	>$work_dir/$genes/hFOB.interaction.arcs
	for snp in ${snps[@]}; do
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $hfob_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF>7{print $0}' >> $work_dir/$genes/hFOB.tmp.arcs
		for count in {0..1}; do
			awk -v count=$count '{print $(1 + 3*count)"\t"$(2 + 3*count)"\t"$(3 + 3*count)"\t"$(7 + 3*count)"\t"$(8 + 3*count)"\t"$(9 + 3*count)"\t"$(13 + count)}' $work_dir/$genes/hFOB.tmp.arcs >> $work_dir/$genes/hFOB.interaction.arcs
		done
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $hfob_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF==7{print $0}' >> $work_dir/$genes/hFOB.interaction.arcs
	done
	rm $work_dir/$genes/hFOB.tmp.arcs
	
	### Make Adipocyte Arcs ###
	#Create temp file with list of intervals for other ends (non-baits) from interaction file (This is for the arcs files)
	>$work_dir/$genes/adipo.tmp.arcs
	>$work_dir/$genes/adipo.interaction.arcs
	for snp in ${snps[@]}; do
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $adipo_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF>7{print $0}' >> $work_dir/$genes/adipo.tmp.arcs
		for count in {0..1}; do
			awk -v count=$count '{print $(1 + 3*count)"\t"$(2 + 3*count)"\t"$(3 + 3*count)"\t"$(7 + 3*count)"\t"$(8 + 3*count)"\t"$(9 + 3*count)"\t"$(13 + count)}' $work_dir/$genes/adipo.tmp.arcs >> $work_dir/$genes/adipo.interaction.arcs
		done
		awk -v snp=$snp '$2==snp {print $12"\t"$10"\t"$16}' $adipo_interaction_file | sed 's/:/\t/g' | sed 's/-/\t/g' | sed 's/&/\t/g'| awk 'NF==7{print $0}' >> $work_dir/$genes/adipo.interaction.arcs
	done
	rm $work_dir/$genes/adipo.tmp.arcs
	
	#Get positions of implicated genes to make sure to get full bounds of gene body
	genes_array=($(awk -v g=$grna '$1 == g {print $3}' $successful_grnas | sort -u | tr '\n' ' '))
	>$work_dir/$genes/imp_gene.bed
	for gene in ${genes_array[@]}; do
		awk -v tmp=$gene '$4==tmp' $gene_bed >> $work_dir/$genes/imp_gene.bed
	done
	
	#Get chromosome number
	chromosome=`awk 'NR==1 {print $1}' $work_dir/$genes/target_snp.bed `
	
	#Get minimum and maximum values for plotting
	min=`awk 'BEGIN{a=100000000000}{if ($2<0+a) a=$2} END{print a}' $work_dir/$genes/target_snp.bed`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($2<0+a) a=$2} END{print a}' $work_dir/$genes/BMP2.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($5<0+a) a=$5} END{print a}' $work_dir/$genes/BMP2.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($2<0+a) a=$2} END{print a}' $work_dir/$genes/hFOB.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($5<0+a) a=$5} END{print a}' $work_dir/$genes/hFOB.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($2<0+a) a=$2} END{print a}' $work_dir/$genes/adipo.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($5<0+a) a=$5} END{print a}' $work_dir/$genes/adipo.interaction.arcs`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($2<0+a) a=$2} END{print a}' $work_dir/$genes/imp_gene.bed`
	min=`awk -v MN=$min 'BEGIN{a=MN}{if ($7<0+a) a=$7} END{print a}' $work_dir/$genes/imp_gene.bed`
	max=`awk 'BEGIN{a=   0}{if ($3>0+a) a=$3} END{print a}' $work_dir/$genes/target_snp.bed`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($3>0+a) a=$3} END{print a}' $work_dir/$genes/BMP2.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($6>0+a) a=$6} END{print a}' $work_dir/$genes/BMP2.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($3>0+a) a=$3} END{print a}' $work_dir/$genes/hFOB.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($6>0+a) a=$6} END{print a}' $work_dir/$genes/hFOB.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($3>0+a) a=$3} END{print a}' $work_dir/$genes/adipo.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($6>0+a) a=$6} END{print a}' $work_dir/$genes/adipo.interaction.arcs`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($3>0+a) a=$3} END{print a}' $work_dir/$genes/imp_gene.bed`
	max=`awk -v MX=$max 'BEGIN{a=MX}{if ($8>0+a) a=$8} END{print a}' $work_dir/$genes/imp_gene.bed`
	
	#calculate increment and plot min and max (ends +/- 5% of spread)
	let plot_min=$(($min - $(($(($max-$min))/20))))
	let plot_max=$(($max + $(($(($max-$min))/20))))
	
	#Now get all genes that need to be plotted
	>$work_dir/$genes/temp.bed
	for gene in ${genes_array[@]}; do
		awk -v tmp=$gene '$4==tmp {print $1"\t"$2"\t"$3"\t"$4"\t1\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' $gene_bed >> $work_dir/$genes/temp.bed
		awk -v tmp=$gene -v min=$plot_min -v max=$plot_max -v chr=$chromosome '($1==chr && ($2>=min && $2<=max) || ($3>=min && $3<=max)) && $4!=tmp' $gene_bed >> $work_dir/$genes/temp.bed
	done
	sort -k1,1 -k2,2n -k3,3n -k4g -k5,5r $work_dir/$genes/temp.bed | uniq -f 5 > $work_dir/$genes/plot_gene.bed
	rm $work_dir/$genes/temp.bed
	

	#Copy tracks.ini file to directory
	cp $work_dir/tracks.ini $work_dir/$genes/tracks.ini

	#Execute plotting function
	echo "/home/conerym/miniconda3/bin/pyGenomeTracks  --tracks $work_dir/$genes/tracks.ini -o $work_dir/$genes/$genes.png --region chr$chromosome:$plot_min-$plot_max --fontSize 8 --dpi 1000" > $work_dir/temp_shell_files/make_plot.$genes.sh
	sed -i '1i#!/bin/bash' $work_dir/temp_shell_files/make_plot.$genes.sh
	sbatch --mem=64G -t 6:00:00 --job-name $genes $work_dir/temp_shell_files/make_plot.$genes.sh

done




