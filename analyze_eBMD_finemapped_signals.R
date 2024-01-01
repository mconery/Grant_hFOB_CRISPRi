################################################################################

#analyze_eBMD_finemapped_signals.R

#The purpose of this script is to identify the targets for the full-scale
#CRISPRi screen based on the interesections of fine-mapped variants with ATAC-
# and H3K27ac ChIP-seq peaks. The script runs with the following modules:

#  R/4.2.3

################################################################################

# 1) Load libraries / set file locations and global variables ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(stringr)
library(pbapply)
library(bedr)
library(tidyverse)
library(cowplot)
library(reshape2)
library(ggpubr)

#Set File Locations 
inp_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection/kanai_fine-mapped_eBMD_variants.atac_chip_peaks.bed"
vep_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection/BMD_variant_VEP_annotations.txt"
vep_cat_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/Grant_hFOB_CRISPRi/vep_categories.prioritized.txt"
gene_ann_loc <- "/mnt/isilon/sfgi/conerym/data/genomes/gencode/gencode.v19.annotation.gtf.gz"
pilot_guides_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/all_targeting_guides.bed"
rna_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection/hFOB_tpm.txt"
plot_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection/plots/"
out_file <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/full_hFOB_screen_bone/target_selection/targets.no_exon.no_exon_cs.no_exon_prox.txt"

#Set Global Variables
repression_radius = 1000 #in bp
expression_thresh = 1 #in tpm
express_gene_range = 1000000 #in bp

# 2) Read in Files ====

#Read in main file 
inp_raw <- read.csv(inp_loc, header = FALSE, sep = "\t")
#Set column names 
colnames(inp_raw) <- c("chromosome", "start", "end", "variant", "rsid", "allele1", "allele2", "minorallele",
                       "cohort", "model_marginal", "method", "trait", "region", "maf", "beta_marginal", 
                       "se_marginal", "chisq_marginal", "pip", "cs_id", "beta_posterior", "sd_posterior", 
                       "LD_HWE", "LD_SV", "peak_type_num", "peak_chr", "peak_start", "peak_end", "peak_name", 
                       "peak_score", "peak_strand", "signalValue", "log10pValue", "log10qValue", "point_source",
                       "peak_type", "overlap")

#Read in vep and gene annotation files
vep_raw <- read.csv(vep_loc, header = TRUE, sep = "\t")
vep_cat_raw <- read.csv(vep_cat_loc, header = TRUE, sep = "\t")
rownames(vep_cat_raw) <- vep_cat_raw$VEP.Annotation
gene_ann_raw <- read.csv(gene_ann_loc, header = FALSE, sep = "\t")
colnames(gene_ann_raw) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
rna_raw <- read.csv(rna_loc, header = TRUE, sep = "\t")

#Read in pilot guides file
pilot_guides_raw <- read.table(pilot_guides_loc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 3) Filter for variants in distal non-coding regions ====

#Double Check Basic Statistics from awk Test
length(unique(inp_raw$rsid)) #3,590 unique rsids in either peak type
length(unique(inp_raw$variant)) #3,591 unique variants in either peak type
length(unique(inp_raw$region)) #265 unique regions
length(unique(inp_raw[which(inp_raw$method == "SUSIE"), "rsid"])) #2,448 rsids mapped by SuSiE
length(unique(inp_raw[which(inp_raw$method == "FINEMAP"), "rsid"])) #2,972 rsids mapped by Fine-Map

#Append coding variant info
inp_raw <- cbind.data.frame(inp_raw, 
                            vep_id=paste0(str_replace(inp_raw$chromosome, "chr", ""), ":", inp_raw$end, "-", 
                                   inp_raw$end + nchar(inp_raw$allele1) - 1))
vep_raw <- vep_raw[!(duplicated(vep_raw$Location)),]
rownames(vep_raw) <- vep_raw$Location
inp_raw <- cbind.data.frame(inp_raw, vep_ann=vep_raw[inp_raw$vep_id,"Consequence"],
                            coding=vep_cat_raw[vep_raw[inp_raw$vep_id,"Consequence"],"Coding"])
#Append cs unique id column onto inp_raw
inp_raw <- cbind.data.frame(inp_raw, unique_cs_id=paste(inp_raw$region, inp_raw$method, inp_raw$cs_id, sep="."))

#Create a dataframe that makes a single row for each unique rsid and combines info on peak types
temp <- unique(inp_raw$rsid)
check_methods = function(rsid, inp_raw, vep_cat_raw){
  temp2=inp_raw[which(inp_raw$rsid == rsid),]
  #First look to combine the fine-mapping info
  if (length(unique(temp2$method)) > 1) {
    check <- c(cbind.data.frame(temp2[which.max(temp2$pip),c("chromosome", "start", "end", "rsid")],
                                unlist(t(colMeans(temp2[,c("pip", "beta_posterior")]))), 
                                temp2[which.max(temp2$pip),c("region", "vep_ann", "coding")],"BOTH"))
    check <- c(check,vapply(sort(unique(temp2$unique_cs_id)), check_cs_coding, FUN.VALUE = numeric(2), inp_raw=inp_raw, vep_cat_raw=vep_cat_raw))
  } else {
    check <- c(cbind.data.frame(temp2[which.max(temp2$pip),c("chromosome", "start", "end", "rsid")],
                                unlist(t(colMeans(temp2[,c("pip", "beta_posterior")]))), 
                                temp2[which.max(temp2$pip),c("region", "vep_ann", "coding", "method")]))
    if (unique(temp2$method) == "FINEMAP") {
      check <- c(check, check_cs_coding(unique(temp2$unique_cs_id), inp_raw=inp_raw, vep_cat_raw=vep_cat_raw), NA, NA)
    } else {
      check <- c(check, NA, NA, check_cs_coding(unique(temp2$unique_cs_id), inp_raw=inp_raw, vep_cat_raw=vep_cat_raw))
    }
  } 
  #Then check the types of chip-seq peak intersections
  if (length(unique(temp2$peak_type)) > 1) {
    check <- c(check, "BOTH")
  } else {
    check <- c(check, temp2$peak_type[1])
  }
  names(check) <- NULL
  return(as.data.frame(t(check)))
}
#Check cs for coding variants
check_cs_coding <- function(unique_cs_id, inp_raw, vep_cat_raw){
  temp3 <- vep_cat_raw[inp_raw[inp_raw$unique_cs_id == unique_cs_id,"vep_ann"],]
  return(c(nrow(temp3), max(temp3[,"Coding"])))
}
temp2 <- pblapply(temp, check_methods, inp_raw=inp_raw, vep_cat_raw = vep_cat_raw)
combined_df <- as.data.frame(matrix(unlist(temp2), ncol = 15, byrow = TRUE))
colnames(combined_df) <- c("chr", "start", "end", "rsid", "mean_pip", "mean_beta_posterior", "region", 
                           "vep_ann", "coding", "method", "finemap_cs_size", "finemap_cs_coding", 
                           "susie_cs_size", "susie_cs_coding", "peak_type")
rownames(combined_df) <- temp
#Append on a label for the number of variants in the region
temp <- as.data.frame(table(combined_df[,c("region")]))
rownames(temp) <- temp$Var1
combined_df <- cbind.data.frame(combined_df, var_per_region=temp[combined_df$region, "Freq"])
#Recode finemap coding columns
combined_df$finemap_cs_coding = as.numeric(combined_df$finemap_cs_coding)
combined_df$susie_cs_coding = as.numeric(combined_df$susie_cs_coding)
#Recode start and end columns
combined_df$start <- as.numeric(combined_df$start)
combined_df$end <- as.numeric(combined_df$end)

#Split the gene annotations table into genes and exons
gene_ann_gene <- gene_ann_raw[which(gene_ann_raw$feature == "gene"),]
gene_ann_exon <- gene_ann_raw[which(gene_ann_raw$feature == "exon"),]
max_zero <- function(value){return(max(value, 0))} #Function to prevent neagative start positions
#Reorder gencode gene annotations into bed file format and pad out start and end positions by 2x the repression radius
temp <- str_locate(string = gene_ann_gene$attribute, "gene_name")[,"start"]
gene_ann_prom_org <- cbind.data.frame(chr=gene_ann_gene[,c("chr")], 
                                      start = vapply(ifelse(gene_ann_gene$strand == "+", 
                                                          gene_ann_gene$start - (2 * repression_radius), 
                                                          gene_ann_gene$end - (2 * repression_radius)),
                                                     FUN = max_zero, FUN.VALUE = numeric(1)),
                                      end = ifelse(gene_ann_gene$strand == "+", 
                                                        gene_ann_gene$start + (2 * repression_radius), 
                                                        gene_ann_gene$end + (2 * repression_radius)),
                                     score=gene_ann_gene[,c("score")], 
                                     name=substr(gene_ann_gene$attribute, temp + 10, str_locate(substr(gene_ann_gene$attribute, temp, nchar(gene_ann_gene$attribute)), ";") + temp - 2),
                                     gene_ann_gene[,c("strand", "feature", "frame", "attribute")])
temp <- str_locate(string = gene_ann_exon$attribute, "gene_name")[,"start"]
gene_ann_exon_org <- cbind.data.frame(chr=gene_ann_exon[,c("chr")], 
                                      start = vapply(gene_ann_exon$start - (2 * repression_radius),
                                                     FUN = max_zero, FUN.VALUE = numeric(1)),
                                      end = gene_ann_exon$end + (2 * repression_radius),
                                      score=gene_ann_exon[,c("score")], 
                                      name=substr(gene_ann_exon$attribute, temp + 10, str_locate(substr(gene_ann_exon$attribute, temp, nchar(gene_ann_exon$attribute)), ";") + temp - 2),
                                      gene_ann_exon[,c("strand", "feature", "frame", "attribute")])
remove(temp)

#Merge gene annotations
gene_ann_prom_merge <- bedr.merge.region(gene_ann_prom_org)
gene_ann_exon_merge <- bedr.merge.region(gene_ann_exon_org)
combined_df_merge <- bedr.merge.region(combined_df)
#Sort bed files
gene_ann_prom_merge <- bedr.sort.region(gene_ann_prom_merge, method = "lexicographical")
gene_ann_exon_merge <- bedr.sort.region(gene_ann_exon_merge, method = "lexicographical")
combined_df_merge <- bedr.sort.region(combined_df_merge, method = "lexicographical")
#Intersect the gencode and the sgrna target files bedtools style
gene_ann_prom_intersect <- bedr(
  input = list(a = combined_df_merge, b = gene_ann_prom_merge), 
  method = "intersect", 
  params = "-c -sorted"
)
gene_ann_exon_intersect <- bedr(
  input = list(a = combined_df_merge, b = gene_ann_exon_merge), 
  method = "intersect", 
  params = "-c -sorted"
)

#De-merge the lines in the intersect bed files
split_intersect_line <- function(line_intersect){
  if (gregexpr(line_intersect[4], ",") == TRUE) {
    temp = unlist(str_split(line_intersect[4], ","))
    return(lapply(temp, recombine_line, line_intersect=line_intersect))
  } else { return(line_intersect) }
}
recombine_line <- function(spot_4, line_intersect) {return(c(line_intersect[1:3], spot_4, line_intersect[5]))}
prom_intersect_df <- as.data.frame(matrix(unlist(apply(gene_ann_prom_intersect, split_intersect_line, MARGIN = 1)), byrow = TRUE, ncol = 5))
exon_intersect_df <- as.data.frame(matrix(unlist(apply(gene_ann_exon_intersect, split_intersect_line, MARGIN = 1)), byrow = TRUE, ncol = 5))
rownames(prom_intersect_df) <- prom_intersect_df$V4
rownames(exon_intersect_df) <- exon_intersect_df$V4

#Append on promoter/exon info to combined_df
combined_df <- cbind.data.frame(combined_df,
                                promoter_prox=ifelse(prom_intersect_df[combined_df$rsid, "V5"] > 0, TRUE, FALSE),
                                exon_prox=ifelse(prom_intersect_df[combined_df$rsid, "V5"] > 0, TRUE, FALSE))

### Count number of targetable variants in different levels of filtering ###
nrow(combined_df[combined_df$coding == 0,]) #3,530 non-exonic variants
nrow(combined_df[combined_df$coding == 0 & 
                   ifelse(is.na(combined_df$finemap_cs_coding), TRUE, combined_df$finemap_cs_coding == 0) &
                   ifelse(is.na(combined_df$susie_cs_coding), TRUE, combined_df$susie_cs_coding == 0),]) #2,955 non-exonic variants not implicated with an exonic variant via either method
nrow(combined_df[combined_df$coding == 0 & 
                   (combined_df$exon_prox == FALSE),]) #2,751 non-exonic variants that are >2kb from an exon
nrow(combined_df[combined_df$coding == 0 & 
                   (combined_df$exon_prox == FALSE) & 
                   ifelse(is.na(combined_df$finemap_cs_coding), TRUE, combined_df$finemap_cs_coding == 0) &
                   ifelse(is.na(combined_df$susie_cs_coding), TRUE, combined_df$susie_cs_coding == 0),]) #2,339 non-exonic variants independent from exonic variants that are >2kb from an exon
nrow(combined_df[combined_df$coding == 0 & 
                   (combined_df$exon_prox == FALSE) & 
                     ifelse(is.na(combined_df$finemap_cs_coding), TRUE, combined_df$finemap_cs_coding == 0) &
                     ifelse(is.na(combined_df$susie_cs_coding), TRUE, combined_df$susie_cs_coding == 0) &
                   (combined_df$promoter_prox == FALSE),]) # This is also 2,339 since there's an exon at the start of the gene!

# 4) Calculate the Number of gRNAs Needed to Hit Targets ====

#Split the combined_df up into different levels of exclusion
combined_no_exon <- combined_df[combined_df$coding == 0,]
combined_no_exon_cs <- combined_df[combined_df$coding == 0 & 
                                     ifelse(is.na(combined_df$finemap_cs_coding), TRUE, combined_df$finemap_cs_coding == 0) &
                                     ifelse(is.na(combined_df$susie_cs_coding), TRUE, combined_df$susie_cs_coding == 0),]
combined_no_exon_cs_prox <- combined_df[combined_df$coding == 0 & 
                                          (combined_df$exon_prox == FALSE) & 
                                          ifelse(is.na(combined_df$finemap_cs_coding), TRUE, combined_df$finemap_cs_coding == 0) &
                                          ifelse(is.na(combined_df$susie_cs_coding), TRUE, combined_df$susie_cs_coding == 0) &
                                          (combined_df$promoter_prox == FALSE),]

#Sort the dataframes
combined_no_exon <- bedr.sort.region(combined_no_exon)
combined_no_exon_cs <- bedr.sort.region(combined_no_exon_cs)
combined_no_exon_cs_prox <- bedr.sort.region(combined_no_exon_cs_prox)
#Merge the dataframes
combined_no_exon_merge <- bedr.merge.region(combined_no_exon, distance = 1000)
combined_no_exon_cs_merge <- bedr.merge.region(combined_no_exon_cs, distance = 1000)
combined_no_exon_cs_prox_merge <- bedr.merge.region(combined_no_exon_cs_prox, distance = 1000)

#Append back on the relevant info from the snp-level dataframe
append_snp_level_info <- function(merged_df_row, combined_df){
  #Split the names column
  rsids <- unlist(str_split(merged_df_row[4], ","))
  #pull the corresponding rows from the snp-level df
  temp <- combined_df[rsids,]
  temp$finemap_cs_size <- as.numeric(temp$finemap_cs_size)
  temp$susie_cs_size <- as.numeric(temp$susie_cs_size)
  method_type <- ifelse(sum(temp$method == "BOTH") > 0, "BOTH", 
                      ifelse(sum(temp$method == "SUSIE") > 0 && sum(temp$method == "FINEMAP") > 0, "BOTH_Separate",
                             temp$method[1]))
  peak_type <- ifelse(sum(temp$peak_type == "BOTH") > 0, "BOTH", 
                      ifelse(sum(temp$peak_type == "chip") > 0 && sum(temp$peak_type == "atac") > 0, "BOTH_Separate",
                             temp$peak_type[1]))
  mean_posterior_beta <- mean(as.numeric(temp$mean_beta_posterior))
  pip_indep <- 1 - prod(1 - as.numeric(temp$mean_pip))
  mean_cs_size <- mean(rowMeans(temp[,c("finemap_cs_size", "susie_cs_size")], na.rm = TRUE), na.rm = TRUE)
  region <- temp$region[1]
  num_snp <- nrow(temp)
  temp2 <- unlist(c(merged_df_row, num_snp, region, method_type, mean_posterior_beta, pip_indep, mean_cs_size, peak_type))
  names(temp2) <- NULL
  return(temp2)
}
target_df_no_exon <- as.data.frame(matrix(unlist(pbapply(combined_no_exon_merge, append_snp_level_info, MARGIN = 1, combined_df = combined_df)), byrow = TRUE, ncol = 11))
target_df_no_exon_cs <- as.data.frame(matrix(unlist(pbapply(combined_no_exon_cs_merge, append_snp_level_info, MARGIN = 1, combined_df = combined_df)), byrow = TRUE, ncol = 11))
target_df_no_exon_cs_prox <- as.data.frame(matrix(unlist(pbapply(combined_no_exon_cs_prox_merge, append_snp_level_info, MARGIN = 1, combined_df = combined_df)), byrow = TRUE, ncol = 11))
#Fix column names
colnames(target_df_no_exon) <- c("chr", "start", "end", "rsid", "num_snp", "locus", "method", "mean_posterior_beta", "indep_pip", "mean_cs_size", "peak_type")
colnames(target_df_no_exon_cs) <- c("chr", "start", "end", "rsid", "num_snp", "locus", "method", "mean_posterior_beta", "indep_pip", "mean_cs_size", "peak_type")
colnames(target_df_no_exon_cs_prox) <- c("chr", "start", "end", "rsid", "num_snp", "locus", "method", "mean_posterior_beta", "indep_pip", "mean_cs_size", "peak_type")
#Set column types
set_target_column_types <- function(target_df){
  target_df$start <- as.numeric(target_df$start)
  target_df$end <- as.numeric(target_df$end)
  target_df$num_snp <- as.numeric(target_df$num_snp)
  target_df$mean_posterior_beta <- as.numeric(target_df$mean_posterior_beta)
  target_df$indep_pip <- as.numeric(target_df$indep_pip)
  return(target_df)
}
target_df_no_exon <- set_target_column_types(target_df_no_exon)
target_df_no_exon_cs <- set_target_column_types(target_df_no_exon_cs)
target_df_no_exon_cs_prox <- set_target_column_types(target_df_no_exon_cs_prox)

#Append on a column for number of gRNAs per target site
target_df_no_exon <- cbind.data.frame(target_df_no_exon, num_grnas=ceiling((target_df_no_exon$end - target_df_no_exon$start)/repression_radius))
target_df_no_exon_cs <- cbind.data.frame(target_df_no_exon_cs, num_grnas=ceiling((target_df_no_exon_cs$end - target_df_no_exon_cs$start)/repression_radius))
target_df_no_exon_cs_prox <- cbind.data.frame(target_df_no_exon_cs_prox, num_grnas=ceiling((target_df_no_exon_cs_prox$end - target_df_no_exon_cs_prox$start)/repression_radius))

# 5) Identify Targets near Expressed genes ====

#Extract gene info from the gencode annotation
temp <- str_locate(string = gene_ann_gene$attribute, "gene_id")[,"start"]
temp2 <- str_locate(string = gene_ann_gene$attribute, "gene_name")[,"start"]
ensg_pos_df <- cbind.data.frame(chr=gene_ann_gene[,c("chr")], 
                                start = gene_ann_gene$start,
                                end = gene_ann_gene$end,
                                id=substr(gene_ann_gene$attribute, temp + 8, 
                                            str_locate(substr(gene_ann_gene$attribute, temp, nchar(gene_ann_gene$attribute)), ";") + temp - 2),
                                name=substr(gene_ann_gene$attribute, temp2 + 10, str_locate(substr(gene_ann_gene$attribute, temp2, nchar(gene_ann_gene$attribute)), ";") + temp2 - 2))
#Calculate average tpm in the differentiated hFOB replicates
hFOBsDiff_tpm_raw <- rna_raw[which(substr(rna_raw$sample, 1, 8) == "hFOBdiff"),]
hFOBsDiff_average_tpm <- aggregate(hFOBsDiff_tpm_raw$tpm, by = list(hFOBsDiff_tpm_raw$gene_id), FUN=mean)
colnames(hFOBsDiff_average_tpm) <- c("gene_id", "tpm")
rownames(hFOBsDiff_average_tpm) <- hFOBsDiff_average_tpm$gene_id
#Merge up the tpm to the gene positions
ensg_pos_df <- cbind.data.frame(ensg_pos_df, 
                                tpm=ifelse(ensg_pos_df$id %in% rownames(hFOBsDiff_average_tpm), 
                                           hFOBsDiff_average_tpm[ensg_pos_df$id,"tpm"],
                                           hFOBsDiff_average_tpm[substr(ensg_pos_df$id, 1, str_locate(pattern = "_", string = ensg_pos_df$name)-1),"tpm"]))
ensg_pos_df_express <- ensg_pos_df[which(!(is.na(ensg_pos_df$tpm))),]
ensg_pos_df_express <- ensg_pos_df_express[which(ensg_pos_df_express$tpm >= expression_thresh),]

#Pad dimensions on the expressed genes and then sort and merge to prepare for bed intersection
fix_bed_order <- function(pilot_guides_row, repression_radius){
  temp <- pilot_guides_row
  temp[2] <- as.numeric(temp[2])
  temp[3] <- as.numeric(temp[3])
  if (temp[2] > temp[3]) {
    temp[2] <- as.numeric(pilot_guides_row[3])
    temp[3] <- as.numeric(pilot_guides_row[2])
  }
  temp[2] <- max(as.numeric(temp[2]) - repression_radius, 0)
  temp[3] <- as.numeric(temp[3]) + repression_radius
  return(temp)
}
ensg_pos_df_order <- as.data.frame(matrix(unlist(pbapply(ensg_pos_df_express, MARGIN = 1, FUN = fix_bed_order, repression_radius = express_gene_range)), ncol = 6, byrow = TRUE))
colnames(ensg_pos_df_order) <- c("chr", "start", "end", "id", "name", "tpm")
ensg_pos_df_order$start <- as.numeric(ensg_pos_df_order$start) 
ensg_pos_df_order$end <- as.numeric(ensg_pos_df_order$end) 
ensg_pos_df_order <- bedr.sort.region(ensg_pos_df_order, method = "lexicographical") #Sort df
ensg_pos_df_merge <- bedr.merge.region(ensg_pos_df_order) #merge df

#Intersect the expressed genes with the target-level dataframes
intersect_w_express_genes <- function(target_df, ensg_pos_df_merge) {
  target_df_express <- bedr(
    input = list(a = target_df, b = ensg_pos_df_merge), 
    method = "intersect", 
    params = "-c -sorted"
  )
  colnames(target_df_express) <- c(colnames(target_df), "near_expressed_gene")
  target_df_express$near_expressed_gene <- ifelse(target_df_express$near_expressed_gene > 0, TRUE, FALSE)
  return(target_df_express)
}
target_df_no_exon <- intersect_w_express_genes(target_df_no_exon, ensg_pos_df_merge = ensg_pos_df_merge)
target_df_no_exon_cs <- intersect_w_express_genes(target_df_no_exon_cs, ensg_pos_df_merge = ensg_pos_df_merge)
target_df_no_exon_cs_prox <- intersect_w_express_genes(target_df_no_exon_cs_prox, ensg_pos_df_merge = ensg_pos_df_merge)

# 6) Identify the signals that were targeted previously ====

#Fix start and end order if misordered and pad ranges of pilot guides
pilot_guides_clean <- as.data.frame(matrix(unlist(pbapply(pilot_guides_raw, FUN = fix_bed_order, MARGIN = 1, repression_radius=repression_radius)), ncol = 4, byrow = TRUE))
pilot_guides_clean[,2] <- as.numeric(pilot_guides_clean[,2])
pilot_guides_clean[,3] <- as.numeric(pilot_guides_clean[,3])
#Merge and sort the pilot guides
pilot_guides_clean <- bedr.sort.region(pilot_guides_clean, method = "lexicographical")
pilot_guides_merge <- bedr.merge.region(pilot_guides_clean)

#Intersect the pilot guides with the new targets
check_pilot_guides <- function(target_df, pilot_guides_merge){
  target_df_intersect <- bedr(
    input = list(a = target_df, b = pilot_guides_merge), 
    method = "intersect", 
    params = "-c -sorted"
  )
  colnames(target_df_intersect) <- c(colnames(target_df), "pilot_inclusion")
  target_df_intersect$pilot_inclusion <- ifelse(target_df_intersect$pilot_inclusion > 0, TRUE, FALSE)
  return(target_df_intersect)
}
target_df_no_exon <- check_pilot_guides(target_df_no_exon, pilot_guides_merge = pilot_guides_merge)
target_df_no_exon_cs <- check_pilot_guides(target_df_no_exon_cs, pilot_guides_merge = pilot_guides_merge)
target_df_no_exon_cs_prox <- check_pilot_guides(target_df_no_exon_cs_prox, pilot_guides_merge = pilot_guides_merge)

# 7) Make plots of the number of Targets/gRNAs Under Different Thresholds ====

#Make a histogram making function
make_histogram <- function(target_df, plot_dir, file_prefix, expressed = TRUE){
  
  #Remove non-expressed genes if flag set
  if (expressed == TRUE) {
    target_df <- target_df[target_df$near_expressed_gene == TRUE,]
  }
  
  ## Make hist of number of targets per region/locus ##
  num_loci <- length(unique(target_df$locus))
  num_targets <- nrow(target_df)
  temp <- table(target_df[,c("locus")])
  locus_df <- cbind.data.frame(locus= names(temp), temp)
  plot_ul <- ggplot(locus_df, aes(x=Freq)) + geom_histogram(binwidth = 1) + 
    xlab("Targets per Locus") + ylab("Number of Loci") +
    annotate("text", x=Inf, y=Inf, label= paste0(num_targets, " targets for ", num_loci, " loci"), hjust=1, vjust=1, size=5) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"),  
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
  
  ## Make hist of number of gRNAs per target ##
  target_df$num_grnas <- as.numeric(target_df$num_grnas)
  num_grnas <- sum(target_df$num_grnas)
  plot_ll <- ggplot(target_df, aes(x=num_grnas, fill=pilot_inclusion)) + geom_histogram(binwidth = 1) + 
    xlab("gRNAs per Target Region") + ylab("Number of Targets") +
    stat_bin(aes(x=num_grnas, y=after_stat(count), label=after_stat(count)), geom="text", vjust=-.5, binwidth = 1, inherit.aes = FALSE)  + 
    scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                        breaks=c(TRUE, FALSE),
                        labels=c("Pilot Target", "Novel Target")) + 
    annotate("text", x=Inf, y=Inf, label= paste0(num_grnas, " gRNAs for\n", num_targets, " targets"), hjust=1, vjust=1, size=5) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 14),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
  
  ## Make hist of PIPs per target ##
  #Make plot
  target_df$indep_pip <- as.numeric(target_df$indep_pip)
  plot_lm <- ggplot(target_df, aes(x=indep_pip, fill=pilot_inclusion)) + geom_histogram(binwidth = 0.1) + 
    xlab("Combined Mean PIP") + ylab("Number of Targets") +
    stat_bin(aes(x=indep_pip, y=after_stat(count), label=after_stat(count)), geom="text", vjust=-.5, binwidth = 0.1, inherit.aes = FALSE)  + 
    scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                      breaks=c(TRUE, FALSE),
                      labels=c("Pilot Target", "Novel Target")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 14),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
  
  ### Make heatmap of targets for pip vs credible set size ###
  #Bin the pips and the cs sizes
  target_df$mean_cs_size <- as.numeric(target_df$mean_cs_size)
  target_df <- cbind.data.frame(target_df, 
                                pip_bin=as.factor(floor(target_df$indep_pip/0.1)*0.1),
                                cs_size=as.factor(ifelse(target_df$mean_cs_size == 1, "1", 
                                               ifelse(target_df$mean_cs_size <= 5, "2-5",
                                                      ifelse(target_df$mean_cs_size <= 10, "6-10",
                                                             ifelse(target_df$mean_cs_size <= 25, "11-25",
                                                                    ifelse(target_df$mean_cs_size <= 50, "26-50","51+")))))))
  plot_um <- ggplot(target_df, aes(x=indep_pip, y=mean_cs_size)) + geom_bin_2d(bins = 20) + 
    xlab("Combined Mean PIP") + ylab("Mean CS Size") + xlim(-0.01,1.01) +
    scale_y_continuous(breaks = seq(0,100,10)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 12),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))

  ### Make heatmap of targets for peak type vs method type ###
  temp <- melt(table(target_df[,c("peak_type", "method")]))
  plot_ur <- ggplot(temp, aes(x=peak_type, y=method)) + geom_tile(aes(fill = value)) + 
    geom_text(aes(label = round(value, 1)), size = 10) +
    xlab("Peak Type") + ylab("Mapping Method") + 
    scale_x_discrete(breaks= c("atac", "chip", "BOTH_Separate", "BOTH"),
                     labels = c("ATAC", "ChIP", "Both\nSeparate", "Both")) + 
    scale_y_discrete(breaks= c("SUSIE", "FINEMAP", "BOTH_Separate", "BOTH"),
                     labels = c("SuSiE", "FINEMAP", "Both\nSeparate", "Both")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 12),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
  
  ## Make hist of PIPs per target colored by peak types ##
  #Make plot
  target_df$indep_pip <- as.numeric(target_df$indep_pip)
  plot_lr <- ggplot(target_df, aes(x=indep_pip, fill=peak_type)) + geom_histogram(binwidth = 0.1) + 
    xlab("Combined Mean PIP") + ylab("Number of Targets") +
    annotate("text", x=Inf, y=Inf, label= paste0(nrow(target_df[which(target_df$peak_type == "BOTH" & target_df$method == "BOTH" & target_df$pilot_inclusion == FALSE),]), 
                                                 " novel and ", nrow(target_df[which(target_df$peak_type == "BOTH" & target_df$method == "BOTH" & target_df$pilot_inclusion == TRUE),]), 
                                                 " pilot\ntargets found in both\npeaks by both methods"), hjust=1, vjust=1, size=5) + 
    stat_bin(aes(x=indep_pip, y=after_stat(count), label=after_stat(count)), geom="text", vjust=-.5, binwidth = 0.1, inherit.aes = FALSE)  + 
    scale_fill_manual(values=rainbow(4), 
                      breaks= c("atac", "chip", "BOTH_Separate", "BOTH"),
                      labels = c("ATAC", "ChIP", "Both\nSeparate", "Both")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 14),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14))
  
  #Stitch plots together
  jpeg(file=paste0(plot_dir, file_prefix, ".combined.jpg"), width = 18, height = 10, units = 'in', res = 200)
  print(plot_grid(plot_ul, plot_um, plot_ur, plot_ll, plot_lm, plot_lr, ncol = 3), scale = 0.9)
  dev.off()
}

#Call function
make_histogram(target_df_no_exon_cs_prox, plot_dir, "no_exon_cs_prox")
make_histogram(target_df_no_exon_cs_prox[which(target_df_no_exon_cs_prox$num_grnas == 1),], plot_dir, "no_exon_cs_prox.1_grna")
make_histogram(target_df_no_exon_cs_prox[which(target_df_no_exon_cs_prox$indep_pip > 0.2),], plot_dir, "no_exon_cs_prox.pip_gt_0.2")
make_histogram(target_df_no_exon_cs_prox[which(target_df_no_exon_cs_prox$num_grnas == 1 & target_df_no_exon_cs_prox$indep_pip > 0.2),], plot_dir, "no_exon_cs_prox.1_grna.pip_gt_0.2")

#Write file
write.table(target_df_no_exon_cs_prox, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 8) Test Enrichment of Both Peak Types in Higher PIP Targets ====

#Make a function that I can call to test enrichment of both peak types in higher pip variants
make_both_peak_enrichment_bar <- function(jump, file_prefix, target_df){
  #Convert to independent pip to numeric
  target_df$indep_pip <- as.numeric(target_df$indep_pip)
  #Create buckets for PIPs
  pip_bucket_df <- cbind.data.frame(pip_bucket_mins=seq(0, 1-jump, jump), pip_bucket_maxes=seq(jump, 1, jump))
  round_selective <- function(inp_num, rdigits=2){ifelse(inp_num %% 1 == 0, trimws(format(inp_num, nsmall = 0)), trimws(format(round(inp_num,digits = rdigits), nsmall=rdigits)))}
  pip_bucket_df <- cbind.data.frame(pip_bucket_df, pip_bin=paste0(vapply(pip_bucket_df$pip_bucket_mins, FUN = round_selective, FUN.VALUE = character(1)), "-", 
                                                                  vapply(pip_bucket_df$pip_bucket_maxes, FUN = round_selective, FUN.VALUE = character(1))))
  #Calculate percentage of targets with both peak types in pip buckets
  calc_both_peak_pct <- function(pip_bucket_row, target_df){
    target_df$indep_pip <- as.numeric(target_df$indep_pip)
    temp <- target_df[which(target_df$indep_pip > pip_bucket_row["pip_bucket_mins"] & target_df$indep_pip <= pip_bucket_row["pip_bucket_maxes"]),]
    both_count <- nrow(temp[which(temp$peak_type == "BOTH"),])
    other_count <- nrow(temp) - both_count
    return(c(both_count, other_count, both_count/nrow(temp)))
  }
  temp <- t(apply(pip_bucket_df, MARGIN = 1, FUN = calc_both_peak_pct, target_df=target_df))
  pip_bucket_df <- cbind.data.frame(pip_bucket_df, temp)
  colnames(pip_bucket_df) <- c("pip_bucket_mins", "pip_bucket_maxes", "pip_bin", "Both_Peak_Targets", "Other_Targets", "Both_Peak_PCT")
  
  #Run statistical tests of whether enrichment is significant
  test_peak_enrichment <- function(big_row, pip_bucket_df){
    temp <- fisher.test(pip_bucket_df[c(big_row-1,big_row), c("Both_Peak_Targets", "Other_Targets")])
    return(t(c(pip_bucket_df[c(big_row-1,big_row),"pip_bin"],format(round(temp$p.value, digits = 2), digits = 2))))
  }
  pval_df <- as.data.frame(rbind.fill.matrix(pblapply(c(2:nrow(pip_bucket_df)), FUN = test_peak_enrichment, pip_bucket_df=pip_bucket_df)))
  pval_df <- cbind.data.frame(rep(1.05, nrow(pval_df)), pval_df)
  colnames(pval_df) <- c("y.position", "group1", "group2", "p")
  
  #Bin the data out in the target_df
  take_min_with <- function(inp_num, min_num){return(min(inp_num, min_num))}
  target_df_append <- cbind.data.frame(target_df, 
                                                       pip_bin=paste0(vapply(vapply((target_df$indep_pip%/%jump)*jump, FUN = take_min_with, FUN.VALUE = numeric(1), min_num=(1%/%jump - 1)*jump), FUN = round_selective, FUN.VALUE = character(1)), "-", 
                                                                      vapply(vapply((target_df$indep_pip%/%jump)*jump + jump, FUN = take_min_with, FUN.VALUE = numeric(1), min_num=1), FUN = round_selective, FUN.VALUE = character(1))))
  target_df_append$pip_bin <- as.factor(target_df_append$pip_bin) 
  target_df_append$peak_type <- factor(target_df_append$peak_type, levels = c("atac", "chip", "BOTH_Separate", "BOTH"))
  
  #Make histogram
  target_df_append$indep_pip <- as.numeric(target_df_append$indep_pip)
  jpeg(file=paste0(plot_dir, file_prefix, ".peak_type.hist.jpg"), res = 200, width = 3, height = 3, units = "in")
  print(ggplot(target_df_append, aes(x=pip_bin, fill=peak_type)) + 
    geom_bar(position = "fill") + 
    #stat_pvalue_manual(pval_df, label = "p = {p}", bracket.shorten = -0.1, label.size = 6) + 
    xlab("Combined Mean PIP") + ylab("Prop. of Targets in Peak Type") +
    scale_fill_manual(values=rainbow(4), 
                      breaks= c("atac", "chip", "BOTH_Separate", "BOTH"),
                      labels = c("ATAC", "ChIP", "Both\nSeparate", "Both")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.background = element_rect(fill = "transparent",colour = NA), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.text =  element_text(size = 5),
          plot.background = element_rect(fill = "transparent",colour = NA), axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5)))
  dev.off()
}

#Call Function
make_both_peak_enrichment_bar(1/3, "tercile", target_df = target_df_no_exon_cs_prox)
make_both_peak_enrichment_bar(1/4, "quartile", target_df = target_df_no_exon_cs_prox)
make_both_peak_enrichment_bar(1/5, "quintile", target_df = target_df_no_exon_cs_prox)

