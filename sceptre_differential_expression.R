################################################################################

#sceptre_differential_expression.R

#The purpose of this script is to test the cells containing a given sgRNA for  
#differential expression of nearby genes using the SCEPTRE methodology. 
#The script is set to run with the following modules:

#  R/4.2.3

#The script will try two different methods for detecting differential 
#expression of genes:

# 1) Compare a sgRNA to all cells without the sgRNA and only receiving 1 sgRNA
# 2) Compare a sgRNA to all cells without the sgRNA

################################################################################

# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(dplyr)
library(cowplot)
library(katlabutils)
library(sceptre)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(stringr)
library(future)
library(plyr)
library(pbapply)
library(bedr)
library(SingleCellExperiment)
library(Matrix)
library(ggplot2)
library(tibble)

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
gencode_genes_loc <- "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/gencode_v19/gencode.v19.annotation.gene_only.bed"
sgrna_target_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/all_targeting_guides.bed"
marker_genes_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/quality_control/hFOB_marker_genes.csv"
implicate_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/gene_implications.tsv"
vep_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/target.vep_annotations.txt"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/differential_expression/sceptre_single_sgrna/"

#Set repression range
rep_range <- 1000000 #Testing within +/- 1Mb
#Set random seed
set.seed(5)

# 1) Read in files ====

#Read in Scanpy QC-ed cell results
Convert(qced_results_loc, dest = "h5seurat", overwrite = TRUE)
qced_results_meta <- read.csv(str_replace(qced_results_loc, "[.]h5ad", "/obs.csv"))
rownames(qced_results_meta) <- qced_results_meta$X
qced_results_meta <- qced_results_meta[,c(2:ncol(qced_results_meta))]
qced_results_seurat <- LoadH5Seurat(str_replace(qced_results_loc, "h5ad", "h5seurat"), meta.data = FALSE, misc = FALSE)
qced_results_seurat <- AddMetaData(object = qced_results_seurat, metadata = qced_results_meta)

#Read in number of guides per cell
guides_per_cell_raw <- read.csv(guides_per_cell_loc, header = TRUE, sep = ",")

#Read in target locations file 
sgrna_target_raw <- read.table(sgrna_target_loc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sgrna_target_raw) <- c("chr", "start", "end", "name+seq")

#Read in genes file
gencode_genes_raw <- read.table(gencode_genes_loc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gencode_genes_raw) <- c("chr", "start", "end", "ENSG+name", "gene_type", "strand")

#Read in marker genes file
marker_genes_raw <- read.table(marker_genes_loc, sep = ",", header = TRUE)

#Read in list of implicated genes
implicate_raw <- read.table(implicate_loc, header = TRUE, stringsAsFactors = FALSE, sep = )

#Read in vep annotations file
vep_raw <- read.csv(vep_loc, sep = "\t", header = TRUE)

# 2) Get List of Genes to Test within Range of sgRNA Targets ====

#Get list of unique sgrna targets
sgrnas <- unique(guides_per_cell_raw$feature_call)
sgrnas <- sgrnas[str_detect(sgrnas, "[|]") == FALSE]
#Divy up the guide rnas into their respective categories
pos_control_sgrnas <- sgrnas[1:2]
neg_control_sgrnas <- sgrnas[3:29]
test_sgrnas <- sgrnas[30:length(sgrnas)]

#Split up the final column of the sgrna target file
split_plus_take_one <- function(string){str_split(string, "[+]")[[1]][1]}
split_plus_take_two <- function(string){str_split(string, "[+]")[[1]][2]}
sgrna_target_append <- cbind.data.frame(sgrna_target_raw, name=vapply(sgrna_target_raw[,"name+seq"], FUN = split_plus_take_one, FUN.VALUE = character(1)))
sgrna_target_append$name <- str_replace_all(sgrna_target_append$name, ",", "-")
rownames(sgrna_target_append) <- sgrna_target_append$name
sgrna_target_append <- sgrna_target_append[order(sgrna_target_append$chr, sgrna_target_append$start),]
#Add on columns for where the repression range should start and end
take_max_w_zero <- function(value){return(max(value,0))}
sgrna_target_append <- cbind.data.frame(sgrna_target_append, 
                                        rep_start=vapply(sgrna_target_append$start - rep_range, FUN = take_max_w_zero, FUN.VALUE = 0),
                                        rep_end=(sgrna_target_append$end + rep_range))

#Drop the sgRNAs that dropped out (namely rs113214136_3)
sgrna_target_filter <- sgrna_target_append[which(sgrna_target_append$name %in% sgrnas),]
#Move the repression ranges to the start and end positions of the sgrna_target_filter df
sgrna_target_filter[,"bind_start"] <- sgrna_target_filter[,"start"]
sgrna_target_filter[,"bind_end"] <- sgrna_target_filter[,"end"]
sgrna_target_filter$start <- sgrna_target_filter$rep_start
sgrna_target_filter$end <- sgrna_target_filter$rep_end

#Split up name+seq column in gencode annotations
gencode_genes_append <- cbind.data.frame(gencode_genes_raw, 
                                         ensembl_gene_id=vapply(gencode_genes_raw[,"ENSG+name"], FUN = split_plus_take_one, FUN.VALUE = character(1)),
                                         external_gene_name=vapply(gencode_genes_raw[,"ENSG+name"], FUN = split_plus_take_two, FUN.VALUE = character(1)))
rownames(gencode_genes_append) <- gencode_genes_append$ensembl_gene_id

#Filter the gencode_genes down for genes that are actually measured in the screen
gencode_genes_filter <- gencode_genes_append[which(str_replace(gencode_genes_append$ensembl_gene_id, "[.][0-9]+", "") %in% rownames((qced_results_seurat))),]

#Reorder the bed file style df columns
sgrna_target_filter <- sgrna_target_filter[,c("chr","start","end","name","rep_start","rep_end","bind_start","bind_end")]
colnames(sgrna_target_filter) <- c("chr","start","end","sgrna","rep_start","rep_end","bind_start","bind_end")
gencode_genes_filter <- gencode_genes_filter[,c("chr","start","end","ensembl_gene_id","external_gene_name","gene_type","strand")]
#Merge bed files
sgrna_target_merge <- bedr.merge.region(sgrna_target_filter)
gencode_genes_merge <- bedr.merge.region(gencode_genes_filter)
#Sort bed files
gencode_genes_merge <- bedr.sort.region(gencode_genes_merge)
sgrna_target_merge <- bedr.sort.region(sgrna_target_merge)
#Intersect the gencode and the sgrna target files bedtools style
gencode_genes_intersect <- bedr(
  input = list(a = sgrna_target_merge, b = gencode_genes_merge), 
  method = "intersect", 
  params = "-wb -sorted"
)
#Get ensembl ids of retained genes
test_genes <- unique(unlist(str_split(str_replace_all(gencode_genes_intersect[,8], "[.][0-9]+", ""),",")))

#Reset names on gencode list
rownames(gencode_genes_filter) <- str_replace(gencode_genes_filter$ensembl_gene_id, "[.][0-9]+", "")

#Create two lists for converting gene symbols to ensg and vice versa
symbol_to_ensg <- rownames(gencode_genes_filter)
names(symbol_to_ensg) <- gencode_genes_filter$external_gene_name
ensg_to_symbol <- gencode_genes_filter$external_gene_name
names(ensg_to_symbol) <- rownames(gencode_genes_filter)

# 3) Prepare Raw Inputs to SCEPTRE Test ====

### Make gRNA-by-cell matrix ###
#Get list of cells in qced results
qced_cells <- colnames(qced_results_seurat)
#Make a function to make a vector of sgrna counts for each cell
get_sgrnas_per_cell <- function(qced_cell, sgrnas, guides_per_cell_raw){
  temp <- guides_per_cell_raw[guides_per_cell_raw$cell_barcode == qced_cell,]
  cell_sgrnas <- unlist(str_split(temp$feature_call, "[|]"))
  num_umis <- as.numeric(unlist(str_split(temp$num_umis, "[|]")))
  names(num_umis) <- cell_sgrnas
  out <- num_umis[sgrnas]
  out <- ifelse(is.na(out), 0, out)
  names(out) <- NULL
  return(t(out))
}
#Call function to make sgrna count matrix
grna_list <- pblapply(qced_cells, get_sgrnas_per_cell, sgrnas=sgrnas, guides_per_cell_raw=guides_per_cell_raw)
grna_hfob_lowmoi <- as.matrix(t(rbind.fill.matrix(grna_list)))
rownames(grna_hfob_lowmoi) <- sgrnas
#Filter for cells with at least 1 guide
cells_with_guide <- qced_cells[colSums(grna_hfob_lowmoi != 0) != 0]
grna_hfob_lowmoi <- grna_hfob_lowmoi[,colSums(grna_hfob_lowmoi != 0) != 0]

### Make response-by-cell matrix ###
sce <- as.SingleCellExperiment(qced_results_seurat)
counts_matrix <- as.matrix(assay(sce, "counts"))
response_hfob_lowmoi <- counts_matrix[,cells_with_guide]
colnames(response_hfob_lowmoi) <- NULL

### Make Cell Covariate Matrix ###
covariate_hfob_lowmoi <- qced_results_seurat@meta.data
covariate_hfob_lowmoi <- covariate_hfob_lowmoi[cells_with_guide,c("total_counts", "n_genes_by_counts", "Pool", "pct_counts_mt")]
rownames(covariate_hfob_lowmoi) <- NULL
colnames(covariate_hfob_lowmoi) <- c("response_n_umis", "response_n_nonzero", "bio_rep", "p_mito")
covariate_hfob_lowmoi$bio_rep <- as.factor(covariate_hfob_lowmoi$bio_rep)
covariate_hfob_lowmoi$p_mito <- covariate_hfob_lowmoi$p_mito/100 #rescale percentages to decimal
#Calculate PCAs on Seurat Object and add top 15 to covariate matrix
qced_results_seurat <- FindVariableFeatures(qced_results_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(qced_results_seurat)
qced_results_seurat <- ScaleData(qced_results_seurat, features = all.genes)
qced_results_seurat <- RunPCA(qced_results_seurat, rev.pca=TRUE)
qced_results_seurat <- ProjectDim(qced_results_seurat, reduction="pca")
temp <- Embeddings(qced_results_seurat, reduction="pca")[cells_with_guide,]
covariate_hfob_lowmoi <- cbind.data.frame(covariate_hfob_lowmoi, temp[,1:15])

### Make gRNA group table ###
grna_group_hfob_lowmoi <- rbind.data.frame(cbind.data.frame(grna_id=test_sgrnas, grna_group=str_replace(test_sgrnas, "_[0-9]+","")),
                                               cbind.data.frame(grna_id=pos_control_sgrnas, grna_group=pos_control_sgrnas),
                                               cbind.data.frame(grna_id=neg_control_sgrnas, grna_group=rep("non-targeting", length(neg_control_sgrnas))))

# 4) Check for Random Assortment of Non-Targeting and Targeting Guides ====

### Run a test for non-random assortment of non-targeting guides and targeting guides ###
check_ntc_assortment <- function(neg_control_sgrna, test_sgrnas, grna_hfob_lowmoi){
  temp <- vapply(test_sgrnas, check_combo_assortment, FUN.VALUE = 0, neg_control_sgrna=neg_control_sgrna, grna_hfob_lowmoi=grna_hfob_lowmoi)
  return(cbind.data.frame(ntc=rep(neg_control_sgrna, length(temp)), test_sgrna=test_sgrnas, fisher_pval=temp))
}
check_combo_assortment <- function(test_sgrna, neg_control_sgrna, grna_hfob_lowmoi){
  temp <- ifelse(grna_hfob_lowmoi[c(test_sgrna, neg_control_sgrna),]==0,0,1)
  check_matrix <- matrix(data = c(sum(colSums(temp) == 2), sum(colSums(temp) == 1 & temp[1,]==1),
                                  sum(colSums(temp) == 1 & temp[2,]==1), sum(colSums(temp) == 0)), nrow = 2, ncol = 2)
  temp2 <- fisher.test(check_matrix)
  return(temp2$p.value)
}
temp <- lapply(neg_control_sgrnas, check_ntc_assortment, grna_hfob_lowmoi=grna_hfob_lowmoi, test_sgrnas=test_sgrnas)
non_random_assort_df <- as.data.frame(rbind.fill.matrix(temp))
#Do a multiple testing correction
non_random_assort_df <- cbind.data.frame(non_random_assort_df, pval_adj=p.adjust(as.numeric(non_random_assort_df$fisher_pval), method = "BH"))
non_random_assort_df <- non_random_assort_df[order(non_random_assort_df$pval_adj),]

#Get list of targeting grna groups
test_grna_groups <- unique(grna_group_hfob_lowmoi[which(grna_group_hfob_lowmoi$grna_group != "non-targeting"),"grna_group"])
#Run a test for non-random assortment of non-targeting guides with groups of targeting guides 
check_ntc_group_assortment <- function(neg_control_sgrna, test_grna_groups, grna_hfob_lowmoi, grna_group_hfob_lowmoi){
  temp <- vapply(test_grna_groups, check_combo_group_assortment, FUN.VALUE = 0, neg_control_sgrna=neg_control_sgrna, grna_hfob_lowmoi=grna_hfob_lowmoi, grna_group_hfob_lowmoi=grna_group_hfob_lowmoi)
  return(cbind.data.frame(ntc=rep(neg_control_sgrna, length(temp)), test_sgrna_group=test_grna_groups, fisher_pval=temp))
}
check_combo_group_assortment <- function(test_grna_group, neg_control_sgrna, grna_hfob_lowmoi, grna_group_hfob_lowmoi){
  test_grnas <- grna_group_hfob_lowmoi[which(grna_group_hfob_lowmoi$grna_group == test_grna_group),"grna_id"]
  temp <- ifelse(grna_hfob_lowmoi[c(test_grnas, neg_control_sgrna),]==0,0,1)
  temp_rows <- nrow(temp)
  temp_top <- temp[1:(temp_rows-1),]
  temp_bottom <- temp[temp_rows,]
  if (temp_rows == 2) {
    top_left <- sum(temp_top * temp_bottom >= 1)
    top_right <- sum(temp_top) - top_left
    bottom_left <- sum(temp_bottom == 1) - top_left
    bottom_right <- length(temp_bottom) - top_left - top_right - bottom_left
  } else {
    top_left <- sum(colSums(temp_top) * temp_bottom >= 1)
    top_right <- sum(colSums(temp_top) >= 1) - top_left
    bottom_left <- sum(temp_bottom == 1) - top_left
    bottom_right <- length(temp_bottom) - top_left - top_right - bottom_left
  }
  check_matrix <- matrix(data = c(top_left, top_right, bottom_left, bottom_right), 
                         nrow = 2, ncol = 2)
  temp2 <- fisher.test(check_matrix)
  return(temp2$p.value)
}
temp <- pblapply(neg_control_sgrnas, check_ntc_group_assortment, grna_hfob_lowmoi=grna_hfob_lowmoi, test_grna_groups=test_grna_groups, grna_group_hfob_lowmoi=grna_group_hfob_lowmoi)
non_random_assort_group_df <- as.data.frame(rbind.fill.matrix(temp))
#Do a multiple testing correction
non_random_assort_group_df <- cbind.data.frame(non_random_assort_group_df, pval_adj=p.adjust(as.numeric(non_random_assort_group_df$fisher_pval), method = "BH"))
non_random_assort_group_df <- non_random_assort_group_df[order(non_random_assort_group_df$pval_adj),]

### Get non-targeting guide RNAs with non-random assortments ###
non_random_guides <- unique(non_random_assort_group_df[which(non_random_assort_group_df$pval_adj < 0.05),"ntc"])
# 15 guides show non-random assortment which necessitates focusing on the single-guide RNA cells

# 5) Establish Pairs of Genes/sgRNAs to Test ====

#Make a map for each guide to the list of genes we are testing it against
check_range <- function(i, chromos, rep_starts, rep_ends, gencode_genes_filter){
  chromo <- chromos[i]
  rep_start <- rep_starts[i]
  rep_end <- rep_ends[i]
  gencode_genes_temp <- gencode_genes_filter[which(gencode_genes_filter$chr == chromo & 
                                                     ((gencode_genes_filter$start >= rep_start & gencode_genes_filter$start <= rep_end) | 
                                                        (gencode_genes_filter$end >= rep_start & gencode_genes_filter$end <= rep_end) | 
                                                        (gencode_genes_filter$start <= rep_start & gencode_genes_filter$end >= rep_end))),"ensembl_gene_id"]
  return(str_replace(gencode_genes_temp, "[.][0-9]+",""))
}
check_genes_for_guide <- function(grna_group, grna_group_hfob_lowmoi, sgrna_target_filter, gencode_genes_filter){
  sgrnas <- grna_group_hfob_lowmoi[which(grna_group_hfob_lowmoi$grna_group == grna_group),"grna_id"]
  chromos <- as.character(sgrna_target_filter[sgrnas,"chr"])
  rep_starts <- as.numeric(sgrna_target_filter[sgrnas,"rep_start"])
  rep_ends <- as.numeric(sgrna_target_filter[sgrnas,"rep_end"])
  temp <- unique(unlist(lapply(seq(1, length(chromos), 1), FUN = check_range, chromos=chromos, rep_starts=rep_starts, rep_ends=rep_ends,gencode_genes_filter)))
  temp2 <- cbind.data.frame(response_id=temp, grna_group=rep(grna_group, length(temp)))
  return(temp2)
}
response_grna_group_pairs_hfob <- as.data.frame(rbind.fill.matrix(pblapply(unique(str_replace(test_sgrnas, "_[0-9]+", "")), FUN = check_genes_for_guide, 
                                                                           grna_group_hfob_lowmoi=grna_group_hfob_lowmoi, sgrna_target_filter=sgrna_target_filter, 
                                                                           gencode_genes_filter=gencode_genes_filter)))
#Append on the positive controls
response_grna_group_pairs_hfob <- rbind.data.frame(response_grna_group_pairs_hfob,
                                                   cbind.data.frame(response_id=symbol_to_ensg[pos_control_sgrnas], grna_group=pos_control_sgrnas))
rownames(response_grna_group_pairs_hfob) <- NULL
#Make columns factors
#response_grna_group_pairs_hfob$response_id <- as.factor(response_grna_group_pairs_hfob$response_id)
#response_grna_group_pairs_hfob$grna_group <- as.factor(response_grna_group_pairs_hfob$grna_group)

#Do a final filter on the response matrix for genes that are actually being tested for some guide
response_hfob_lowmoi <- response_hfob_lowmoi[unique(response_grna_group_pairs_hfob$response_id),]

# 6) Filter For Single-sgRNA Cells ====

#Remove the cells with multiple sgRNAs from consideration
response_hfob_lowmoi <- response_hfob_lowmoi[,colSums(grna_hfob_lowmoi != 0) == 1]
covariate_hfob_lowmoi <- covariate_hfob_lowmoi[colSums(grna_hfob_lowmoi != 0) == 1,]
grna_hfob_lowmoi <- grna_hfob_lowmoi[,colSums(grna_hfob_lowmoi != 0) == 1]
#Remove any guides no longer represented
grna_hfob_lowmoi <- grna_hfob_lowmoi[rowSums(grna_hfob_lowmoi) > 0,]
grna_group_hfob_lowmoi <- grna_group_hfob_lowmoi[grna_group_hfob_lowmoi$grna_id %in% rownames(grna_hfob_lowmoi),]
response_grna_group_pairs_hfob <- response_grna_group_pairs_hfob[response_grna_group_pairs_hfob$grna_group %in% grna_group_hfob_lowmoi$grna_group,]

# 7) Set the Formula Object ====

#Set the formula object for use with pcas
formula_object <- formula(~log(response_n_umis) + 
                            log(response_n_nonzero) +
                            bio_rep + 
                            p_mito + 
                            PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + 
                              PC_6 + PC_7 + PC_8 + PC_9 + PC_10 +
                              PC_11 + PC_12 + PC_13 + PC_14 + PC_15)

# 8) Run and Assess Calibration Check ====

#Run basic check
calibration_result <- run_sceptre_lowmoi(
  response_matrix = response_hfob_lowmoi,
  grna_matrix = grna_hfob_lowmoi,
  covariate_data_frame = covariate_hfob_lowmoi,
  grna_group_data_frame = grna_group_hfob_lowmoi,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_hfob,
  calibration_check = TRUE # calibration_check TRUE
) 

#Plot calibration result
jpeg(paste0(out_dir, "sceptre_calibration.hfob.jpeg"))
plot_calibration_result(calibration_result)
dev.off()
#Write calibration results to file
write.table(calibration_result, file = paste0(out_dir, "non_targeting_test_results.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 9) Run Discovery Analyses ====

### Run Discovery on single-guide cells ###
discovery_result <- run_sceptre_lowmoi(
  response_matrix = response_hfob_lowmoi,
  grna_matrix = grna_hfob_lowmoi,
  covariate_data_frame = covariate_hfob_lowmoi,
  grna_group_data_frame = grna_group_hfob_lowmoi,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_hfob,
  calibration_check = FALSE
) 

#Compare calibration and discovery results
jpeg(paste0(out_dir, "sceptre_discovery_calibration_compare.hfob.jpeg"))
compare_calibration_and_discovery_results(
  calibration_result = calibration_result,
  discovery_result = discovery_result
)
dev.off()

#Make volcano
jpeg(paste0(out_dir, "discovery_volcano.hfob.jpeg"))
make_volcano_plot(discovery_result = discovery_result)
dev.off()

#obtain discovery set and write to file
discovery_set <- obtain_discovery_set(discovery_result)
discovery_set$response_id <- ensg_to_symbol[discovery_set$response_id]
write.table(discovery_set, file = paste0(out_dir, "discovery_results.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 10) Examine Results for V2G-Mapped Connections ====

### Make dataframe of grna pairs ###
#create unique identifiers for each implicated connection
join_genes_to_proxies <- function(i, implicate_raw, symbol_to_ensg){return(unlist(lapply(str_split(implicate_raw[i,"genes_any"], ","), FUN = paste_proxy, proxy_char=implicate_raw[i,"proxies"], symbol_to_ensg=symbol_to_ensg)))}
paste_proxy <- function(gene, proxy_char, symbol_to_ensg){return(ifelse(is.na(symbol_to_ensg[gene]), NA, str_replace_all(paste0(symbol_to_ensg[gene], ".", proxy_char),pattern = ",", replacement = "-")))}
implicated_connections <- unlist(lapply(1:nrow(implicate_raw), FUN = join_genes_to_proxies, implicate_raw=implicate_raw, symbol_to_ensg=symbol_to_ensg))
names(implicated_connections) <- NULL
implicated_connections <- implicated_connections[is.na(implicated_connections)==FALSE]
#Convert it to a response pairs dataframe
split_period_take_one <- function(string){str_split(string, "[.]")[[1]][1]}
split_period_take_two <- function(string){str_split(string, "[.]")[[1]][2]}
response_grna_group_pairs_implicated <- cbind.data.frame(
  response_id=vapply(implicated_connections, split_period_take_one, character(1)),
  grna_group=vapply(implicated_connections, split_period_take_two, character(1)))
rownames(response_grna_group_pairs_implicated) <- NULL
#Make response matrix for genes that were implicated
response_hfob_implicated <- counts_matrix[,cells_with_guide]
response_hfob_implicated <- response_hfob_implicated[unique(response_grna_group_pairs_implicated$response_id),rownames(covariate_hfob_lowmoi)]
colnames(response_hfob_implicated) <- NULL

### Run calibration and discovery ###
calibration_result_implicated <- run_sceptre_lowmoi(
  response_matrix = response_hfob_implicated,
  grna_matrix = grna_hfob_lowmoi,
  covariate_data_frame = covariate_hfob_lowmoi,
  grna_group_data_frame = grna_group_hfob_lowmoi,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_implicated,
  calibration_check = TRUE # calibration_check TRUE
) 
discovery_result_implicated <- run_sceptre_lowmoi(
  response_matrix = response_hfob_implicated,
  grna_matrix = grna_hfob_lowmoi,
  covariate_data_frame = covariate_hfob_lowmoi,
  grna_group_data_frame = grna_group_hfob_lowmoi,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_implicated,
  calibration_check = FALSE
) 

### obtain discovery set and write to file ###
#Get vep annotations for targets
proxy_to_vep <- vep_raw$Consequence
names(proxy_to_vep) <- vep_raw$X.Uploaded_variation
get_vep_status <- function(grna_group, proxy_to_vep){return(paste0(lapply(unlist(str_split(grna_group, "-")), FUN = paste_proxy, proxy_to_vep=proxy_to_vep),collapse = "/"))}
paste_proxy <- function(proxy_char, proxy_to_vep){return(ifelse(is.na(proxy_to_vep[proxy_char]), NA, proxy_to_vep[proxy_char]))}
discovery_result_implicated <- cbind.data.frame(discovery_result_implicated, 
                                                vep_annotation=vapply(discovery_result_implicated$grna_group, FUN = get_vep_status, FUN.VALUE = character(1), proxy_to_vep=proxy_to_vep),
                                                p_value_BH=p.adjust(discovery_result_implicated$p_value, method = "BH"))
#Write to file
discovery_result_implicated$response_id <- ensg_to_symbol[discovery_result_implicated$response_id]
write.table(discovery_result_implicated, file = paste0(out_dir, "discovery_results.implicated_connections_only.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

