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

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
gencode_genes_loc <- "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.annotation.gene_only.bed"
sgrna_target_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/all_targeting_guides.bed"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/differential_expression/sceptre_single_sgrna/"

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

# 2) Get list of genes to test within range of sgrna targets ====

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

# 3) Prepare Inputs to SCEPTRE Test ====

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

### Make gRNA group table ###
grna_group_hfob_lowmoi <- rbind.data.frame(cbind.data.frame(grna_id=test_sgrnas, grna_group=str_replace(test_sgrnas, "_[0-9]","")),
                                               cbind.data.frame(grna_id=pos_control_sgrnas, grna_group=pos_control_sgrnas),
                                               cbind.data.frame(grna_id=neg_control_sgrnas, grna_group=rep("non-targeting", length(neg_control_sgrnas))))

# 4) Establish Pairs of Genes/sgRNAs to Test ====

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
response_grna_group_pairs_hfob <- as.data.frame(rbind.fill.matrix(pblapply(unique(str_replace(test_sgrnas, "_[0-9]", "")), FUN = check_genes_for_guide, 
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

# 5) Set the Formula Object ====

#Set the object using the names from the tutorial
formula_object <- formula(~log(response_n_umis) + 
                            log(response_n_nonzero) +
                            bio_rep + 
                            p_mito)

# 6) Run and Assess Calibration Check ====

#Run check
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

# 7) Run Leave One-Out Test for Non-Targeting Control Outliers ====

#Check if directory exists and make it if not
if (!(dir.exists(paste0(out_dir, "leave_one_out")))) {
  dir.create(paste0(out_dir, "leave_one_out"))
}

#Loop over non-targeting controls and rerun calibration analysis
for (sgrna in neg_control_sgrnas) {
  set.seed(5)
  #Remove the selected grna from the grna-by-cell matrix and the grouping file
  grna_hfob_lowmoi_loo <- grna_hfob_lowmoi[!(rownames(grna_hfob_lowmoi) %in% sgrna),]
  grna_group_hfob_lowmoi_loo <- grna_group_hfob_lowmoi[grna_group_hfob_lowmoi$grna_id != sgrna,]
  #Get new gene expression matrix, covariate data frame, and grna matrix for cells (The order is important here!)
  response_hfob_lowmoi_loo <- response_hfob_lowmoi[,colSums(grna_hfob_lowmoi_loo != 0) != 0]
  covariate_hfob_lowmoi_loo <- covariate_hfob_lowmoi[colSums(grna_hfob_lowmoi_loo != 0) != 0,]
  grna_hfob_lowmoi_loo <- grna_hfob_lowmoi_loo[,colSums(grna_hfob_lowmoi_loo != 0) != 0]
  #Run check
  calibration_result_loo <- run_sceptre_lowmoi(
    response_matrix = response_hfob_lowmoi_loo,
    grna_matrix = grna_hfob_lowmoi_loo,
    covariate_data_frame = covariate_hfob_lowmoi_loo,
    grna_group_data_frame = grna_group_hfob_lowmoi_loo,
    formula_object = formula_object,
    response_grna_group_pairs = response_grna_group_pairs_hfob,
    calibration_check = TRUE # calibration_check TRUE
  ) 
  #Plot calibration leave-one-out result
  jpeg(paste0(out_dir, "leave_one_out/", "sceptre_calibration.hfob.loo.", sgrna, ".jpeg"))
  print(plot_calibration_result(calibration_result_loo))
  dev.off()
}

# 8) Run Discovery Analysis ====

#Run discovery
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