################################################################################

#test_differential_expression.R

#The purpose of this script is to test the cells containing a given sgRNA for  
#differential expression of nearby genes. The script is set to run with the 
#following modules:

#  R/4.2.3

#The script will try two main different methods for detecting differential 
#expression of genes:

# 1) Compare a sgRNA to all cells without the sgRNA and only receiving 1 sgRNA
# 2) Compare a sgRNA to all cells without the sgRNA

################################################################################

# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(MAST)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(stringr)
library(future)
library(plyr)
library(pbapply)
library(bedr)
library(biomaRt)

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
gencode_genes_loc <- "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.annotation.gene_only.bed"
sgrna_target_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/all_targeting_guides.bed"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/differential_expression/"

#Set parallelization scheme
plan("multicore")
options(future.globals.maxSize = 4 * 1000 * 1024^2)

#Set repression range
rep_range <- 1000000 #Testing within +/- 1Mb

# 1) Read in variables ====

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 2) {
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  test_type <- args[1]
  cell_type <- as.integer(args[2])
} else {
  print("Usage: %> Rscript test_differential_expression_null.R test_type cell_type");
  quit(save="no");
}

# 2) Read in files ====

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

# 3) Prepare to Run Tests ====

#Get list of unique sgrna targets
sgrnas <- unique(guides_per_cell_raw$feature_call)
sgrnas <- sgrnas[str_detect(sgrnas, "[|]") == FALSE]
#Make maps for each sgRNA to the indices of cells that contain it both as a single sgRNA and as a multiple
find_cells_with_guide <- function(sgrna, guides_per_cell_raw){
  return(guides_per_cell_raw$cell_barcode[str_detect(guides_per_cell_raw$feature_call, sgrna)])
}
guide_to_cells <- lapply(sgrnas, FUN = find_cells_with_guide, guides_per_cell_raw=guides_per_cell_raw)
names(guide_to_cells) <- sgrnas
#Get list of single guide cells and all cells
all_cells <- guides_per_cell_raw$cell_barcode
single_guide_cells <- guides_per_cell_raw$cell_barcode[str_detect(guides_per_cell_raw$feature_call, "[|]") == FALSE]
#Make a map for the single guide cells
find_cells_with_single_guide <- function(sgrna, guides_per_cell_raw){
  return(guides_per_cell_raw$cell_barcode[str_detect(guides_per_cell_raw$feature_call, sgrna) & (str_detect(guides_per_cell_raw$feature_call, "[|]") == FALSE)])
}
guide_to_single_guide_cells <- lapply(sgrnas, FUN = find_cells_with_single_guide, guides_per_cell_raw=guides_per_cell_raw)
names(guide_to_single_guide_cells) <- sgrnas

#Divy up the guide rnas into their respective categories
pos_control_sgrnas <- sgrnas[1:2]
neg_control_sgrnas <- sgrnas[3:29]
test_sgrnas <- sgrnas[30:length(sgrnas)]

#Log normalize the data for the mast test
qced_results_mast <- NormalizeData(qced_results_seurat, normalization.method = "LogNormalize")

# 4) Get list of genes to test within range of sgrna targets ====

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
gencode_genes_filter <- gencode_genes_append[which(str_replace(gencode_genes_append$ensembl_gene_id, "[.][0-9]+", "") %in% rownames((qced_results_mast))),]

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

#Filter down the qced results for just the genes that are nearby the target genes 
qced_results_ntc_test <- qced_results_mast[test_genes,] 

#Create two lists for converting gene symbols to ensg and vice versa
symbol_to_ensg <- rownames(gencode_genes_filter)
names(symbol_to_ensg) <- gencode_genes_filter$external_gene_name
ensg_to_symbol <- gencode_genes_filter$external_gene_name
names(ensg_to_symbol) <- rownames(gencode_genes_filter)

# 5) Create Null Distribution Testing Functions ====

#Create function for testing each sgRNA using in-out method
test_random_mast <- function(sgrna, qced_results_mast, guide_to_cells, cell_group, test_features = NULL){
  num_cells_with_guide <- length(guide_to_cells[sgrna][[1]])
  cells_with_guide <- sample(cell_group, num_cells_with_guide, replace = FALSE)#Select random cells
  comparison_cells <- cell_group[!(cell_group %in% cells_with_guide)]
  qced_results_test <- qced_results_mast
  qced_results_test <- SetIdent(qced_results_test, cells = cells_with_guide, value = sgrna)
  qced_results_test <- SetIdent(qced_results_test, cells = comparison_cells, value = "out_group")
  test <- FindMarkers(qced_results_test, slot = "data",
                      ident.1 = sgrna, ident.2 = "out_group",
                      test.use='MAST', logfc.threshold=0, 
                      min.cells.feature=0, min.cells.group=0, 
                      latent.vars = "Pool", features = test_features)
  test <- as.data.frame(test)
  test <- cbind.data.frame(sgrna=rep(sgrna, nrow(test)), gene=rownames(test), test)
  return(test)
}

#Make a wrapper function
test_random_mast_wrapper <- function(iteration, neg_control_sgrnas, qced_results_ntc_test, guide_to_cells, cell_group){
  #Lapply test function 
  null_results <- lapply(neg_control_sgrnas, 
                           test_random_mast, 
                           qced_results_mast=qced_results_ntc_test, 
                           guide_to_cells=guide_to_cells,
                           cell_group=cell_group)
  names(null_results) <- neg_control_sgrnas
  #Combine results into a matrix
  null_df <- rbind.fill.matrix(null_results)
  #Append itertion colunn onto dataframe
  null_df <- cbind.data.frame(iteration=rep(iteration, nrow(null_df)), null_df)
  #print iteration
  print(paste0("Iteration ", iteration, " complete."))
  #Retun dataframe
  return(null_df)
}

# 6) Run Null Distributions test ====

#Get list of non-targeting cells
single_sgrna_ntc_cells <- unique(unlist(guide_to_single_guide_cells[neg_control_sgrnas]))
ntc_only_cells <- unlist(guide_to_cells[neg_control_sgrnas])
ntc_only_cells <- ntc_only_cells[!(ntc_only_cells %in% unlist(guide_to_cells[pos_control_sgrnas]))]
ntc_only_cells <- unique(ntc_only_cells[!(ntc_only_cells %in% unlist(guide_to_cells[test_sgrnas]))])

#Run null tests
if (test_type == "in_out" && cell_type == "single_guide") {
  #Run null tests for each non-targeting sgRNA using in-out method in single-sgRNA cells
  null_in_out_single_sgrna_results <- lapply(seq(1,100,1), 
                                             FUN = test_random_mast_wrapper, 
                                             neg_control_sgrnas=neg_control_sgrnas,
                                             qced_results_ntc_test=qced_results_ntc_test, 
                                             guide_to_cells=guide_to_single_guide_cells,
                                             cell_group=single_guide_cells)
  names(null_in_out_single_sgrna_results) <- neg_control_sgrnas
  #Combine results into a matrix
  null_in_out_single_sgrna_df <- rbind.fill.matrix(null_in_out_single_sgrna_results)
  #Save results as an RDS object
  saveRDS(null_in_out_single_sgrna_df, file = paste0(out_dir, "null_in_out_single_sgrna_df.rds"))
} else if (test_type = "in_out" && cell_type == "all"){
  #Run null tests for each non-targeting sgRNA using in-out method in all cells
  null_in_out_all_results <- lapply(seq(1,100,1), 
                                    FUN = test_random_mast_wrapper, 
                                    neg_control_sgrnas=neg_control_sgrnas,
                                    qced_results_ntc_test=qced_results_ntc_test, 
                                    guide_to_cells=guide_to_cells,
                                    cell_group=all_cells)
  names(null_in_out_all_results) <- neg_control_sgrnas
  #Combine results into a matrix
  null_in_out_all_df <- rbind.fill.matrix(null_in_out_all_results)
  #Save results as an RDS object
  saveRDS(null_in_out_all_df, file = paste0(out_dir, "null_in_out_all_df.rds"))
} else if (test_type = "ntc" && cell_type == "single_guide") {
  #Run null tests for each non-targeting sgRNA using ntc method in single-sgRNA cells
  null_ntc_single_sgrna_results <- lapply(seq(1,100,1), 
                                          FUN = test_random_mast_wrapper, 
                                          neg_control_sgrnas=neg_control_sgrnas,
                                          qced_results_ntc_test=qced_results_ntc_test, 
                                          guide_to_cells=guide_to_single_guide_cells,
                                          cell_group=single_sgrna_ntc_cells)
  names(null_ntc_single_sgrna_results) <- neg_control_sgrnas
  #Combine results into a matrix
  null_ntc_single_sgrna_df <- rbind.fill.matrix(null_ntc_single_sgrna_results)
  #Save results as an RDS object
  saveRDS(null_ntc_single_sgrna_df, file = paste0(out_dir, "null_ntc_single_sgrna_df.rds"))
} else if (test_type = "ntc" && cell_type == "all"){
  #Run null tests for each non-targeting sgRNA using ntc method in all cells
  null_ntc_all_results <- lapply(seq(1,100,1), 
                                 FUN = test_random_mast_wrapper, 
                                 neg_control_sgrnas=neg_control_sgrnas,
                                 qced_results_ntc_test=qced_results_ntc_test, 
                                 guide_to_cells=guide_to_cells,
                                 cell_group=ntc_only_cells)
  names(null_ntc_all_results) <- neg_control_sgrnas
  #Combine results into a matrix
  null_ntc_all_df <- rbind.fill.matrix(null_ntc_all_results)
  #Save results as an RDS object
  saveRDS(null_ntc_all_df, file = paste0(out_dir, "null_ntc_all_df.rds"))
} else {
  "ERROR: Invalid cell or test type given."
}







