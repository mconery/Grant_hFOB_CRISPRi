################################################################################

#test_differential_expression.R

#The purpose of this script is to test the cells containing a given sgRNA for  
#differential expression of nearby genes. The script is set to run with the 
#following modules:

#  R/4.2.2
#  hdf5/1.10.1
#  Pandoc/2.10

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
library(zellkonverter)
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

# 1) Read in files ====

#Read in Scanpy QC-ed cell results
qced_results_raw <- readH5AD(qced_results_loc)
qced_results_seurat <- adata_Seurat <- as.Seurat(qced_results_raw, counts = "X", data = NULL)

#Read in number of guides per cell
guides_per_cell_raw <- read.csv(guides_per_cell_loc, header = TRUE, sep = ",")

#Read in target locations file 
sgrna_target_raw <- read.table(sgrna_target_loc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sgrna_target_raw) <- c("chr", "start", "end", "name+seq")

#Read in genes file
gencode_genes_raw <- read.table(gencode_genes_loc, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gencode_genes_raw) <- c("chr", "start", "end", "ENSG+name", "gene_type", "strand")

# 2) Prepare to Run Tests ====

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

# 3) Get list of genes to test within range of sgrna targets ====

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

# 4) Define Testing Functions ====

#Create function for testing each sgRNA using in-out method in single-sgRNA cells
test_mast_in_out_single_sgrna_cells <- function(sgrna, qced_results_mast, guide_to_single_guide_cells, single_guide_cells, test_features = NULL){
  cells_with_guide <- guide_to_single_guide_cells[sgrna][[1]]
  comparison_cells <- single_guide_cells[!(single_guide_cells %in% cells_with_guide)]
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

#Create function for testing each sgRNA using in-out method in all cells
test_mast_in_out_all_cells <- function(sgrna, qced_results_mast, guide_to_cells, all_cells, test_features = NULL){
  cells_with_guide <- guide_to_cells[sgrna][[1]]
  comparison_cells <- all_cells[!(all_cells %in% cells_with_guide)]
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

#Create function for testing each sgRNA using ntc method in single-sgrna cells
test_mast_ntc_single_sgrna_cells <- function(sgrna, qced_results_mast, guide_to_single_guide_cells, neg_control_sgrnas, test_features = NULL){
  cells_with_guide <- guide_to_single_guide_cells[sgrna][[1]]
  remain_neg_control_sgrnas <- neg_control_sgrnas[neg_control_sgrnas != sgrna]
  comparison_cells <- unique(unlist(guide_to_single_guide_cells[remain_neg_control_sgrnas]))
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

#Create function for testing each sgRNA using ntc method in all cells
test_mast_ntc_all_cells <- function(sgrna, qced_results_mast, guide_to_cells, neg_control_sgrnas, test_features = NULL){
  cells_with_guide <- guide_to_cells[sgrna][[1]]
  remain_neg_control_sgrnas <- neg_control_sgrnas[neg_control_sgrnas != sgrna]
  comparison_cells <- unique(unlist(guide_to_cells[remain_neg_control_sgrnas]))
  comparison_cells <- comparison_cells[!(comparison_cells %in% cells_with_guide)]
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

# 5) Run false discovery MAST Tests (Commented Out Since Previously Run) ====

##Run false discovery tests for each non-targeting sgRNA using in-out method in single-sgRNA cells
#non_targeting_in_out_single_sgrna_results <- lapply(neg_control_sgrnas, 
#                                                    test_mast_in_out_single_sgrna_cells, 
#                                                    qced_results_mast=qced_results_ntc_test, 
#                                                    guide_to_single_guide_cells=guide_to_single_guide_cells,
#                                                    single_guide_cells=single_guide_cells)
#names(non_targeting_in_out_single_sgrna_results) <- neg_control_sgrnas
##Combine results into a matrix
#non_targeting_in_out_single_sgrna_df <- rbind.fill.matrix(non_targeting_in_out_single_sgrna_results)
##Save results as an RDS object
#saveRDS(non_targeting_in_out_single_sgrna_df, file = paste0(out_dir, "non_targeting_in_out_single_sgrna_df.rds"))
#
##Run false discovery tests for each non-targeting sgRNA using in-out method in all cells
#non_targeting_in_out_all_results <- lapply(neg_control_sgrnas, 
#                                           test_mast_in_out_single_sgrna_cells, 
#                                           qced_results_mast=qced_results_ntc_test, 
#                                           guide_to_single_guide_cells=guide_to_single_guide_cells,
#                                           single_guide_cells=single_guide_cells)
#names(non_targeting_in_out_all_results) <- neg_control_sgrnas
##Combine results into a matrix
#non_targeting_in_out_all_df <- rbind.fill.matrix(non_targeting_in_out_all_results)
##Save results as an RDS object
#saveRDS(non_targeting_in_out_all_df, file = paste0(out_dir, "non_targeting_in_out_all_df.rds"))
#
##Run false discovery tests for each non-targeting sgRNA using ntc method in single-sgRNA cells
#non_targeting_ntc_single_sgrna_results <- lapply(neg_control_sgrnas, 
#                                                    test_mast_ntc_single_sgrna_cells, 
#                                                    qced_results_mast=qced_results_ntc_test, 
#                                                    guide_to_single_guide_cells=guide_to_single_guide_cells,
#                                                    neg_control_sgrnas=neg_control_sgrnas)
#names(non_targeting_ntc_single_sgrna_results) <- neg_control_sgrnas
##Combine results into a matrix
#non_targeting_ntc_single_sgrna_df <- rbind.fill.matrix(non_targeting_ntc_single_sgrna_results)
##Save results as an RDS object
#saveRDS(non_targeting_ntc_single_sgrna_df, file = paste0(out_dir, "non_targeting_ntc_single_sgrna_df.rds"))
#
##Run false discovery tests for each non-targeting sgRNA using ntc method in all cells
#non_targeting_ntc_all_results <- lapply(neg_control_sgrnas, 
#                                         test_mast_ntc_all_cells, 
#                                         qced_results_mast=qced_results_ntc_test, 
#                                         guide_to_cells=guide_to_cells,
#                                         neg_control_sgrnas=neg_control_sgrnas)
#names(non_targeting_ntc_all_results) <- neg_control_sgrnas
##Combine results into a matrix
#non_targeting_ntc_all_df <- rbind.fill.matrix(non_targeting_ntc_all_results)
##Save results as an RDS object
#saveRDS(non_targeting_ntc_all_df, file = paste0(out_dir, "non_targeting_ntc_all_df.rds"))
# 6) Run Positive Control Tests (Commented Out Since Previously Run) ====

##Test positive control genes in single-guide rna cells with in-out approach
#pos_control_in_out_single_sgrna_results <- lapply(pos_control_sgrnas, 
#                                                  test_mast_in_out_single_sgrna_cells, 
#                                                  qced_results_mast=qced_results_mast[symbol_to_ensg[pos_control_sgrnas],], 
#                                                  guide_to_single_guide_cells=guide_to_single_guide_cells,
#                                                  single_guide_cells=single_guide_cells)
#names(pos_control_in_out_single_sgrna_results) <- pos_control_sgrnas
##Combine results into a matrix
#pos_control_in_out_single_sgrna_df <- rbind.fill.matrix(pos_control_in_out_single_sgrna_results)
##Overwrite gene names and filter for relevant tests
#pos_control_in_out_single_sgrna_df[,"gene"] <- ensg_to_symbol[pos_control_in_out_single_sgrna_df[,"gene"]]
#pos_control_in_out_single_sgrna_df <- pos_control_in_out_single_sgrna_df[which(pos_control_in_out_single_sgrna_df[,"gene"] == pos_control_in_out_single_sgrna_df[,"sgrna"]),]
##Save RDS object
#saveRDS(pos_control_in_out_single_sgrna_df, file = paste0(out_dir, "pos_control_in_out_single_sgrna_df.rds"))
#
##Test positive control genes in all cells with in-out approach
#pos_control_in_out_all_results <- lapply(pos_control_sgrnas, 
#                                                  test_mast_in_out_all_cells, 
#                                                  qced_results_mast=qced_results_mast[symbol_to_ensg[pos_control_sgrnas],], 
#                                                  guide_to_cells=guide_to_cells,
#                                                  all_cells=all_cells)
#names(pos_control_in_out_all_results) <- pos_control_sgrnas
##Combine results into a matrix
#pos_control_in_out_all_df <- rbind.fill.matrix(pos_control_in_out_all_results)
##Overwrite gene names and filter for relevant tests
#pos_control_in_out_all_df[,"gene"] <- ensg_to_symbol[pos_control_in_out_all_df[,"gene"]]
#pos_control_in_out_all_df <- pos_control_in_out_all_df[which(pos_control_in_out_all_df[,"gene"] == pos_control_in_out_all_df[,"sgrna"]),]
##Save results as an RDS object
#saveRDS(pos_control_in_out_all_df, file = paste0(out_dir, "pos_control_in_out_all_df.rds"))
#
##Test positive control genes in single-guide rna cells with non-targeting approach
#pos_control_ntc_single_sgrna_results <- lapply(pos_control_sgrnas, 
#                                                  test_mast_ntc_single_sgrna_cells, 
#                                                  qced_results_mast=qced_results_mast[symbol_to_ensg[pos_control_sgrnas],], 
#                                                  guide_to_single_guide_cells=guide_to_single_guide_cells,
#                                                  neg_control_sgrnas=neg_control_sgrnas)
#names(pos_control_ntc_single_sgrna_results) <- pos_control_sgrnas
##Combine results into a matrix
#pos_control_ntc_single_sgrna_df <- rbind.fill.matrix(pos_control_ntc_single_sgrna_results)
##Overwrite gene names and filter for relevant tests
#pos_control_ntc_single_sgrna_df[,"gene"] <- ensg_to_symbol[pos_control_ntc_single_sgrna_df[,"gene"]]
#pos_control_ntc_single_sgrna_df <- pos_control_ntc_single_sgrna_df[which(pos_control_ntc_single_sgrna_df[,"gene"] == pos_control_ntc_single_sgrna_df[,"sgrna"]),]
##Save RDS object
#saveRDS(pos_control_ntc_single_sgrna_df, file = paste0(out_dir, "pos_control_ntc_single_sgrna_df.rds"))
#
##Test positive control genes in all cells with non-targeting approach
#pos_control_ntc_all_results <- lapply(pos_control_sgrnas, 
#                                         test_mast_ntc_all_cells, 
#                                         qced_results_mast=qced_results_mast[symbol_to_ensg[pos_control_sgrnas],], 
#                                         guide_to_cells=guide_to_cells,
#                                         neg_control_sgrnas=neg_control_sgrnas)
#names(pos_control_ntc_all_results) <- pos_control_sgrnas
##Combine results into a matrix
#pos_control_ntc_all_df <- rbind.fill.matrix(pos_control_ntc_all_results)
##Overwrite gene names and filter for relevant tests
#pos_control_ntc_all_df[,"gene"] <- ensg_to_symbol[pos_control_ntc_all_df[,"gene"]]
#pos_control_ntc_all_df <- pos_control_ntc_all_df[which(pos_control_ntc_all_df[,"gene"] == pos_control_ntc_all_df[,"sgrna"]),]
##Save results as an RDS object
#saveRDS(pos_control_ntc_all_df, file = paste0(out_dir, "pos_control_ntc_all_df.rds"))

# 7) Run Enhancer Targeting Test (Commented Out Since Previously Run) ====

##Make a map for each guide to the list of genes we are testing it against
#check_genes_for_guide <- function(sgrna_target_row, gencode_genes_filter){
#  chromo <- as.character(sgrna_target_row["chr"])
#  rep_start <- as.numeric(sgrna_target_row["rep_start"])
#  rep_end <- as.numeric(sgrna_target_row["rep_end"])
#  gencode_genes_temp <- gencode_genes_filter[which(gencode_genes_filter$chr == chromo),]
#  gencode_genes_temp <- gencode_genes_temp[which((gencode_genes_temp$start >= rep_start & gencode_genes_temp$start <= rep_end) | 
#                                                 (gencode_genes_temp$end >= rep_start & gencode_genes_temp$end <= rep_end) | 
#                                                 (gencode_genes_temp$start <= rep_start & gencode_genes_temp$end >= rep_end)),]
#  return(str_replace(gencode_genes_temp$ensembl_gene_id, "[.][0-9]+",""))
#}
#guide_to_genes <- pbapply(sgrna_target_filter, FUN = check_genes_for_guide, MARGIN = 1, gencode_genes_filter=gencode_genes_filter)
#
##Create a wrapper function so that we can apply the in-out tests to the targeting guides
#run_targeting_in_out_test <- function(sgrna, test_mast_in_out, qced_results_mast, guide_to_cells, guide_to_genes, cell_list){
#  qced_results_test <- qced_results_mast[unlist(guide_to_genes[sgrna]),]
#  return(test_mast_in_out(sgrna, qced_results_test, guide_to_cells, cell_list))
#}
#
##Run in-out test for single-sgrna cells
#targeting_in_out_single_sgrna_results <- lapply(test_sgrnas, 
#                                                FUN=run_targeting_in_out_test,
#                                                test_mast_in_out=test_mast_in_out_single_sgrna_cells, 
#                                                qced_results_mast=qced_results_mast, 
#                                                guide_to_cells=guide_to_single_guide_cells,
#                                                guide_to_genes=guide_to_genes,
#                                                cell_list=single_guide_cells)
#names(targeting_in_out_single_sgrna_results) <- test_sgrnas
##Combine results into a matrix
#targeting_in_out_single_sgrna_df <- rbind.fill.matrix(targeting_in_out_single_sgrna_results)
##Overwrite gene names
#targeting_in_out_single_sgrna_df[,"gene"] <- ensg_to_symbol[targeting_in_out_single_sgrna_df[,"gene"]]
##Save RDS object
#saveRDS(targeting_in_out_single_sgrna_df, file = paste0(out_dir, "targeting_in_out_single_sgrna_df.rds"))
#
##Run in-out test for all cells
#targeting_in_out_all_results <- lapply(test_sgrnas, 
#                                                FUN=run_targeting_in_out_test,
#                                                test_mast_in_out=test_mast_in_out_all_cells, 
#                                                qced_results_mast=qced_results_mast, 
#                                                guide_to_cells=guide_to_cells,
#                                                guide_to_genes=guide_to_genes,
#                                                cell_list=all_cells)
#names(targeting_in_out_all_results) <- test_sgrnas
##Combine results into a matrix
#targeting_in_out_all_df <- rbind.fill.matrix(targeting_in_out_all_results)
##Overwrite gene names
#targeting_in_out_all_df[,"gene"] <- ensg_to_symbol[targeting_in_out_all_df[,"gene"]]
##Save RDS object
#saveRDS(targeting_in_out_all_df, file = paste0(out_dir, "targeting_in_out_all_df.rds"))
#
##Create a wrapper function so that we can apply the ntc tests to the targeting guides
#run_targeting_ntc_test <- function(sgrna, test_mast_in_out, qced_results_mast, guide_to_cells, guide_to_genes, neg_control_sgrnas){
#  qced_results_test <- qced_results_mast[unlist(guide_to_genes[sgrna]),]
#  return(test_mast_in_out(sgrna, qced_results_test, guide_to_cells, neg_control_sgrnas))
#}                                         
#
##Run ntc test for single-sgrna cells
#targeting_ntc_single_sgrna_results <- lapply(test_sgrnas, 
#                                                FUN=run_targeting_ntc_test,
#                                                test_mast_in_out=test_mast_ntc_single_sgrna_cells, 
#                                                qced_results_mast=qced_results_mast, 
#                                                guide_to_cells=guide_to_single_guide_cells,
#                                                guide_to_genes=guide_to_genes,
#                                                neg_control_sgrnas=neg_control_sgrnas)
#names(targeting_ntc_single_sgrna_results) <- test_sgrnas
##Combine results into a matrix
#targeting_ntc_single_sgrna_df <- rbind.fill.matrix(targeting_ntc_single_sgrna_results)
##Overwrite gene names
#targeting_ntc_single_sgrna_df[,"gene"] <- ensg_to_symbol[targeting_ntc_single_sgrna_df[,"gene"]]
##Save RDS object
#saveRDS(targeting_ntc_single_sgrna_df, file = paste0(out_dir, "targeting_ntc_single_sgrna_df.rds"))
#
##Run ntc test for all cells
#targeting_ntc_all_results <- lapply(test_sgrnas, 
#                                             FUN=run_targeting_ntc_test,
#                                             test_mast_in_out=test_mast_ntc_all_cells, 
#                                             qced_results_mast=qced_results_mast, 
#                                             guide_to_cells=guide_to_cells,
#                                             guide_to_genes=guide_to_genes,
#                                             neg_control_sgrnas=neg_control_sgrnas)
#names(targeting_ntc_all_results) <- test_sgrnas
##Combine results into a matrix
#targeting_ntc_all_df <- rbind.fill.matrix(targeting_ntc_all_results)
##Overwrite gene names
#targeting_ntc_all_df[,"gene"] <- ensg_to_symbol[targeting_ntc_all_df[,"gene"]]
##Save RDS object
#saveRDS(targeting_ntc_all_df, file = paste0(out_dir, "targeting_ntc_all_df.rds"))


# 8) Create Null Distribution Testing Functions ====

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

# 9) Run Null Distributions test

#Run null tests for each non-targeting sgRNA using in-out method in single-sgRNA cells
null_in_out_single_sgrna_results <- lapply(seq(1,1000,1), 
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

#Run null tests for each non-targeting sgRNA using in-out method in all cells
null_in_out_all_results <- lapply(seq(1,1000,1), 
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

#Get list of non-targeting cells
single_sgrna_ntc_cells <- unique(unlist(guide_to_single_guide_cells[neg_control_sgrnas]))
ntc_only_cells <- unlist(guide_to_cells[neg_control_sgrnas])
ntc_only_cells <- ntc_only_cells[!(ntc_only_cells %in% unlist(guide_to_cells[pos_control_sgrnas]))]
ntc_only_cells <- unique(ntc_only_cells[!(ntc_only_cells %in% unlist(guide_to_cells[test_sgrnas]))])

#Run null tests for each non-targeting sgRNA using ntc method in single-sgRNA cells
null_ntc_single_sgrna_results <- lapply(seq(1,1000,1), 
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

#Run null tests for each non-targeting sgRNA using ntc method in all cells
null_ntc_all_results <- lapply(seq(1,1000,1), 
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

# 10) Calculate false discovery Test Rates ====

##Read in RDS objects
#non_targeting_in_out_all_df <- readRDS(paste0(out_dir, "non_targeting_in_out_all_df.rds"))
#non_targeting_in_out_single_sgrna_df <- readRDS(paste0(out_dir, "non_targeting_in_out_single_sgrna_df.rds"))
#non_targeting_ntc_all_df <- readRDS(paste0(out_dir, "non_targeting_ntc_all_df.rds"))
#non_targeting_ntc_single_sgrna_df <- readRDS(paste0(out_dir, "non_targeting_ntc_single_sgrna_df.rds"))

##Recast matrix test results to dataframes
#non_targeting_in_out_all_df <- as.data.frame(non_targeting_in_out_all_df)
#non_targeting_in_out_single_sgrna_df <- as.data.frame(non_targeting_in_out_single_sgrna_df)
#non_targeting_ntc_all_df <- as.data.frame(non_targeting_ntc_all_df)
#non_targeting_ntc_single_sgrna_df <- as.data.frame(non_targeting_ntc_single_sgrna_df)

##Recast p-value columns as numerics
#non_targeting_in_out_all_df$p_val <- as.numeric(non_targeting_in_out_all_df$p_val)
#non_targeting_in_out_single_sgrna_df$p_val <- as.numeric(non_targeting_in_out_single_sgrna_df$p_val)
#non_targeting_ntc_all_df$p_val <- as.numeric(non_targeting_ntc_all_df$p_val)
#non_targeting_ntc_single_sgrna_df$p_val <- as.numeric(non_targeting_ntc_single_sgrna_df$p_val)

##Benjamini-hochberg correct p-values
#non_targeting_in_out_all_df$p_val_adj <- p.adjust(non_targeting_in_out_all_df$p_val, method="BH")
#non_targeting_in_out_single_sgrna_df$p_val_adj <- p.adjust(non_targeting_in_out_single_sgrna_df$p_val, method="BH")
#non_targeting_ntc_all_df$p_val_adj <- p.adjust(non_targeting_ntc_all_df$p_val, method="BH")
#non_targeting_ntc_single_sgrna_df$p_val_adj <- p.adjust(non_targeting_ntc_single_sgrna_df$p_val, method="BH")

##Make qqplots
#make_non_targeting_qq <- function(non_targeting_df){
#  
#}

