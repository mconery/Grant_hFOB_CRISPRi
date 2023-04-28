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

# 0) Call libraries and set directories ====

#Call libraries
library(MAST)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(biomaRt)
library(stringr)
library(zellkonverter)
library(future)
library(DESeq2)
library(plyr)

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/differential_expression/"

#Set parallelization scheme
plan("multicore")
options(future.globals.maxSize = 4 * 1000 * 1024^2)

# 1) Read in files ====

#Read in Scanpy QC-ed cell results
qced_results_raw <- readH5AD(qced_results_loc)
qced_results_seurat <- adata_Seurat <- as.Seurat(qced_results_raw, counts = "X", data = NULL)

#Read in number of guides per cell
guides_per_cell_raw <- read.csv(guides_per_cell_loc, header = TRUE, sep = ",")

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

# 3) Run false discovery MAST Tests vs Out Group ====

#Run false discovery tests for each non-targeting sgRNA in single-sgRNA cells
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
non_targeting_in_out_single_sgrna_results <- lapply(neg_control_sgrnas, 
                                               test_mast_in_out_single_sgrna_cells, 
                                               qced_results_mast=qced_results_mast, 
                                               guide_to_single_guide_cells=guide_to_single_guide_cells,
                                               single_guide_cells=single_guide_cells)
names(non_targeting_in_out_single_sgrna_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_in_out_single_sgrna_df <- rbind.fill.matrix(non_targeting_in_out_single_sgrna_results)
#Save results as an RDS object
saveRDS(non_targeting_in_out_single_sgrna_df, file = paste0(out_dir, "non_targeting_in_out_single_sgrna_df.rds"))

#Run the false discovery tests for each non-targeting sgRNA in all cells
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
non_targeting_in_out_all_results <- lapply(neg_control_sgrnas, 
                                                    test_mast_in_out_single_sgrna_cells, 
                                                    qced_results_mast=qced_results_mast, 
                                                    guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                    single_guide_cells=single_guide_cells)
names(non_targeting_in_out_all_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_in_out_all_df <- rbind.fill.matrix(non_targeting_in_out_all_results)
#Save results as an RDS object
saveRDS(non_targeting_in_out_all_df, file = paste0(out_dir, "non_targeting_in_out_all_df.rds"))

# 4) Run false discovery MAST Tests vs Remaining Non-Targeting Controls ====

#Run the single-guide test first
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
non_targeting_ntc_single_sgrna_results <- lapply(neg_control_sgrnas, 
                                                    test_mast_ntc_single_sgrna_cells, 
                                                    qced_results_mast=qced_results_mast, 
                                                    guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                    neg_control_sgrnas=neg_control_sgrnas)
names(non_targeting_ntc_single_sgrna_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_ntc_single_sgrna_df <- rbind.fill.matrix(non_targeting_ntc_single_sgrna_results)
#Save results as an RDS object
saveRDS(non_targeting_ntc_single_sgrna_df, file = paste0(out_dir, "non_targeting_ntc_single_sgrna_df.rds"))

#Run the multi-guide test 
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
non_targeting_ntc_all_results <- lapply(neg_control_sgrnas, 
                                         test_mast_ntc_all_cells, 
                                         qced_results_mast=qced_results_mast, 
                                         guide_to_cells=guide_to_cells,
                                         neg_control_sgrnas=neg_control_sgrnas)
names(non_targeting_ntc_all_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_ntc_all_df <- rbind.fill.matrix(non_targeting_ntc_all_results)
#Save results as an RDS object
saveRDS(non_targeting_ntc_all_df, file = paste0(out_dir, "non_targeting_ntc_all_df.rds"))
