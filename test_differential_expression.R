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
library(biomaRt)
library(stringr)
library(zellkonverter)
library(future)
library(DESeq2)
library(plyr)
library(pbapply)

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
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

#Set biomart usage variables
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl", host="https://grch37.ensembl.org")
attributes<-c("hgnc_symbol", 'ensembl_gene_id', 'external_gene_name', "chromosome_name",
              "start_position", "end_position", "strand", "band")

#Get genomic regions as variables
regions <- paste(str_replace(sgrna_target_filter[test_sgrnas,"chr"], "chr", ""),
                 sgrna_target_filter[test_sgrnas,"rep_start"], 
                 sgrna_target_filter[test_sgrnas,"rep_end"], sep = ":")

# extract all regions at once:
target_genes_raw <- as.data.frame(getBM(attributes=attributes,
                          filters=c("chromosomal_region"),
                          values=regions,mart=ensembl))
#Filter target genes for those in seurat matrix and set rownames of target genes to the ensembl ids for easier conversion
target_genes_filter <- target_genes_raw[which(target_genes_raw$ensembl_gene_id %in% rownames(qced_results_mast)),]
target_genes_filter <- target_genes_filter[which(!(target_genes_filter$hgnc_symbol == "SEBOX" & target_genes_filter$ensembl_gene_id == "ENSG00000109072")),] #Manual filter for a gene duplicated with two hgnc symbols but otherwise identical
target_genes_filter <- target_genes_filter[which(!(target_genes_filter$hgnc_symbol != "MIR451B" & target_genes_filter$ensembl_gene_id == "ENSG00000264066")),] #Manual filter for a mirna set retaining only the one where hgnc symbol matches the external gene name
rownames(target_genes_filter) <- target_genes_filter$ensembl_gene_id

#Filter down the qced results for just the genes that are nearby the target genes 
qced_results_ntc_test <- qced_results_mast[target_genes_filter$ensembl_gene_id,] 

# 4) Run false discovery MAST Tests vs Out Group ====

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
                                               qced_results_mast=qced_results_ntc_test, 
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
                                                    qced_results_mast=qced_results_ntc_test, 
                                                    guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                    single_guide_cells=single_guide_cells)
names(non_targeting_in_out_all_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_in_out_all_df <- rbind.fill.matrix(non_targeting_in_out_all_results)
#Save results as an RDS object
saveRDS(non_targeting_in_out_all_df, file = paste0(out_dir, "non_targeting_in_out_all_df.rds"))

# 5) Run false discovery MAST Tests vs Remaining Non-Targeting Controls ====

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
                                                    qced_results_mast=qced_results_ntc_test, 
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
                                         qced_results_mast=qced_results_ntc_test, 
                                         guide_to_cells=guide_to_cells,
                                         neg_control_sgrnas=neg_control_sgrnas)
names(non_targeting_ntc_all_results) <- neg_control_sgrnas
#Combine results into a matrix
non_targeting_ntc_all_df <- rbind.fill.matrix(non_targeting_ntc_all_results)
#Save results as an RDS object
saveRDS(non_targeting_ntc_all_df, file = paste0(out_dir, "non_targeting_ntc_all_df.rds"))

# 6) Run Positive Control Tests ====

#Get ensemble ids for positive control genes
pos_control_conversion <- as.data.frame(getBM(attributes=c("hgnc_symbol", 'ensembl_gene_id'),
                                        filters=c("hgnc_symbol"),
                                        values= c("RAB1A", "SYVN1"),mart=ensembl))
rownames(pos_control_conversion) <- pos_control_conversion$hgnc_symbol
#Make a reverse version for going back the other way
pos_control_inversion <- pos_control_conversion
rownames(pos_control_inversion) <- pos_control_inversion$ensembl_gene_id

#Save results as an RDS object#Test positive control genes in single-guide rna cells with in-out approach
pos_control_in_out_single_sgrna_results <- lapply(pos_control_sgrnas, 
                                                  test_mast_in_out_single_sgrna_cells, 
                                                  qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                                  guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                  single_guide_cells=single_guide_cells)
names(pos_control_in_out_single_sgrna_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_in_out_single_sgrna_df <- rbind.fill.matrix(pos_control_in_out_single_sgrna_results)
#Overwrite gene names and filter for relevant tests
pos_control_in_out_single_sgrna_df[,"gene"] <- pos_control_inversion[pos_control_in_out_single_sgrna_df[,"gene"], "hgnc_symbol"]
pos_control_in_out_single_sgrna_df <- pos_control_in_out_single_sgrna_df[which(pos_control_in_out_single_sgrna_df[,"gene"] == pos_control_in_out_single_sgrna_df[,"sgrna"]),]
#Save RDS object
saveRDS(pos_control_in_out_single_sgrna_df, file = paste0(out_dir, "pos_control_in_out_single_sgrna_df.rds"))

#Test positive control genes in all cells with in-out approach
pos_control_in_out_all_results <- lapply(pos_control_sgrnas, 
                                                  test_mast_in_out_all_cells, 
                                                  qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                                  guide_to_cells=guide_to_cells,
                                                  all_cells=all_cells)
names(pos_control_in_out_all_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_in_out_all_df <- rbind.fill.matrix(pos_control_in_out_all_results)
#Overwrite gene names and filter for relevant tests
pos_control_in_out_all_df[,"gene"] <- pos_control_inversion[pos_control_in_out_all_df[,"gene"], "hgnc_symbol"]
pos_control_in_out_all_df <- pos_control_in_out_all_df[which(pos_control_in_out_all_df[,"gene"] == pos_control_in_out_all_df[,"sgrna"]),]
#Save results as an RDS object
saveRDS(pos_control_in_out_all_df, file = paste0(out_dir, "pos_control_in_out_all_df.rds"))

#Test positive control genes in single-guide rna cells with in-out approach
pos_control_in_out_single_sgrna_results <- lapply(pos_control_sgrnas, 
                                                  test_mast_in_out_single_sgrna_cells, 
                                                  qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                                  guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                  single_guide_cells=single_guide_cells)
names(pos_control_in_out_single_sgrna_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_in_out_single_sgrna_df <- rbind.fill.matrix(pos_control_in_out_single_sgrna_results)
#Overwrite gene names and filter for relevant tests
pos_control_in_out_single_sgrna_df[,"gene"] <- pos_control_inversion[pos_control_in_out_single_sgrna_df[,"gene"], "hgnc_symbol"]
pos_control_in_out_single_sgrna_df <- pos_control_in_out_single_sgrna_df[which(pos_control_in_out_single_sgrna_df[,"gene"] == pos_control_in_out_single_sgrna_df[,"sgrna"]),]
#Save RDS object
saveRDS(pos_control_in_out_single_sgrna_df, file = paste0(out_dir, "pos_control_in_out_single_sgrna_df.rds"))

#Test positive control genes in all cells with in-out approach
pos_control_in_out_all_results <- lapply(pos_control_sgrnas, 
                                                  test_mast_in_out_all_cells, 
                                                  qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                                  guide_to_cells=guide_to_cells,
                                                  all_cells=all_cells)
names(pos_control_in_out_all_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_in_out_all_df <- rbind.fill.matrix(pos_control_in_out_all_results)
#Overwrite gene names and filter for relevant tests
pos_control_in_out_all_df[,"gene"] <- pos_control_inversion[pos_control_in_out_all_df[,"gene"], "hgnc_symbol"]
pos_control_in_out_all_df <- pos_control_in_out_all_df[which(pos_control_in_out_all_df[,"gene"] == pos_control_in_out_all_df[,"sgrna"]),]
#Save results as an RDS object
saveRDS(pos_control_in_out_all_df, file = paste0(out_dir, "pos_control_in_out_all_df.rds"))

#Test positive control genes in single-guide rna cells with non-targeting approach
pos_control_ntc_single_sgrna_results <- lapply(pos_control_sgrnas, 
                                                  test_mast_ntc_single_sgrna_cells, 
                                                  qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                                  guide_to_single_guide_cells=guide_to_single_guide_cells,
                                                  neg_control_sgrnas=neg_control_sgrnas)
names(pos_control_ntc_single_sgrna_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_ntc_single_sgrna_df <- rbind.fill.matrix(pos_control_ntc_single_sgrna_results)
#Overwrite gene names and filter for relevant tests
pos_control_ntc_single_sgrna_df[,"gene"] <- pos_control_inversion[pos_control_ntc_single_sgrna_df[,"gene"], "hgnc_symbol"]
pos_control_ntc_single_sgrna_df <- pos_control_ntc_single_sgrna_df[which(pos_control_ntc_single_sgrna_df[,"gene"] == pos_control_ntc_single_sgrna_df[,"sgrna"]),]
#Save RDS object
saveRDS(pos_control_ntc_single_sgrna_df, file = paste0(out_dir, "pos_control_ntc_single_sgrna_df.rds"))

#Test positive control genes in all cells with non-targeting approach
pos_control_ntc_all_results <- lapply(pos_control_sgrnas, 
                                         test_mast_ntc_all_cells, 
                                         qced_results_mast=qced_results_mast[pos_control_conversion[pos_control_sgrnas,"ensembl_gene_id"],], 
                                         guide_to_cells=guide_to_cells,
                                         neg_control_sgrnas=neg_control_sgrnas)
names(pos_control_ntc_all_results) <- pos_control_sgrnas
#Combine results into a matrix
pos_control_ntc_all_df <- rbind.fill.matrix(pos_control_ntc_all_results)
#Overwrite gene names and filter for relevant tests
pos_control_ntc_all_df[,"gene"] <- pos_control_inversion[pos_control_ntc_all_df[,"gene"], "hgnc_symbol"]
pos_control_ntc_all_df <- pos_control_ntc_all_df[which(pos_control_ntc_all_df[,"gene"] == pos_control_ntc_all_df[,"sgrna"]),]
#Save results as an RDS object
saveRDS(pos_control_ntc_all_df, file = paste0(out_dir, "pos_control_ntc_all_df.rds"))

# 7) Run Enhancer Targeting Test ====

#Make a map for each guide to the list of genes we are testing it against
check_genes_for_guide <- function(sgrna_target_row, target_genes_filter){
  chromo <- as.integer(str_replace(as.character(sgrna_target_row["chr"]),"chr", ""))
  rep_start <- as.numeric(sgrna_target_row["rep_start"])
  rep_end <- as.numeric(sgrna_target_row["rep_end"])
  target_genes_temp <- target_genes_filter[which(target_genes_filter$chromosome_name == chromo),]
  target_genes_temp <- target_genes_temp[which((target_genes_temp$start_position >= rep_start & target_genes_temp$start_position <= rep_end) | 
                                                 (target_genes_temp$end_position >= rep_start & target_genes_temp$end_position <= rep_end) | 
                                                 (target_genes_temp$start_position <= rep_start & target_genes_temp$end_position >= rep_end)),]
  return(target_genes_temp$ensembl_gene_id)
}
guide_to_genes <- pbapply(sgrna_target_filter, FUN = check_genes_for_guide, MARGIN = 1, target_genes_filter=target_genes_filter)

#Create a wrapper function so that we can apply the in-out tests to the targeting guides
run_targeting_in_out_test <- function(sgrna, test_mast_in_out, qced_results_mast, guide_to_cells, guide_to_genes, cell_list){
  qced_results_test <- qced_results_mast[unlist(guide_to_genes[sgrna]),]
  return(test_mast_in_out(sgrna, qced_results_test, guide_to_cells, cell_list))
}

#Run in-out test for single-sgrna cells
targeting_in_out_single_sgrna_results <- lapply(test_sgrnas, 
                                                FUN=run_targeting_in_out_test,
                                                test_mast_in_out=test_mast_in_out_single_sgrna_cells, 
                                                qced_results_mast=qced_results_mast, 
                                                guide_to_cells=guide_to_single_guide_cells,
                                                guide_to_genes=guide_to_genes,
                                                cell_list=single_guide_cells)
names(targeting_in_out_single_sgrna_results) <- test_sgrnas
#Combine results into a matrix
targeting_in_out_single_sgrna_df <- rbind.fill.matrix(targeting_in_out_single_sgrna_results)
#Overwrite gene names
targeting_in_out_single_sgrna_df[,"gene"] <- target_genes_filter[targeting_in_out_single_sgrna_df[,"gene"], "external_gene_name"]
#Save RDS object
saveRDS(targeting_in_out_single_sgrna_df, file = paste0(out_dir, "targeting_in_out_single_sgrna_df.rds"))

#Run in-out test for all cells
targeting_in_out_all_results <- lapply(test_sgrnas, 
                                                FUN=run_targeting_in_out_test,
                                                test_mast_in_out=test_mast_in_out_all_cells, 
                                                qced_results_mast=qced_results_mast, 
                                                guide_to_cells=guide_to_cells,
                                                guide_to_genes=guide_to_genes,
                                                cell_list=all_cells)
names(targeting_in_out_all_results) <- test_sgrnas
#Combine results into a matrix
targeting_in_out_all_df <- rbind.fill.matrix(targeting_in_out_all_results)
#Overwrite gene names
targeting_in_out_all_df[,"gene"] <- target_genes_filter[targeting_in_out_all_df[,"gene"], "external_gene_name"]
#Save RDS object
saveRDS(targeting_in_out_all_df, file = paste0(out_dir, "targeting_in_out_all_df.rds"))

#Create a wrapper function so that we can apply the ntc tests to the targeting guides
run_targeting_ntc_test <- function(sgrna, test_mast_in_out, qced_results_mast, guide_to_cells, guide_to_genes, neg_control_sgrnas){
  qced_results_test <- qced_results_mast[unlist(guide_to_genes[sgrna]),]
  return(test_mast_in_out(sgrna, qced_results_test, guide_to_cells, neg_control_sgrnas))
}                                         

#Run ntc test for single-sgrna cells
targeting_ntc_single_sgrna_results <- lapply(test_sgrnas, 
                                                FUN=run_targeting_ntc_test,
                                                test_mast_in_out=test_mast_ntc_single_sgrna_cells, 
                                                qced_results_mast=qced_results_mast, 
                                                guide_to_cells=guide_to_single_guide_cells,
                                                guide_to_genes=guide_to_genes,
                                                neg_control_sgrnas=neg_control_sgrnas)
names(targeting_ntc_single_sgrna_results) <- test_sgrnas
#Combine results into a matrix
targeting_ntc_single_sgrna_df <- rbind.fill.matrix(targeting_ntc_single_sgrna_results)
#Overwrite gene names
targeting_ntc_single_sgrna_df[,"gene"] <- target_genes_filter[targeting_ntc_single_sgrna_df[,"gene"], "external_gene_name"]
#Save RDS object
saveRDS(targeting_ntc_single_sgrna_df, file = paste0(out_dir, "targeting_ntc_single_sgrna_df.rds"))

#Run ntc test for all cells
targeting_ntc_all_results <- lapply(test_sgrnas, 
                                             FUN=run_targeting_ntc_test,
                                             test_mast_in_out=test_mast_ntc_all_cells, 
                                             qced_results_mast=qced_results_mast, 
                                             guide_to_cells=guide_to_cells,
                                             guide_to_genes=guide_to_genes,
                                             neg_control_sgrnas=neg_control_sgrnas)
names(targeting_ntc_all_results) <- test_sgrnas
#Combine results into a matrix
targeting_ntc_all_df <- rbind.fill.matrix(targeting_ntc_all_results)
#Overwrite gene names
targeting_ntc_all_df[,"gene"] <- target_genes_filter[targeting_ntc_all_df[,"gene"], "external_gene_name"]
#Save RDS object
saveRDS(targeting_ntc_all_df, file = paste0(out_dir, "targeting_ntc_all_df.rds"))


# 8) Compare false discovery Test Rates ====

#Read in RDS objects
non_targeting_in_out_all_df <- readRDS(paste0(out_dir, "non_targeting_in_out_all_df.rds"))
non_targeting_in_out_single_sgrna_df <- readRDS(paste0(out_dir, "non_targeting_in_out_single_sgrna_df.rds"))
non_targeting_ntc_all_df <- readRDS(paste0(out_dir, "non_targeting_ntc_all_df.rds"))
non_targeting_ntc_single_sgrna_df <- readRDS(paste0(out_dir, "non_targeting_ntc_single_sgrna_df.rds"))

#Recast matrix test results to dataframes
non_targeting_in_out_all_df <- as.data.frame(non_targeting_in_out_all_df)
non_targeting_in_out_single_sgrna_df <- as.data.frame(non_targeting_in_out_single_sgrna_df)
non_targeting_ntc_all_df <- as.data.frame(non_targeting_ntc_all_df)
non_targeting_ntc_single_sgrna_df <- as.data.frame(non_targeting_ntc_single_sgrna_df)

#Recast p-value columns as numerics
non_targeting_in_out_all_df$p_val <- as.numeric(non_targeting_in_out_all_df$p_val)
non_targeting_in_out_single_sgrna_df$p_val <- as.numeric(non_targeting_in_out_single_sgrna_df$p_val)
non_targeting_ntc_all_df$p_val <- as.numeric(non_targeting_ntc_all_df$p_val)
non_targeting_ntc_single_sgrna_df$p_val <- as.numeric(non_targeting_ntc_single_sgrna_df$p_val)


