################################################################################

#test_markers_for_screen_hits.R

#The purpose of this script is to test the osteoblast marker genes for 
#perturbations with the sgRNA that resulted in one or more gene knockdowns in 
#the CRISPRi screen. To do this we use the the SCEPTRE methodology. 
#The script is set to run with the following modules:

#  R/4.2.3

################################################################################

# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(dplyr)
library(cowplot)
library(katlabutils)
library(sceptre) #v0.3.0
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
library(ggrepel)

#Set directories and file locations
qced_results_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/quality_control/aggr/aggr.post_qc.h5ad"
guides_per_cell_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/protospacer_calls_per_cell.csv"
sgrna_target_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/all_targeting_guides.bed"
gencode_genes_loc <- "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.annotation.gene_only.bed"
marker_genes_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/quality_control/hFOB_marker_genes.csv"
screen_hits_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/differential_expression/sceptre_single_sgrna/discovery_results.tsv" 
lentiviral_pool_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/lentiviral_rna_seq_raw_read_counts.tsv"
guide_sequence_ref_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/cellranger_outputs/aggr/crispr_analysis/feature_reference.csv"
vep_loc <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/target.vep_annotations.txt"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/pilot_hFOB_screen_bone/differential_expression/marker_gene_test/"

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
rownames(marker_genes_raw) <- marker_genes_raw$ENSG

#Read in screen hit results
screen_hits_raw <- read.table(screen_hits_loc, sep = "\t", header = TRUE)
#Remove RAB1A control
screen_hits_filt <- screen_hits_raw[screen_hits_raw$response_id != "RAB1A",]

#Read in guide sequence reference
guide_sequence_ref_raw <- read.csv(guide_sequence_ref_loc, header = TRUE, stringsAsFactors = FALSE, sep = ",")

#Read in file of lentiviral read counts
lentiviral_pool_raw <- read.csv(lentiviral_pool_loc, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Read in vep annotations file
vep_raw <- read.csv(vep_loc, sep = "\t", header = TRUE)

# 2) Get List of Marker Genes to Test and Pair with Perturbed Targets ====

#Get list of unique sgrna targets
sgrnas <- unique(guides_per_cell_raw$feature_call)
#Divy up the guide rnas into their respective categories
pos_control_sgrnas <- sgrnas[1:2]
neg_control_sgrnas <- sgrnas[3:29]
test_sgrnas <- sgrnas[30:length(sgrnas)]

#Get unique list of perturbed sgRNA groups
hit_grna_groups <- unique(screen_hits_filt$grna_group) 
#Combine into a dataframe
response_grna_group_pairs_marker <- cbind.data.frame(response_id=rep(marker_genes_raw$ENSG, length(hit_grna_groups)),
                                                     grna_group=rep(hit_grna_groups, each=nrow(marker_genes_raw)))

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
grna_lowmoi_marker <- as.matrix(t(rbind.fill.matrix(grna_list)))
rownames(grna_lowmoi_marker) <- sgrnas
#Filter for cells with only 1 guide
cells_with_single_guide <- qced_cells[colSums(grna_lowmoi_marker != 0) == 1]
grna_lowmoi_marker <- grna_lowmoi_marker[,colSums(grna_lowmoi_marker != 0) == 1]

### Make response-by-cell matrix ###
sce <- as.SingleCellExperiment(qced_results_seurat)
counts_matrix <- as.matrix(assay(sce, "counts"))
response_lowmoi_marker <- counts_matrix[,cells_with_single_guide]
colnames(response_lowmoi_marker) <- NULL

### Make Cell Covariate Matrix ###
covariate_lowmoi_marker <- qced_results_seurat@meta.data
covariate_lowmoi_marker <- covariate_lowmoi_marker[cells_with_single_guide,c("total_counts", "n_genes_by_counts", "Pool", "pct_counts_mt")]
rownames(covariate_lowmoi_marker) <- NULL
colnames(covariate_lowmoi_marker) <- c("response_n_umis", "response_n_nonzero", "bio_rep", "p_mito")
covariate_lowmoi_marker$bio_rep <- as.factor(covariate_lowmoi_marker$bio_rep)
covariate_lowmoi_marker$p_mito <- covariate_lowmoi_marker$p_mito/100 #rescale percentages to decimal
#Append on number of gRNA UMIs per cells
covariate_lowmoi_marker[,"grna_n_umis"] <- colSums(grna_lowmoi_marker)
#Calculate PCAs on Seurat Object and add top 50 to covariate matrix (We can use fewer later)
qced_results_seurat <- FindVariableFeatures(qced_results_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(qced_results_seurat)
qced_results_seurat <- ScaleData(qced_results_seurat, features = all.genes)
qced_results_seurat <- RunPCA(qced_results_seurat, rev.pca=TRUE)
qced_results_seurat <- ProjectDim(qced_results_seurat, reduction="pca")
temp <- Embeddings(qced_results_seurat, reduction="pca")[cells_with_single_guide,]
covariate_lowmoi_marker <- cbind.data.frame(covariate_lowmoi_marker, temp[,1:50])

### Make gRNA group table ###
grna_group_lowmoi_marker <- rbind.data.frame(cbind.data.frame(grna_id=test_sgrnas, grna_group=str_replace(test_sgrnas, "_[0-9]+","")),
                                               cbind.data.frame(grna_id=pos_control_sgrnas, grna_group=pos_control_sgrnas),
                                               cbind.data.frame(grna_id=neg_control_sgrnas, grna_group=rep("non-targeting", length(neg_control_sgrnas))))

# 4) Set the Formula Object ====

#Set the formula object for use with pcas in the cis analysis
formula_object <- formula(~log(response_n_umis) + 
                            log(response_n_nonzero) +
                            bio_rep + 
                            p_mito + 
                            PC_1 + PC_2 + PC_3)

# 5) Run and Assess Calibration Check ====

#Run basic cis check
calibration_result <- run_sceptre_lowmoi(
  response_matrix = response_lowmoi_marker,
  grna_matrix = grna_lowmoi_marker,
  covariate_data_frame = covariate_lowmoi_marker,
  grna_group_data_frame = grna_group_lowmoi_marker,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_marker,
  calibration_check = TRUE # calibration_check TRUE
) 

#Plot cis calibration result
jpeg(paste0(out_dir, "sceptre_calibration.marker_test.jpeg"), width = 10800, height=10800, res=1000)
plot_calibration_result(calibration_result)
dev.off()
#Write cis calibration results to file
write.table(calibration_result, file = paste0(out_dir, "non_targeting_test_results.marker_test.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 6) Run Discovery Analyses ====

### Run Discovery on single-guide cells ###
discovery_result <- run_sceptre_lowmoi(
  response_matrix = response_lowmoi_marker,
  grna_matrix = grna_lowmoi_marker,
  covariate_data_frame = covariate_lowmoi_marker,
  grna_group_data_frame = grna_group_lowmoi_marker,
  formula_object = formula_object,
  response_grna_group_pairs = response_grna_group_pairs_marker,
  calibration_check = FALSE
) 

#Compare calibration and discovery results
jpeg(paste0(out_dir, "sceptre_discovery_calibration_compare.marker_test.jpeg"), width = 10800, height=10800, res=1000)
compare_calibration_and_discovery_results(
  calibration_result = calibration_result,
  discovery_result = discovery_result) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_blank(),
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16), axis.title=element_text(size=20), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=15)) 
dev.off()

#Make volcano
jpeg(paste0(out_dir, "discovery_volcano.marker_test.jpeg"), width = 10800, height=10800, res=1000)
make_volcano_plot(discovery_result = discovery_result) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_blank(),
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=16), axis.title=element_text(size=20), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=16), 
        legend.text = element_text(size=15)) 
dev.off()

#Append BH-Adjusted p-values to the output
discovery_result <- discovery_result %>% mutate(p_value_BH=p.adjust(p_value, method = "BH"))
#obtain discovery set and write to file
discovery_set <- obtain_discovery_set(discovery_result)
discovery_set$response_id <- marker_genes_raw[discovery_set$response_id,"Symbol"]
write.table(discovery_set, file = paste0(out_dir, "discovery_set.marker_test.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(discovery_result, file = paste0(out_dir, "discovery_results.marker_test.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#Count discovery test results
discovery_result %>% dplyr::filter(is.na(p_value)) %>% nrow #20 test had insufficient cells to run (all the BGLAP tests)
discovery_result %>% dplyr::filter(is.na(p_value) == FALSE) %>% nrow #160 tests ran (all 8 other marker genes) 
discovery_result %>% dplyr::filter(is.na(p_value) == FALSE & log_2_fold_change < 0) %>% nrow #101/160 (63.1%) of test showed marker gene repression
#And 59/160 (36.9%) tests showed increased marker gene expression

# 7) Make publication grade figures ====
#Much of this is reproduced code from sceptre v.0.9.0

#Set theme function
get_my_theme <- function(element_text_size = 36) {
  ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(color = "black"),
                                       axis.text.x = ggplot2::element_text(color="black", size=element_text_size*.75),
                                       axis.text.y = ggplot2::element_text(color="black", size=element_text_size*.75),
                                       axis.title = ggplot2::element_text(color = "black", size = element_text_size),
                                       panel.grid.major = ggplot2::element_blank(),
                                       panel.grid.minor = ggplot2::element_blank(),
                                       panel.border = ggplot2::element_blank(),
                                       panel.background = ggplot2::element_blank(),
                                       plot.title = ggplot2::element_blank())
}

### Make the QQ-Plot First ###
#Set code for qqplot bands
StatQQBand <- ggplot2::ggproto("StatQQBand", ggplot2::Stat,
                               required_aes = c("y"), # the sample will come through the y aesthetic
                               # code below implements logic that you only want to have one band per panel,
                               # rather than having one band per group
                               setup_data = function(data, params) {
                                 if ("group" %in% names(data)) {
                                   # calculate the maximum across panels of the number of unique numbers of
                                   # points there are in each group
                                   max_unique_pts_per_group <- data |>
                                     dplyr::group_by(PANEL, group) |>
                                     dplyr::summarise(pts_per_group = dplyr::n()) |>
                                     dplyr::summarise(unique_pts_per_group = length(unique(pts_per_group))) |>
                                     dplyr::summarise(max(unique_pts_per_group)) |>
                                     dplyr::pull()
                                   if (max_unique_pts_per_group > 1) {
                                     # throw an error if there is a panel with differing numbers of points per group
                                     stop("Within each panel, you must have the same number of points per QQ plot!")
                                   } else {
                                     data |>
                                       dplyr::select(y, PANEL, group) |> # remove attributes like color, shape, etc.
                                       dplyr::filter(group == min(group)) # keep only one of the groups to plot
                                   }
                                 } else {
                                   data
                                 }
                               },
                               compute_group = function(data, scales, distribution = "unif", max_pts_to_plot = 500, ci_level = 0.95) {
                                 # get the quantile function from the stats package
                                 quantile_fun <- eval(parse(text = sprintf("stats::q%s", distribution)))
                                 # set x and y transformations to identity if they are not given
                                 if (is.null(scales$x$trans)) {
                                   scales$x$trans <- scales::identity_trans()
                                 }
                                 if (is.null(scales$y$trans)) {
                                   scales$y$trans <- scales::identity_trans()
                                 }
                                 # compute the upper and lower confidence bands
                                 data |>
                                   # the given y is already transformed, so transform it back first
                                   dplyr::mutate(y = scales$y$trans$inverse(y)) |>
                                   # apply the formulas for the confidence bands
                                   dplyr::mutate(
                                     r = rank(y, ties.method = "first"),
                                     x = quantile_fun(stats::ppoints(dplyr::n())[r]),
                                     ymin = quantile_fun(stats::qbeta(p = (1 - ci_level) / 2, shape1 = r, shape2 = dplyr::n() + 1 - r)),
                                     ymax = quantile_fun(stats::qbeta(p = (1 + ci_level) / 2, shape1 = r, shape2 = dplyr::n() + 1 - r))
                                   ) |>
                                   # transform
                                   dplyr::mutate(
                                     x = scales$x$trans$transform(x),
                                     y = scales$y$trans$transform(y),
                                     ymin = scales$y$trans$transform(ymin),
                                     ymax = scales$y$trans$transform(ymax)
                                   ) |>   # downsample
                                   dplyr::mutate(x_int = ggplot2::cut_interval(x, n = max_pts_to_plot),
                                                 y_int = ggplot2::cut_interval(y, n = max_pts_to_plot)) |>
                                   dplyr::group_by(x_int, y_int) |>
                                   dplyr::slice_sample(n = 1) |>
                                   dplyr::ungroup()
                               }
)
stat_qq_band <- function(mapping = NULL, data = NULL, geom = "ribbon",
                         position = "identity", show.legend = FALSE,
                         inherit.aes = TRUE, distribution = "unif",
                         max_pts_to_plot = 500,
                         ci_level = 0.95, ...) {
  ggplot2::layer(
    stat = StatQQBand, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(alpha = 0.25, distribution = distribution,
                  max_pts_to_plot = max_pts_to_plot, ci_level = ci_level, ...)
  )
}

#Set publication-grade figure making function and call it
make_publication_qq <- function(calibration_result, discovery_result, p_thresh = max(discovery_set$p_value),
                                                      transform_scale = TRUE, include_legend = TRUE,
                                                      include_y_axis_text = TRUE, point_size = 4,
                                                      transparency = 0.8) {
  lab <- c(rep(factor("Negative control"), nrow(calibration_result)),
           rep(factor("Target sgRNA"), nrow(discovery_result))) |>
    factor(levels = c("Target sgRNA", "Negative control"))
  df <- data.frame(p_value = c(calibration_result$p_value, discovery_result$p_value), lab = lab)
  
  p_out <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(y = p_value, col = lab)) +
    stat_qq_points(size = point_size, alpha = transparency) +
    stat_qq_band(data = df[df$lab == "Negative control",], mapping = ggplot2::aes(y = p_value, col = lab)) +
    ggplot2::labs(x = "Expected Null P-Value", y = "Observed P-Value") +
    ggplot2::geom_abline(col = "black") +
    get_my_theme() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = if (include_legend) "bottom" else "none",
                   axis.title.y = if (include_y_axis_text) ggplot2::element_text() else ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("dodgerblue3", "firebrick2"))
  
  if (!transform_scale) {
    p_out <- p_out +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_reverse() +
      ggplot2::ggtitle("QQ plot (bulk)")
  } else {
    p_out <- p_out +
      ggplot2::scale_x_continuous(trans = revlog_trans(10), breaks = c(0.001,0.01,0.1,1)) +
      ggplot2::scale_y_continuous(trans = revlog_trans(10), breaks = c(0.001,0.01,0.1,1)) +
      ggplot2::ggtitle("QQ plot (tail)") +
      (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL)
  }
  
  if (include_legend) {
    p_out <- p_out +
      ggplot2::theme(legend.position = c(0.8, 0.05),
                     legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                     legend.title = ggplot2::element_blank(),
                     legend.text = element_text(size = 24)) +
      ggplot2::guides(color = ggplot2::guide_legend(
        keywidth = 0.0,
        keyheight = 0.3,
        default.unit = "inch",
        override.aes = list(size = 5)))
  }
  
  return(p_out)
}
jpeg(paste0(out_dir, "dicovery_qq.publication_grade.marker_test.jpeg"), width = 10800, height=10800, res=1000)
make_publication_qq(calibration_result = calibration_result,
                    discovery_result = discovery_result)
dev.off()

### Make the Volcano Plot ###
get_screen_hit_genes <- function(grna_group_check, screen_hits_filt){
  return(paste0(screen_hits_filt[screen_hits_filt$grna_group == grna_group_check, "response_id"],collapse = ","))
}
make_volcano_plot <- function(discovery_result, p_thresh, x_limits = c(-1.5, 1.5), transparency = 0.5, point_size = 4) {
  p_lower_lim <- 1e-20
  temp_df <- discovery_result |> dplyr::mutate(reject = p_value_BH <= p_thresh,
                                               p_value_BH = ifelse(p_value_BH < p_lower_lim, p_lower_lim, p_value_BH),
                                               log_2_fold_change = ifelse(log_2_fold_change > x_limits[2], x_limits[2], log_2_fold_change),
                                               log_2_fold_change = ifelse(log_2_fold_change < x_limits[1], x_limits[1], log_2_fold_change),
                                               gene_lab = ifelse(p_value_BH <= p_thresh & !(is.na(p_value_BH)), marker_genes_raw[response_id,"Symbol"],"")) |>
    dplyr::mutate(gene_lab = ifelse(gene_lab %in% pos_control_sgrnas, "", gene_lab))
  hit_genes <- sapply(temp_df$grna_group, FUN = get_screen_hit_genes, screen_hits_filt = screen_hits_filt)
  temp_df <- temp_df %>% dplyr::mutate(hit_genes=hit_genes) %>%
    dplyr::mutate(gene_lab = ifelse(gene_lab != "", paste0(hit_genes, "/", gene_lab), ""))
  out <- ggplot2::ggplot(data = temp_df,
                         mapping = ggplot2::aes(x = log_2_fold_change, y = p_value_BH, col = reject)) +
    ggplot2::geom_point(alpha = transparency, size = point_size) +
    ggrepel::geom_text_repel(mapping = ggplot2::aes(label = gene_lab), color = "dodgerblue3", box.padding = 0.5, max.overlaps = Inf, size = 7.5, force = 5) + 
    ggplot2::scale_y_continuous(trans = revlog_trans(10), expand = c(0.02, 0), breaks = c(1, 0.1)) +
    get_my_theme() + ggplot2::xlab("Log Fold Change") + ggplot2::ylab("Adj. P-Value") +
    (if (!is.na(p_thresh)) ggplot2::geom_hline(yintercept = p_thresh, linetype = "dashed") else NULL) +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_color_manual(values = c("gray", "dodgerblue3")) +
    ggplot2::theme(plot.margin = margin(0,0.4,0,0, "cm")) + ggplot2::ggtitle("Discovery volcano plot")
  return(out)
}

#Call Volcano function
jpeg(paste0(out_dir, "dicovery_volcano.publication_grade.marker_test.jpeg"), width = 10800, height=10800, res=1000)
make_volcano_plot(discovery_result = discovery_result, p_thresh = 0.1)
dev.off()
