################################################################################

#cluster_signals_and_check_overlap.R

#The purpose of this script is twofold: 
# 1) Cluster the signals identified from multi-trait fine-mapping by weights
# 2) Identify the intersections between our screen target sites and signals

#This code runs with the R/4.2 module

################################################################################


# 0) Call libraries and set variables and directories ====

#Call libraries
library(scales)
library(ggtext)
library(glue)
library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap
library(reticulate)
library(umap)
library(igraph)
library(pbapply)

#Set directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/summary_files/"
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/"
plot_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/multi-trait_fine-mapping/clustering/"

#Set target bed file
target_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/target_selection/all_targeting_guides.bed"
#Set sceptre results file
sceptre_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/sceptre_results/discovery_results.pca_corrected.15_pcs.tsv"

#Set file_prefix for the filtering we want to use
#I did it this way because I wanted to trial multiple thresholds for cutting off signals at
file_prefix <- "purity_0.5.activity_0.5.gwas_5e-08" 

# 1) Read in files ====

#List files in input directory
input_files <- list.files(inp_dir)

#Separate file types
weight_files <- input_files[grepl("bmd_signal_weights", input_files)]
activity_files <- input_files[grepl("bmd_signal_activity", input_files)]
bed_files <- input_files[grepl("bmd_variants", input_files)]

#Get file prefixes
file_prefixes = str_replace(str_replace(weight_files, pattern="[.]tsv", replacement=""), pattern = "bmd_signal_weights[.]", replacement = "")

#Check number of signals in each file
count_bmd_signals <- function(temp_file_prefix, inp_dir){
  temp <- read.table(paste0(inp_dir, activity_files[file_prefixes == temp_file_prefix]), header = TRUE, row.names = 1)
  print(temp_file_prefix)
  print(nrow(temp))
}
lapply(file_prefixes, count_bmd_signals, inp_dir=inp_dir)

#Read in files and cluster
inp_raw <- read.table(paste0(inp_dir, activity_files[file_prefixes == file_prefix]), header = TRUE, row.names = 1)
#Remove completely empty columns
inp_filt <- inp_raw[,colSums(inp_raw) != 0] 
colnames(inp_filt) <- colnames(inp_filt) %>% toTitleCase() %>%
  str_replace_all(pattern = "_", replacement = " ") %>% 
  str_replace_all(pattern = "[.]", replacement = "-") %>%
  str_replace(pattern = "Alp", "ALP") %>% 
  str_replace(pattern = "Alt", "ALT") %>% 
  str_replace(pattern = "BF", "Bone Fracture") %>% 
  str_replace(pattern = "Bmi", "BMI") %>% 
  str_replace(pattern = "Ldl", "LDL") %>% 
  str_replace(pattern = "Hdl", "HDL") %>%
  str_replace(pattern = "Egfr Creat", "eGFR Creatinine") %>%
  str_replace(pattern = "Ggt", "GGT") %>% 
  str_replace(pattern = "Hba1c", "HbA1c") %>%
  str_replace(pattern = "Vitamin d", "Vitamin D") %>%
  str_replace(pattern = "Smoking Ever Never", "Smoking Ever/Never") %>%
  str_replace(pattern = "Lbs", "(lbs)") %>% 
  str_replace(pattern = "Bone Mineral Density", "BMD (Pan-UKBB)") %>%
  str_replace(pattern = "Whole body", replacement = "Whole-body")

#Output the number of signals per trait
temp <- cbind.data.frame(y = colSums(inp_filt > 0.5), name = colnames(inp_filt))
jpeg(paste0(plot_dir, file_prefix, ".signals_per_trait.jpeg"), width = 10800, height=6000, res=1000)
temp %>% ggplot(aes(x = name, y = y)) + geom_col() + 
  ylab("Count of Signals") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 16, angle=45, hjust = 1, vjust = 1), 
        legend.position = "none", axis.text.y = element_text(color="black", size=24), 
        axis.title.y = element_text(size = 24), axis.title.x = element_blank())
dev.off()
#Output the number of traits per signal
temp <- cbind.data.frame(y = rowSums(inp_filt > 0.5))
jpeg(paste0(plot_dir, file_prefix, ".traits_per_signal.trait.jpeg"), width = 10800, height=10800, res=1000)
temp %>% ggplot(aes(x = y)) + geom_histogram(binwidth = 1) + 
  xlab("Number of Non-BMD Traits") + ylab("Count of Signals") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=45, hjust = 1, vjust = 1), 
        legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title = element_text(size = 24))
dev.off()
#Calc number of BMD only-traits and traits with highest counts
sum(rowSums(inp_filt > 0.5) != 0)
sort(colSums(inp_filt > 0.5), decreasing = TRUE)[1:5]

#Clean up trait names
colnames(inp_filt) <- toTitleCase(colnames(inp_filt)) %>% str_replace_all(pattern = "_", replacement = " ") %>% 
  str_replace_all(pattern = "[.]", replacement = "-") %>%
  str_replace(pattern = "Alp", replacement = "ALP") %>% 
  str_replace(pattern = "Alt", replacement = "ALT") %>%
  str_replace(pattern = "BF", replacement = "Bone Fracture") %>% 
  str_replace(pattern = "Bmi", replacement = "BMI") %>% 
  str_replace(pattern = "Ldl", replacement = "LDL") %>% 
  str_replace(pattern = "Hdl", replacement = "HDL") %>%
  str_replace(pattern = "Egfr Creat", replacement = "eGFR Creatinine") %>%
  str_replace(pattern = "Ggt", replacement = "GGT") %>% 
  str_replace(pattern = "Hba1c", replacement = "HbA1c") %>%
  str_replace(pattern = "Vitamin d", replacement = "Vitamin D") %>%
  str_replace(pattern = "Smoking Ever Never", replacement = "Smoking Ever/Never") %>%
  str_replace(pattern = "Lbs", replacement = "(lbs)") %>% 
  str_replace(pattern = "Bone Mineral Density", replacement = "BMD (Pan-UKBB)") %>%
  str_replace(pattern = "Whole body", replacement = "Whole-body")

#Make heatmaps
jpeg(paste0(plot_dir, file_prefix, ".heatmap.jpeg"), width = 7200, height = 10800, res = 1000)
pheatmap(t(inp_filt), 
         show_colnames = FALSE, color=colorRampPalette(c("navy", "white", "red"),)(50), fontsize_row = 16,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
dev.off()

pca_result <- prcomp(inp_filt, scale. = TRUE)
#Plot pcas
jpeg(paste0(plot_dir, file_prefix, ".pca_plot.jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
plot(pca_result$x[, 1], pca_result$x[, 2], col = "blue", pch = 16, main = "Principal Components Plot", xlab = "PC1", ylab = "PC2")
# Add grid lines (optional)
abline(h = 0, v = 0, col = "gray", lty = 2)
dev.off()

# Example data (replace this with your actual data)
set.seed(5)
# Perform UMAP transformation
umap_result <- umap(inp_filt, n_neighbors = 15, n_components = 2)
# Plot the UMAP results with KNN clusters
plot_umap_colored_trait <- function(trait, umap_result, inp_filt, plot_dir, file_prefix){
  file_trait <- trait %>% str_replace_all(pattern = " ", replacement = "_") %>% str_replace_all(pattern = "/", replacement = "_")
  trait_color = inp_filt %>% mutate(trait_col = ifelse(inp_filt[,trait] > 0.5, 1, 0)) %>% select(trait_col)
  jpeg(paste0(plot_dir, "umap.", file_trait, ".", file_prefix, ".jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
    plot(umap_result$layout, col = as.numeric(trait_color$trait_col) + 1, pch = 16, main = paste0(trait, " UMAP"), xlab = "UMAP 1", ylab = "UMAP 2")
  dev.off()
}
lapply(colnames(inp_filt), plot_umap_colored_trait, umap_result=umap_result, inp_filt=inp_filt, plot_dir=plot_dir, file_prefix=file_prefix)

#Output a supplemental table
write.table(inp_filt, file=paste0(inp_dir, "supplement.multi-trait_results.", file_prefix, ".tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


# 3) Check for signal intersections with the file of targets ====

#Read in the target file
target_bed <- read.csv(target_file, header = FALSE, sep = "\t")
colnames(target_bed) <- c("chr", "start", "end", "target")
#Filter out the positive controls
target_bed <- target_bed %>% filter(str_count(target, "rs") > 0)

#Read in the bed file of signals
bed_raw <- read.table(paste0(inp_dir, bed_files[file_prefixes == file_prefix]), header = TRUE, sep = "\t", row.names = 1)

#Intersect the bed files to see if any fine-mapped variants were found in the repressed targeted regions
intersect_target_bed <- function(target_bed_row, bed_raw, repression_range=1000){
  #Get position of repression range
  target_chr = as.numeric(str_replace(target_bed_row[1], "chr", ""))
  target_start = as.numeric(target_bed_row[2])
  target_end = as.numeric(target_bed_row[3])
  target_name = target_bed_row[4]
  #Filter for SNPs in target repression range
  temp <- bed_raw[which(bed_raw$chr == target_chr & target_start - repression_range <= bed_raw$start & target_end + repression_range >= bed_raw$end),]
  #Check if any fine-mapped variants were found in the repressed region
  if (nrow(temp) == 0) {
    return(NULL)
  } else {
    temp2 <- suppressWarnings(cbind.data.frame(matrix(unlist(rep(target_bed_row, nrow(temp))), ncol = length(target_bed_row), byrow = TRUE), temp))
    return(temp2)
  }
}
#Apply the function over the targeting sgRNAs and then remove any non-intersections before forming into a dataframe
intersection_list <- pbapply(target_bed, FUN = intersect_target_bed, bed_raw=bed_raw, MARGIN = 1)
intersection_list <- intersection_list[lapply(intersection_list, FUN = is.null) == FALSE]
targeted_signal_df <- bind_rows(intersection_list)
colnames(targeted_signal_df) <- c("target_chr", "target_start", "target_end", "target_name", colnames(bed_raw))
#Create a grna group column in the data frame of targeted signals to group the targets by
targeted_signal_df <- targeted_signal_df %>% mutate(grna_group = str_replace(target_name, "_.+$", ""))
#Also append a signal column
targeted_signal_df <- targeted_signal_df %>% mutate(signal = str_replace(signal.snp_id, "((\\.[^\\.]+){3}\\.).*", "\\1")) %>% mutate(signal = str_replace(signal, "\\.$", ""))
#Identify the number of unique signals that were targeted in the screen and how many targets that lines up with
targeted_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow #41 targets
targeted_signal_df %>% dplyr::select(signal) %>% unique %>% nrow #42 signals

#Read in the file of sceptre results to identify how many signals had genes identified if any
sceptre_raw <- read.csv(sceptre_file, sep = "\t", header = TRUE)
#Merge the dataframes by the grna group names
targeted_perturbed_signal_df <- targeted_signal_df %>% inner_join(y = sceptre_raw, by = "grna_group")
#Identify the number of unique signals that were actually perturbed in the screen and how many targets that lines up with
targeted_perturbed_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow #6 targets
targeted_perturbed_signal_df %>% dplyr::select(signal) %>% unique %>% nrow #6 signals
#Also count the number of genes
targeted_perturbed_signal_df %>% dplyr::select(response_id) %>% unique %>% nrow #8 genes could be involved

