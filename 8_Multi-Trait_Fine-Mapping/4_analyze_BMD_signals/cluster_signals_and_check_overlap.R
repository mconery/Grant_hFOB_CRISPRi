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
library(tools)

#Set directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/summary_files/"
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/"
plot_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/multi-trait_fine-mapping/"

#Set target bed file
target_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/target_selection/all_targeting_guides.bed"
#Set sceptre results file
sceptre_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/sceptre_results/discovery_results.pca_corrected.15_pcs.tsv"

#Set file_prefix for the filtering we want to use
file_prefix <- "purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered" 

# 1) Read in files ====

#List files in input directory
input_files <- list.files(inp_dir)

#Separate file types
weight_files <- input_files[grepl("bmd_signal_weights", input_files)]
activity_files <- input_files[grepl("bmd_signal_activity", input_files)]
signal_files <- input_files[grepl("bmd_signals", input_files)]
bed_files <- input_files[grepl("bmd_variants", input_files)]

#Get file prefixes
file_prefixes = str_replace(str_replace(signal_files, pattern="[.]tsv", replacement=""), pattern = "bmd_signals[.]", replacement = "")

#Make a function that corrects the names for the traits
correct_trait_names <- function(trait_names){
  trait_names %>% str_replace_all(pattern = "_", replacement = " ") %>% 
    str_replace_all(pattern = "[.]", replacement = "-") %>%
    toTitleCase %>%
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
    str_replace(pattern = "Whole Body", replacement = "Whole-Body") %>% 
    return
}
#Make a function that checks the distributions of non-BMD traits per signal and signals per non-BMD Trait
count_bmd_signals <- function(temp_file_prefix, inp_dir, out_plot_dir=plot_dir){
  temp <- read.table(paste0(inp_dir, signal_files[file_prefixes == temp_file_prefix]), header = TRUE, row.names = 1)
  print(temp_file_prefix)
  print(nrow(temp))
  #Calc number of BMD only-traits and traits with highest counts
  temp %>% dplyr::select(num_traits) %>% dplyr::filter(num_traits > 0) %>% nrow %>% print
  temp %>% dplyr::select(num_traits) %>% dplyr::filter(num_traits == 0) %>% nrow %>% print
  temp %>% dplyr::select(num_traits) %>% dplyr::filter(num_traits == 1) %>% nrow %>% print
  temp2 <- temp %>% dplyr::select(traits) %>% dplyr::filter(traits != "None") %>% as.vector() %>% unlist %>% str_split(pattern = ",") %>%
    unlist %>% table %>% sort(decreasing = TRUE)
  print(temp2[1:5])
  #Output the number of signals per trait
  temp3 <- cbind.data.frame(y = temp2, name = names(temp2))
  temp3$name <- correct_trait_names(temp3$name)
  t_plot <- temp3 %>% ggplot(aes(x = name, y = y.Freq)) + geom_col() + 
    ylab("Count of Signals") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 16, angle=45, hjust = 1, vjust = 1), 
          legend.position = "none", axis.text.y = element_text(color="black", size=24), 
          axis.title.y = element_text(size = 24), axis.title.x = element_blank())
  tiff(paste0(out_plot_dir, temp_file_prefix, ".signals_per_trait.tiff"), width = 10800, height=6000, res=1000)
    print(t_plot)
  dev.off()
  #Output the number of traits per signal
  temp3 <- cbind.data.frame(y = temp$num_traits)
  #Make plot
  t_plot <- temp3 %>% ggplot(aes(x = y)) + geom_histogram(binwidth = 1) + 
    xlab("Number of Non-BMD Traits") + ylab("Count of Signals") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 24, angle=45, hjust = 1, vjust = 1), 
          legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title = element_text(size = 24))
  tiff(paste0(out_plot_dir, temp_file_prefix, ".traits_per_signal.trait.tiff"), width = 10800, height=10800, res=1000)
    print(t_plot)
  dev.off()
}
#Call function
hold = lapply(file_prefixes, count_bmd_signals, inp_dir=inp_dir)
remove(hold)

#Make a function to cluster and make plots for the various activity files
cluster_by_activity <- function(activity_file, out_plot_dir=paste0(plot_dir, "clustering/"), input_dir=inp_dir, activity_thresh=0.95){
  out_file_prefix <- activity_file %>% str_replace(pattern = "[.]tsv", replacement = "") %>% str_replace(pattern = "bmd_signal_activity[.]", replacement = "")
  out_file_folder <- paste0("clustering.", 
                            ifelse(str_detect(out_file_prefix, "highest_residual_filtered"), "highest_residual_filtered.", ""),
                            ifelse(str_detect(out_file_prefix, "binarized"), "binarized/", "/"))
  #Read in files and cluster
  inp_raw <- read.table(paste0(input_dir, activity_file), header = TRUE, row.names = 1)
  #Remove completely empty columns
  inp_filt <- inp_raw[,colSums(inp_raw) != 0] 
  #Rename column names
  colnames(inp_filt) <- colnames(inp_filt) %>% correct_trait_names
  #Create needed directory
  dir.create(paste0(out_plot_dir, out_file_folder), recursive = TRUE, mode = "0777", showWarnings = FALSE)
  #Make heatmaps
  tiff(paste0(out_plot_dir, out_file_folder, out_file_prefix, ".heatmap.tiff"), width = 7200, height = 10800, res = 1000)
  print(pheatmap(t(inp_filt), 
           show_colnames = FALSE, color=colorRampPalette(c("#C5C6D0", "black"),)(50), fontsize_row = 16,
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean"))
  dev.off()
  
  # Set seed
  set.seed(5)
  # Perform UMAP transformation
  umap_result <- umap(inp_filt, n_neighbors = 15, n_components = 2)
  #Plot uncolored UMAP clusters
  tiff(paste0(out_plot_dir, out_file_folder, "umap.uncolored.tiff"), width = 6, height = 11.5, units = 'in', res = 500)
  plot(umap_result$layout, pch = 16, main = "", xlab = "UMAP 1", ylab = "UMAP 2")
  dev.off()
  # Plot the UMAP results with KNN clusters colored by trait
  plot_umap_colored_trait <- function(trait, umap_result, inp_filt, plot_dir, file_prefix){
    file_trait <- trait %>% str_replace_all(pattern = " ", replacement = "_") %>% str_replace_all(pattern = "/", replacement = "_")
    trait_color = inp_filt %>% mutate(trait_col = ifelse(inp_filt[,trait] > 0.95, 1, 0)) %>% select(trait_col)
    tiff(paste0(plot_dir, "umap.", file_trait, ".tiff"), width = 6, height = 11.5, units = 'in', res = 1000)
    plot(umap_result$layout, col = as.numeric(trait_color$trait_col) + 1, pch = 16, main = paste0(trait, " UMAP"), xlab = "UMAP 1", ylab = "UMAP 2")
    dev.off()
  }
  lapply(colnames(inp_filt), plot_umap_colored_trait, umap_result=umap_result, inp_filt=inp_filt, plot_dir=paste0(out_plot_dir, out_file_folder), file_prefix=out_file_prefix)
  
  #Output a supplemental table
  write.table(inp_filt, file=paste0(inp_dir, "supplement.multi-trait_results.", out_file_prefix, ".tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  
}
#Call function
lapply(c("bmd_signal_activity.purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered.binarized.tsv",
         "bmd_signal_activity.purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered.tsv"), cluster_by_activity)

#Make a function to cluster and make plots for the various weight files
cluster_by_weight <- function(weight_file, out_plot_dir=paste0(plot_dir, "clustering/"), input_dir=inp_dir, activity_thresh=0.95){
  out_file_prefix <- weight_file %>% str_replace(pattern = "[.]tsv", replacement = "") %>% str_replace(pattern = "bmd_signal_weights[.]", replacement = "")
  out_file_folder <- paste0("clustering.", 
                            ifelse(str_detect(out_file_prefix, "highest_residual_filtered"), "highest_residual_filtered.", ""),
                            "binarized/")
  #Read in files and cluster
  inp_raw <- read.table(paste0(input_dir, weight_file), header = TRUE, row.names = 1)
  binary_activity_raw <- read.table(paste0(input_dir, "bmd_signal_activity.", out_file_prefix, ".binarized.tsv"), header = TRUE, row.names = 1)
  #Binarize the weight file
  weight_binarized <- as.data.frame(ifelse(binary_activity_raw == 0,0,1) * ifelse(inp_raw > 0, 1, -1))
  #Remove completely empty columns
  weight_filt <- weight_binarized[,colSums(weight_binarized != 0) != 0] 
  #Rename column names
  colnames(weight_filt) <- colnames(weight_filt) %>% correct_trait_names
  #Make heatmaps
  tiff(paste0(out_plot_dir, out_file_folder, out_file_prefix, ".weight.heatmap.tiff"), width = 7200, height = 10800, res = 1000)
  print(pheatmap(t(weight_filt), 
                 show_colnames = FALSE, color=colorRampPalette(c("navy", "#C5C6D0", "red"),)(50), fontsize_row = 16,
                 clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean"))
  dev.off()
  
  # Set seed
  set.seed(5)
  # Perform UMAP transformation
  umap_result <- umap(weight_filt, n_neighbors = 15, n_components = 2)
  #Plot uncolored UMAP clusters
  tiff(paste0(out_plot_dir, out_file_folder, "umap.weight.tiff"), width = 6, height = 11.5, units = 'in', res = 500)
  plot(umap_result$layout, pch = 16, main = "", xlab = "UMAP 1", ylab = "UMAP 2")
  dev.off()
  
  #Make a barplot of the signals by trait split by positive and negative effects
  t_plot <- weight_filt %>% mutate(signal=rownames(weight_filt)) %>% pivot_longer(!signal, names_to = "trait", values_to = "activity") %>%
    group_by(trait, activity) %>% summarize(count=n()) %>% filter(activity != 0) %>% mutate(count=activity*count) %>%
    ggplot(aes(x = trait, y = count, fill = activity)) + geom_bar(position = "stack", stat="identity") + 
    ylab("Count of Signals") + 
    scale_fill_gradientn(colours = colorRampPalette(c("navy", "red"),)(100)) + 
    scale_y_continuous(breaks = c(-20,-10,0,10,20), labels = c(20,10,0,10,20)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 16, angle=45, hjust = 1, vjust = 1), 
          legend.position = "none", axis.text.y = element_text(color="black", size=16), 
          axis.title.y = element_text(size = 16), axis.title.x = element_blank())
  tiff(paste0(out_plot_dir, out_file_folder, out_file_prefix, ".signals_per_trait.weighted.tiff"), width = 10800, height=6000, res=1000)
  print(t_plot)
  dev.off()
  
  #Make the barplot again, but this time make it vertical
  t_plot <- weight_filt %>% mutate(signal=rownames(weight_filt)) %>% pivot_longer(!signal, names_to = "trait", values_to = "activity") %>%
    group_by(trait, activity) %>% summarize(count=n()) %>% filter(activity != 0) %>% mutate(count=activity*count) %>%
    ggplot(aes(y = trait, x = count, fill = activity)) + geom_bar(position = "stack", stat="identity") + 
    xlab("Count of Signals") + 
    scale_fill_gradientn(colours = colorRampPalette(c("navy", "red"),)(100)) + 
    scale_x_continuous(breaks = c(-20,-10,0,10,20), labels = c(20,10,0,10,20)) + 
    scale_y_discrete(limits=rev) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 16, hjust = 1), 
          legend.position = "none", axis.text.x = element_text(color="black", size=16), 
          axis.title.x = element_text(size = 16), axis.title.y = element_blank())
  tiff(paste0(out_plot_dir, out_file_folder, out_file_prefix, ".signals_per_trait.weighted.vertical.tiff"), height = 10800, width=6000, res=1000)
  print(t_plot)
  dev.off()
  
  #Write table to file
  write.table(weight_filt, file=paste0(inp_dir, "supplement.multi-trait_results.", out_file_prefix, ".weights.tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
}
lapply(c("bmd_signal_weights.purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered.tsv")
       , cluster_by_weight)

# 2) Read in Signal Files and Process for Supplement ====

#Read in the signal-level results file and the bed file of varaint-level info about signals
bed_raw <- read.table(paste0(inp_dir, bed_files[file_prefixes == file_prefix]), header = TRUE, sep = "\t", row.names = 1)
signal_raw <- read.table(paste0(inp_dir, signal_files[file_prefixes == file_prefix]), header = TRUE, sep = "\t", row.names = 1)

#Clean up trait names
bed_raw$active_traits <- bed_raw$active_traits %>% str_replace_all(pattern = ",", replacement = ", ") %>% correct_trait_names
signal_raw$traits <- signal_raw$traits %>% str_replace_all(pattern = ",", replacement = ", ") %>% correct_trait_names

#Write files
write.table(signal_raw, file = paste0(inp_dir, "supplement.signals.", file_prefix, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(bed_raw, file = paste0(inp_dir, "supplement.variants.", file_prefix, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Read in unfiltered signals files
unfiltered_signal_raw <- read.table(paste0(inp_dir, "bmd_signals.purity_0.0.activity_0.5.gwas_1.0.tsv"), header = TRUE, sep = "\t", row.names = 1)
unfiltered_signal_raw$traits <- unfiltered_signal_raw$traits %>% str_replace_all(pattern = ",", replacement = ", ") %>% correct_trait_names
write.table(unfiltered_signal_raw, file = paste0(inp_dir, "supplement.signals.purity_0.0.activity_0.5.gwas_1.0.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3) Check for signal intersections with the file of targets ====

#Read in the target file
target_bed <- read.csv(target_file, header = FALSE, sep = "\t")
colnames(target_bed) <- c("chr", "start", "end", "target")
#Filter out the positive controls
target_bed <- target_bed %>% filter(str_count(target, "rs") > 0)

#Intersect the bed files to see if any fine-mapped variants were found in the repressed targeted regions
intersect_target_bed <- function(target_bed_row, bed_filt, repression_range=1000){
  #Get position of repression range
  target_chr = as.numeric(str_replace(target_bed_row[1], "chr", ""))
  target_start = as.numeric(target_bed_row[2])
  target_end = as.numeric(target_bed_row[3])
  target_name = target_bed_row[4]
  #Filter for SNPs in target repression range
  temp <- bed_filt[which(bed_filt$chr == target_chr & target_start - repression_range <= bed_filt$start & target_end + repression_range >= bed_filt$end),]
  #Check if any fine-mapped variants were found in the repressed region
  if (nrow(temp) == 0) {
    return(NULL)
  } else {
    temp2 <- suppressWarnings(cbind.data.frame(matrix(unlist(rep(target_bed_row, nrow(temp))), ncol = length(target_bed_row), byrow = TRUE), temp))
    return(temp2)
  }
}
#Apply the function over the targeting sgRNAs and then remove any non-intersections before forming into a dataframe
intersection_list <- pbapply(target_bed, FUN = intersect_target_bed, bed_filt=bed_raw, MARGIN = 1)
intersection_list <- intersection_list[lapply(intersection_list, FUN = is.null) == FALSE]
targeted_signal_df <- bind_rows(intersection_list)
colnames(targeted_signal_df) <- c("target_chr", "target_start", "target_end", "target_name", colnames(bed_raw))
#Create a grna group column in the data frame of targeted signals to group the targets by
targeted_signal_df <- targeted_signal_df %>% mutate(grna_group = str_replace(target_name, "_.+$", ""))
#Also append a signal column
targeted_signal_df <- targeted_signal_df %>% mutate(signal = str_replace(signal.snp_id, "((\\.[^\\.]+){3}\\.).*", "\\1")) %>% mutate(signal = str_replace(signal, "\\.$", ""))
#Identify the number of unique signals that were targeted in the screen and how many targets that lines up with
targeted_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow #Targeted 29 signals
targeted_signal_df %>% dplyr::select(signal) %>% unique %>% nrow #The 29 signals account for 30 Targets

#Make distribution of PIPs plot
fine_map_hist <- targeted_signal_df %>% dplyr::select(rsid, total_pip) %>% distinct() %>% ggplot(aes(x = total_pip)) + geom_histogram(binwidth = 0.1) + 
  xlab("PIP") + ylab("Fine-Mapped Variants") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=45, hjust = 1, vjust = 1), 
        legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title = element_text(size = 24))
tiff(paste0(plot_dir, "CAFEH.fine-mapped_pips.hist.tif"), width = 10800, height=10800, res=1000)
print(fine_map_hist)
dev.off()

#Read in the file of sceptre results to identify how many signals had genes identified if any
sceptre_raw <- read.csv(sceptre_file, sep = "\t", header = TRUE)
#Merge the dataframes by the grna group names
targeted_perturbed_signal_df <- targeted_signal_df %>% inner_join(y = sceptre_raw, by = "grna_group")
#Identify the number of unique signals that were actually perturbed in the screen and how many targets that lines up with
targeted_perturbed_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow #5 targets
targeted_perturbed_signal_df %>% dplyr::select(signal) %>% unique %>% nrow #5 signals
#Also count the number of genes
targeted_perturbed_signal_df %>% dplyr::select(response_id) %>% unique %>% nrow #7 genes could be involved

#Make an output dataframe
targeted_perturbed_signal_df %>% group_by(grna_group) %>% summarise(genes=paste(unique(response_id), collapse = ","),
                                                                    log2_fold_change=paste(unique(round(log_2_fold_change, digits = 3)), collapse = ","),
                                                                    perturb_pvalue_BH=paste(unique(format(p_value_BH, digits = 3, scientific=TRUE)), collapse = ","),
                                                                    fine_mapped_signal=paste(unique(signal), collapse = ","),
                                                                    fine_mapped_snps=paste(unique(rsid), collapse = ","),
                                                                    pip=paste(unique(round(total_pip, digits = 3)), collapse = ","),
                                                                    other_traits=paste(unique(active_traits), collapse = ",")) %>%
  write.table(file = paste0(inp_dir, "supplement.screen_intersected_perturbed.fine-mapped_targets.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
