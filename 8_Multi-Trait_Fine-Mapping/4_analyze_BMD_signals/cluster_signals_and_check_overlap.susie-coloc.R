################################################################################

#cluster_signals_and_check_overlap.susie-coloc.R

#The purpose of this script is twofold: 
# 1) Cluster the signals identified from SuSie-Coloc fine-mapping by PP4s
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
library(VennDiagram)

#Set directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/SuSiE-Coloc_Results/"
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/SuSiE-Coloc_Results/"
plot_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/multi-trait_fine-mapping/SuSiE-Coloc/"

#Set target bed file
target_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/target_selection/all_targeting_guides.bed"
#Set sceptre results file
sceptre_file <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/sceptre_results/discovery_results.tsv"

#Set filtering criteria
purity_thresh=0.1
bmd_sig_thresh=5e-8
pp4_thresh=0.95

# 1) Read in Signal, PP4, and Bed Files and Filter them ====

#Set name of signal file and pp4 file
signal_raw <- read.csv(paste0(inp_dir, "master.susie-coloc.tsv"), header = TRUE, sep = "\t")
pp4_raw <- read.csv(paste0(inp_dir, "master.signed_pp4s.tsv"), header = TRUE, sep = "\t", row.names = 1)
bed_raw <- read.csv(paste0(inp_dir, "master.variants.bed"), header = TRUE, sep = "\t")

#Filter the signal table for the pure high significance signals
signal_filt <- signal_raw %>%
  mutate(PP4 = ifelse(is.na(other_traits) | other_traits == "", "", PP4))
#Write the filtered signal table to file  
write.table(signal_filt, paste0(inp_dir, "master.filtered.susie-coloc.tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#Filter the pp4 table down for the same signals
pp4_raw[rownames(pp4_raw) %in% signal_filt$signal_id,] %>% 
  write.table(paste0(inp_dir, "master.filtered.signed_pp4s.tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

#Count number of traits with a colocalization
pp4_binarized = pp4_raw >= pp4_thresh | pp4_raw <= -pp4_thresh
((pp4_binarized %>% colSums) > 0) %>% sum 

#Get name of most pleiotropic signal
pp4_binarized %>% rowSums %>% which.max %>% names
pp4_binarized %>% rowSums %>% max

#Filter down the variant table as well
get_signal <- function(signal_variant){sub("\\.[^.]*$", "", signal_variant)}
bed_filt <- bed_raw %>% filter(vapply(bed_raw$signal.snp_id, get_signal, FUN.VALUE = character(1)) %in%  signal_filt$signal_id)
check_thresh_traits <- function(PP4s, pp4_cut=0.95){PP4s > pp4_cut}
get_good_trait_values <- function(i, trait_values, PP4s=bed_filt$colocalizing_pp4s){paste0(str_split(trait_values, pattern = ",")[[i]][check_thresh_traits(unlist(str_split(PP4s[i], ",")))], collapse = ",")}
bed_filt <- bed_filt %>% mutate(colocalizing_traits = vapply(seq(1, nrow(bed_filt), 1), FUN = get_good_trait_values, trait_values=colocalizing_traits, FUN.VALUE = character(1)),
                                colocalizing_traits = vapply(seq(1, nrow(bed_filt), 1), FUN = get_good_trait_values, trait_values=colocalizing_pp4s, FUN.VALUE = character(1)))
#Write it to file
write.table(paste0(inp_dir, "master.filtered.variants.bed"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 2) Make Clustering and Signal Count Plots ====

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
count_bmd_signals <- function(inp_file, out_plot_dir=plot_dir, pp4_thresh=0.95){
  temp <- read.csv(inp_file, header = TRUE, sep = "\t")
  #Recount traits based on pp4_thresh
  count_thresh_traits <- function(PP4s, pp4_cut=pp4_thresh){sum(PP4s > pp4_cut)}
  check_thresh_traits <- function(PP4s, pp4_cut=pp4_thresh){PP4s > pp4_cut}
  get_good_trait_values <- function(i, trait_values, PP4s=temp$PP4){paste0(str_split(trait_values, pattern = ",")[[i]][check_thresh_traits(unlist(str_split(PP4s[i], ",")))], collapse = ",")}
  temp <- temp %>% mutate(num_other_traits_cut = lapply(str_split(PP4, ","), count_thresh_traits),
                          other_traits_cut = vapply(seq(1, nrow(temp), 1), FUN = get_good_trait_values, trait_values=other_traits, FUN.VALUE = character(1)))
  print(nrow(temp))
  #Calc number of BMD only-traits and traits with highest counts
  temp %>% dplyr::select(num_other_traits_cut) %>% dplyr::filter(num_other_traits_cut > 0) %>% nrow %>% print
  temp %>% dplyr::select(num_other_traits_cut) %>% dplyr::filter(num_other_traits_cut == 0) %>% nrow %>% print
  temp2 <- temp %>% dplyr::select(other_traits_cut) %>% dplyr::filter(other_traits_cut != "") %>% as.vector() %>% unlist %>% str_split(pattern = ",") %>%
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
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.signals_per_trait.tiff"), width = 10800, height=6000, res=1000)
    print(t_plot)
  dev.off()
  #Output the number of traits per signal
  temp3 <- cbind.data.frame(y = temp$num_other_traits)
  #Make plot
  t_plot <- temp3 %>% ggplot(aes(x = y)) + geom_histogram(binwidth = 1) + 
    xlab("Number of Non-BMD Traits") + ylab("Count of Signals") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 24, angle=45, hjust = 1, vjust = 1), 
          legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title = element_text(size = 24))
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.traits_per_signal.trait.tiff"), width = 10800, height=10800, res=1000)
    print(t_plot)
  dev.off()
}
#Call function
count_bmd_signals(paste0(inp_dir, "master.filtered.susie-coloc.tsv"))

#Make a function to cluster by PP4s
cluster_by_activity <- function(activity_file, out_plot_dir=plot_dir, input_dir=inp_dir, activity_thresh=0.95, binarize=FALSE){
  #Read in files and cluster
  inp_raw <- read.table(activity_file, header = TRUE)
  #Remove completely empty columns
  inp_filt <- inp_raw[,colSums(inp_raw) != 0] 
  #Take absolute values
  inp_filt <- abs(inp_filt)
  #Binarize if appropriate 
  if (binarize) {
    inp_filt <- ifelse(inp_filt > activity_thresh, 1, 0) %>% as.data.frame()
    name_addendum <- paste0("binarize_", activity_thresh, ".")
  } else {
    name_addendum <- ""
  }
  #Rename column names
  colnames(inp_filt) <- colnames(inp_filt) %>% correct_trait_names
  #Create needed directory
  dir.create(paste0(out_plot_dir), recursive = TRUE, mode = "0777", showWarnings = FALSE)
  #Make heatmaps
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "pp4s.heatmap.tiff"), width = 7200, height = 10800, res = 1000)
  print(pheatmap(t(inp_filt), 
           show_colnames = FALSE, color=colorRampPalette(c("#C5C6D0", "black"),)(50), fontsize_row = 16,
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean"))
  dev.off()
  
  # Set seed
  set.seed(5)
  # Perform UMAP transformation
  umap_result <- umap(inp_filt, n_neighbors = 15, n_components = 2)
  #Plot uncolored UMAP clusters
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "pp4s.umap.tiff"), width = 6, height = 11.5, units = 'in', res = 500)
  plot(umap_result$layout, pch = 16, main = "", xlab = "UMAP 1", ylab = "UMAP 2")
  dev.off()
  # Plot the UMAP results with KNN clusters colored by trait
  plot_umap_colored_trait <- function(trait, umap_result, inp_filt, plot_dir, file_prefix){
    file_trait <- trait %>% str_replace_all(pattern = " ", replacement = "_") %>% str_replace_all(pattern = "/", replacement = "_")
    trait_color = inp_filt %>% mutate(trait_col = ifelse(inp_filt[,trait] > activity_thresh, 1, 0)) %>% select(trait_col)
    tiff(paste0(plot_dir, "SuSiE-Coloc.", name_addendum, "pp4s.umap.", file_trait, ".tiff"), width = 6, height = 11.5, units = 'in', res = 1000)
    plot(umap_result$layout, col = as.numeric(trait_color$trait_col) + 1, pch = 16, main = paste0(trait, " UMAP"), xlab = "UMAP 1", ylab = "UMAP 2")
    dev.off()
  }
  lapply(colnames(inp_filt), plot_umap_colored_trait, umap_result=umap_result, inp_filt=inp_filt, plot_dir=out_plot_dir, file_prefix="SuS")
  
  #Output a supplemental table
  write.table(inp_filt, file=paste0(inp_dir, "supplement.", name_addendum, "pp4s.SuSie-Coloc.tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
}
#Call function
cluster_by_activity(paste0(inp_dir, "master.filtered.signed_pp4s.tsv"))
cluster_by_activity(paste0(inp_dir, "master.filtered.signed_pp4s.tsv"), binarize=TRUE)

#Make a function to cluster signals by signed pp4s
cluster_by_weight <- function(activity_file, out_plot_dir=plot_dir, input_dir=inp_dir, activity_thresh=0.95, binarize=FALSE){
  #Read in files and cluster
  inp_raw <- read.table(activity_file, header = TRUE)
  #Remove completely empty columns
  inp_filt <- inp_raw[,colSums(inp_raw) != 0] 
  #Binarize if appropriate 
  if (binarize) {
    inp_filt <- ifelse(inp_filt > activity_thresh, 1, ifelse(inp_filt < -activity_thresh, -1, 0)) %>% as.data.frame()
    name_addendum <- paste0("binarize_", activity_thresh, ".")
  } else {
    name_addendum <- ""
  }
  #Rename column names
  colnames(inp_filt) <- colnames(inp_filt) %>% correct_trait_names
  #Create needed directory
  dir.create(paste0(out_plot_dir), recursive = TRUE, mode = "0777", showWarnings = FALSE)
  #Make heatmaps
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "signed_pp4s.heatmap.tiff"), width = 7200, height = 10800, res = 1000)
  print(pheatmap(t(inp_filt), 
                 show_colnames = FALSE, color=colorRampPalette(c("navy", "#C5C6D0", "red"),)(50), fontsize_row = 16,
                 clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean"))
  dev.off()
  
  # Set seed
  set.seed(5)
  # Perform UMAP transformation
  umap_result <- umap(inp_filt, n_neighbors = 15, n_components = 2)
  #Plot uncolored UMAP clusters
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "Signed_pp4s.umap.tiff"), width = 6, height = 11.5, units = 'in', res = 500)
  plot(umap_result$layout, pch = 16, main = "", xlab = "UMAP 1", ylab = "UMAP 2")
  dev.off()
  # Plot the UMAP results with KNN clusters colored by trait
  plot_umap_colored_trait <- function(trait, umap_result, inp_filt, plot_dir, file_prefix){
    file_trait <- trait %>% str_replace_all(pattern = " ", replacement = "_") %>% str_replace_all(pattern = "/", replacement = "_")
    trait_color = inp_filt %>% mutate(trait_col = ifelse(inp_filt[,trait] > activity_thresh, 1, 0)) %>% select(trait_col)
    tiff(paste0(plot_dir, "SuSiE-Coloc.", name_addendum, "signed_pp4s.umap.", file_trait, ".tiff"), width = 6, height = 11.5, units = 'in', res = 1000)
    plot(umap_result$layout, col = as.numeric(trait_color$trait_col) + 1, pch = 16, main = paste0(trait, " UMAP"), xlab = "UMAP 1", ylab = "UMAP 2")
    dev.off()
  }
  lapply(colnames(inp_filt), plot_umap_colored_trait, umap_result=umap_result, inp_filt=inp_filt, plot_dir=out_plot_dir, file_prefix="SuS")
  
  #Make a barplot of the signals by trait split by positive and negative effects
  t_plot <- inp_filt %>% mutate(signal=rownames(inp_filt)) %>% pivot_longer(!signal, names_to = "trait", values_to = "activity") %>%
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
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "signals_per_trait.weighted.tiff"), width = 10800, height=6000, res=1000)
  print(t_plot)
  dev.off()
  
  #Make the barplot again, but this time make it vertical
  t_plot <- inp_filt %>% mutate(signal=rownames(inp_filt)) %>% pivot_longer(!signal, names_to = "trait", values_to = "activity") %>%
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
  tiff(paste0(out_plot_dir, "SuSiE-Coloc.", name_addendum, "signals_per_trait.weighted.vertical.tiff"), height = 10800, width=6000, res=1000)
  print(t_plot)
  dev.off()
  
  #Write table to file
  write.table(inp_filt, file=paste0(inp_dir, "supplement.", name_addendum, "signed_pp4s.SuSie-Coloc.tsv"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
}
cluster_by_weight(paste0(inp_dir, "master.filtered.signed_pp4s.tsv"))
cluster_by_weight(paste0(inp_dir, "master.filtered.signed_pp4s.tsv"), binarize=TRUE)

# 3) Read in Signal Files and Process for Supplement ====

#Clean up trait names
bed_filt$active_traits <- bed_filt$colocalizing_traits %>% str_replace_all(pattern = ",", replacement = ", ") %>% correct_trait_names
signal_filt$traits <- signal_filt$other_traits %>% str_replace_all(pattern = ",", replacement = ", ") %>% correct_trait_names

#Write files
write.table(signal_raw, file = paste0(inp_dir, "supplement.signals.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(bed_raw, file = paste0(inp_dir, "supplement.variants.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 4) Check for signal intersections with the file of targets ====

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
targeted_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow 
targeted_signal_df %>% dplyr::select(signal) %>% unique %>% nrow 

#Read in the file of sceptre results to identify how many signals had genes identified if any
sceptre_raw <- read.csv(sceptre_file, sep = "\t", header = TRUE)
#Merge the dataframes by the grna group names
targeted_perturbed_signal_df <- targeted_signal_df %>% inner_join(y = sceptre_raw, by = "grna_group")
#Identify the number of unique signals that were actually perturbed in the screen and how many targets that lines up with
targeted_perturbed_signal_df %>% dplyr::select(grna_group) %>% unique %>% nrow #8 targets
targeted_perturbed_signal_df %>% dplyr::select(signal) %>% unique %>% nrow #8 signals
#Also count the number of genes
targeted_perturbed_signal_df %>% dplyr::select(response_id) %>% unique %>% nrow #10 genes could be involved

#Make an output dataframe
targeted_perturbed_signal_df %>% group_by(grna_group) %>% summarise(genes=paste(unique(response_id), collapse = ","),
                                                                    log2_fold_change=paste(unique(round(log_2_fold_change, digits = 3)), collapse = ","),
                                                                    perturb_pvalue_BH=paste(unique(format(p_value_BH, digits = 3, scientific=TRUE)), collapse = ","),
                                                                    fine_mapped_signal=paste(unique(signal), collapse = ","),
                                                                    fine_mapped_snps=paste(unique(rsid), collapse = ","),
                                                                    pip=paste(unique(round(total_pip, digits = 3)), collapse = ","),
                                                                    colocalizing_traits=paste(unique(colocalizing_traits), collapse = ",")) %>%
  write.table(file = paste0(inp_dir, "supplement.screen_intersected_perturbed.susie-fine-mapped_targets.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  

# 5) Compare Overlap of CAFEH Signals ====

#Read in CAFEH Files
cafeh_bed <- read.csv(paste0("C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/summary_files/supplement.variants.purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered.txt"), header = TRUE, sep = "\t")
cafeh_signals <- read.csv("C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/multi-trait_fine-mapping/summary_files/supplement.signals.purity_0.5.activity_0.95.gwas_5e-08.highest_residual_filtered.txt", header = TRUE, sep = "\t")
#Remove spurious signals from cafeh_bed file
cafeh_bed <- cafeh_bed %>% mutate(cafeh_signal=sub("\\.[^.]*$", "", signal.snp_id)) %>% dplyr::filter(cafeh_signal %in% cafeh_signals$signal_id)
#Find intersection of variant sets
cafeh_intersect <- cafeh_bed %>% inner_join(bed_raw, by = c("chr", "start", "end"))
#Extract signal IDs
cafeh_intersect <- cafeh_intersect %>% mutate(cafeh_signal=sub("\\.[^.]*$", "", signal.snp_id.x), 
                                              susie_signal=sub("\\.[^.]*$", "", signal.snp_id.y)) %>%
  dplyr::select(cafeh_signal, susie_signal) %>% distinct()
# 371 unique pairings
#Count intersection
cafeh_intersect$cafeh_signal %>% unique %>% length
cafeh_intersect$susie_signal %>% unique %>% length
#Make Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = nrow(cafeh_signals),          # size of set A
  area2 = nrow(signal_filt),           # size of set B
  cross.area = cafeh_intersect$susie_signal %>% unique %>% length,      # size of intersection A ∩ B
  fill = c("red", "blue"),
  filename = paste0(plot_dir, "CAFEH-SuSiE_Venn.png"),  # Output file path
  height = 2000,                   # Height in pixels
  width = 2000,                    # Width in pixels
  resolution = 300
)

#Compare credible set sizes of the intersecting signals
overlap_susie_sizes <- signal_filt %>% filter(signal_id %in% cafeh_intersect$susie_signal) %>% dplyr::select(signal_id, set_size)
overlap_cafeh_sizes <- cafeh_signals %>% filter(signal_id %in% cafeh_intersect$cafeh_signal) %>% dplyr::select(signal_id, num_snps)
#Link sizes together
size_compare <- cafeh_intersect %>% inner_join(overlap_cafeh_sizes, by = c("cafeh_signal" = "signal_id")) %>%
  inner_join(overlap_susie_sizes, by = c("susie_signal" = "signal_id"))
size_test <- wilcox.test(size_compare$num_snps, size_compare$set_size, paired = TRUE, alternative = "two.sided")
mean(overlap_susie_sizes$set_size)
mean(overlap_cafeh_sizes$num_snps)