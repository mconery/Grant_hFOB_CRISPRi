################################################################################

#make_assay_marker_figures.R

#The purpose of this script is to make figures and tables for the osteoblast
#and adipocyte marker gene qPCR experiments that were run in parallel with the 
#siRNA knockdown assays.

################################################################################


# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(tidyverse)
library(ggpubr)

#Set directories and file locations
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/osteoblast_assays/"
hmsc_osteo_file <- paste0(inp_dir, "hMSC-BMP2_osteo_markers.tsv")
hmsc_adipo_file <- paste0(inp_dir, "hMSC-BMP2_adipo_markers.tsv")
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/osteoblast_assays/"

#Set random seed
set.seed(5)

# 1) Read in files ====

#Read in files
hmsc_osteo_raw <- read.csv(hmsc_osteo_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_adipo_raw <- read.csv(hmsc_adipo_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# 2) Normalize Data to Control Level of Each Plate ====

#Make a function that normalizes each plate relative to its control for the hMSC assays and call it
normalize_hmsc <- function(inp_raw, treat_name){
  #Get names of genes
  genes <- inp_raw$Gene %>% unique
  return(bind_rows(lapply(genes, FUN = normalize_hmsc_gene, inp_raw=inp_raw, treat_name=treat_name)))
}
normalize_hmsc_gene <- function(gene, inp_raw, treat_name){
  #Filter for the gene
  inp_raw <- inp_raw %>% dplyr::filter(Gene == gene)
  #Bind on a column representing unique plate
  inp_raw <- inp_raw %>% mutate(Unique_plate = paste(Donor, Plate, Gene, sep = "."))
  unique_plates = unique(inp_raw$Unique_plate)
  bind_rows(lapply(unique_plates, normalize_hmsc_plate, inp_raw=inp_raw, treat_name=treat_name))
}
normalize_hmsc_plate <- function(unique_plate, inp_raw, treat_name){
  temp <- inp_raw %>% dplyr::filter(Unique_plate == unique_plate)
  temp$Value <- temp$Value/temp[temp$siRNA == "CON" & temp$Treatment == treat_name,"Value"]
  return(temp)
}
hmsc_osteo_normal <- normalize_hmsc(hmsc_osteo_raw, treat_name = "BMP2")
hmsc_adipo_normal <- normalize_hmsc(hmsc_adipo_raw, treat_name = "Adipo")

# 3) Make functions for generating comparison df ====

#Define function
calc_comparison_df <- function(inp_raw){
  #Get unique list of siRNAs/grnas names
  sirnas_genes <- inp_raw %>% dplyr::filter(siRNA != "CON") %>% dplyr::mutate(siRNA_Gene = paste0(siRNA, "_", Gene)) %>% dplyr::select(siRNA_Gene) %>% unique()
  sirnas_genes <- sirnas_genes[,1]
  #Calculate p-values for all the sirnas and bind into a dataframe
  sirna_p_list <- lapply(sirnas_genes, calc_comparison_sirna_gene, inp_raw=inp_raw)
  sirna_df <- dplyr::bind_rows(sirna_p_list)
  #Run test for controls
  treats = unique(inp_raw$Treatment)
  #Need to test multiple plates worth of controls
  genes <- inp_raw %>% dplyr::select(Gene) %>% unique
  genes <- genes[,1]
  temp <- t(vapply(genes, calc_treat_test_control, FUN.VALUE = numeric(2), inp_raw=inp_raw))
  temp <- cbind.data.frame(genes, t(paste0("CON", "-", treats)), temp)
  colnames(temp) <- c("Gene", "group1", "group2", "reps", "p")
  sirna_df <- rbind.data.frame(sirna_df, temp)
  return(sirna_df)
}
#Calc comparisons for a given siRNA
calc_comparison_sirna_gene <- function(sirna_gene, inp_raw){
  #Split the sirna and gene
  sirna_gene <- str_split(sirna_gene, "_") %>% unlist
  sirna <- sirna_gene[1]
  gene <- sirna_gene[2]
  #Filter for the control and given sirna only
  temp <- inp_raw %>% dplyr::filter(siRNA %in% c(sirna, "CON", "Control") & Gene == gene)
  #Check plate of sirna (This step assumes that all reps of a given sirna were run on the same plate)
  if ("Plate" %in% colnames(temp)) {
    check_plate = temp %>% filter(siRNA == sirna) %>% select(Plate) %>% unique %>% as.numeric
    temp <- temp %>% filter(Plate == check_plate)
  }
  #Check if Measurements are present
  if ("Measurement" %in% colnames(temp)) {
    temp <- temp %>% dplyr::group_by(siRNA, Replicate) %>% dplyr::summarize(Value=mean(Value))
  }
  #Check if Multiple passages are present for donors (Passage is mutually exclusive with Measurement in input file)
  if ("Passage" %in% colnames(temp)){
    temp <- temp %>% dplyr::group_by(siRNA, Donor, Treatment) %>% dplyr::summarize(Value=mean(Value)) %>% as.data.frame
  }
  #Execute Tests
  treats <- temp %>% select(Treatment) %>% unique
  treats <- treats[,1]
  #Execute test of comparing siRNA knockdown to control in each treatment state
  treat_p_vals <- t(vapply(treats, FUN = calc_sirna_test_treat, FUN.VALUE = numeric(2), temp=temp, sirna=sirna))
  #Execute test of comparing treated to untreated for the sirna and append to dataframe
  treat_untreat_p_val <- calc_treat_test_sirna(sirna, temp)
  #Create output df for the function after checking for whether plates need to be appended to the controls
  if ("Plate" %in% colnames(temp)) {
    plate = temp %>% dplyr::select("Plate") %>% unique %>% as.character
    sirna_df <- cbind.data.frame(Gene = gene, group1 = paste(sirna, treats, sep = "-"), group2 = paste0("CON", "-", treats), reps=treat_p_vals[,1], p=treat_p_vals[,2])
  } else {
    sirna_df <- cbind.data.frame(Gene = gene, group1 = paste(sirna, treats, sep = "-"), group2 = paste("CON", treats, sep = "-"), reps=treat_p_vals[,1], p=treat_p_vals[,2])
  }
  sirna_df <- rbind.data.frame(sirna_df, c(gene, paste(sirna, treats, sep = "-"), treat_untreat_p_val))
  return(sirna_df)
}
#Make functions for executing tests at sirna level when treatment present
calc_sirna_test_treat <- function(treat, temp, sirna){
  temp_two <- temp %>% filter(Treatment == treat)
  #Filter for donor controls with matched sirna treatment
  sirna_donors <- temp_two %>% dplyr::filter(siRNA == sirna & is.na(Value) == FALSE) %>% dplyr::select(Donor)
  sirna_donors <- sirna_donors[,1]
  sirna_donors <- temp_two %>% dplyr::filter(siRNA == "CON" & Donor %in% sirna_donors & is.na(Value) == FALSE) %>% dplyr::select(Donor)
  sirna_donors <- sirna_donors[,1]
  temp_two <- temp_two %>% dplyr::filter(Donor %in% sirna_donors) %>% dplyr::arrange(Donor)
  control_vals <- temp_two %>% filter(siRNA == "CON") %>% select(Value)
  control_vals <- control_vals[,1]
  sirna_vals <- temp_two %>% filter(siRNA == sirna) %>% select(Value)
  sirna_vals <- sirna_vals[,1]
  if (length(sirna_vals) > 0) {
    test_result <- t.test(sirna_vals, control_vals, alternative = "less", paired = TRUE)
  } else {
    test_result = 0
  }
  return(c(length(sirna_vals), test_result$p.value))
}
calc_treat_test_sirna <- function(sirna, temp){
  temp_two <- temp %>% filter(siRNA == sirna)
  #Filter for donor controls with matched sirna treatment
  sirna_donors <- temp_two %>% dplyr::filter(Treatment != "CON" & is.na(Value) == FALSE) %>% dplyr::select(Donor) %>% unique
  sirna_donors <- sirna_donors[,1]
  sirna_donors <- temp_two %>% dplyr::filter(Treatment == "CON" & Donor %in% sirna_donors & is.na(Value) == FALSE) %>% dplyr::select(Donor) %>% unique
  sirna_donors <- sirna_donors[,1]
  temp_two <- temp_two %>% dplyr::filter(Donor %in% sirna_donors) %>% dplyr::arrange(Donor)
  control_vals <- temp_two %>% dplyr::filter(Treatment == "CON") %>% dplyr::select(Value)
  control_vals <- control_vals[,1]
  treat_vals <- temp_two %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(Value)
  treat_vals <- treat_vals[,1]
  #Introduce error check for constant vectors
  if (all(control_vals == median(control_vals)) && all(treat_vals == median(treat_vals))) {
    return(c(length(treat_vals), NA))
  } else {
    test_result <- t.test(control_vals, treat_vals, alternative = "less", paired = TRUE)
    return(c(length(treat_vals), test_result$p.value))
  }
}
#Make function that tests treatment relative to control for control siRNA
calc_treat_test_control <- function(gene, inp_raw){
  temp <- inp_raw %>% dplyr::filter(Gene == gene & siRNA == "CON")
  #Average plates together
  temp <- temp %>% dplyr::group_by(siRNA, Treatment, Donor) %>% dplyr::summarize(Value=mean(Value, na.rm = TRUE)) %>% as.data.frame
  return(calc_treat_test_sirna("CON", temp))
}

# 4) Call functions and Correct for Multiple Testing ====

#Call functions
hmsc_osteo_comparison_df <- calc_comparison_df(hmsc_osteo_normal)
hmsc_adipo_comparison_df <- calc_comparison_df(hmsc_adipo_normal)

#Remove results with fewer than three replicates
remove_few_rep_results <- function(comparison_df, rep_cutoff=3){
  temp <- comparison_df %>% dplyr::mutate(p = ifelse(reps < rep_cutoff, NA, p))
  return(temp)
}
hmsc_osteo_comparison_df <- remove_few_rep_results(hmsc_osteo_comparison_df)
hmsc_adipo_comparison_df <- remove_few_rep_results(hmsc_adipo_comparison_df)

#Correct for multiple testing separately for the three test types in the hMSC osteo experiment (and split by gene)
hmsc_osteo_comparison_df <- hmsc_osteo_comparison_df %>% mutate(p.adj = 1)
#Get Genes
genes <- hmsc_osteo_comparison_df$Gene %>% unique
for (gene in genes) {
  hmsc_osteo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_osteo_comparison_df$p[str_detect(pattern = "-CON",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene], method = "BH")
  hmsc_osteo_comparison_df$p.adj[str_detect(pattern = "-BMP2",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_osteo_comparison_df$p[str_detect(pattern = "-BMP2",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene], method = "BH")
  hmsc_osteo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_osteo_comparison_df$p[str_detect(pattern = "-CON",hmsc_osteo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_osteo_comparison_df$group2) & hmsc_osteo_comparison_df$Gene == gene], method = "BH")
}
#Correct for multiple testing separately for the three test types in the hMSC adipo experiment (and split by gene)
hmsc_adipo_comparison_df <- hmsc_adipo_comparison_df %>% mutate(p.adj = 1)
#Get Genes
genes <- hmsc_adipo_comparison_df$Gene %>% unique
for (gene in genes) {
  hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene], method = "BH")
  hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-Adipo",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-Adipo",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene], method = "BH")
  hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene] <-
    p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_adipo_comparison_df$group2) & hmsc_adipo_comparison_df$Gene == gene], method = "BH")
}

# 5) Make Box Plots of Results ====

#Get order of names for plots
sirna_order <- sort(unique(as.character(hmsc_osteo_normal$siRNA[str_detect(pattern = "CON", hmsc_osteo_normal$Name) == FALSE])))

#Make a function that makes a plot for each dataframe 
make_combined_plot <- function(hmsc_normal, hmsc_comparison_df, out_dir, treat_type){
  #Get gene names
  genes <- hmsc_normal$Gene %>% unique
  #Apply facet function
  facet_dfs <- lapply(genes, make_facet_df, hmsc_normal = hmsc_normal, hmsc_comparison_df=hmsc_comparison_df, out_dir=out_dir, treat_type=treat_type)
  #Combine facet dfs into a single dataframe
  combined_plot_df <- bind_rows(facet_dfs)
  combined_plot_df$subplot <- factor(combined_plot_df$subplot, levels = genes)
  #Make combined plot in tiff format
  tiff(paste0(out_dir, paste0("marker_gene_qpcr.", treat_type, ".tiff")), width = 12000, height=12000, res=1000)
  print(ggplot(combined_plot_df, aes(x = siRNA, y = Value, color=signif)) +
    geom_boxplot(linewidth=1.03, width=0.5) + 
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
    facet_grid(subplot ~ ., scales = "free", space = "free_x") +
    scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05     ", "Not Tested     ")) + 
    ylab("Expression Fold Change") + 
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
                   axis.title.x = element_blank(), legend.position = "bottom",
                   axis.title.y = element_text(size = 24),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
                   legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                   legend.title = ggplot2::element_blank(),
                   legend.text = element_text(size = 20),
                   axis.text.y = element_text(size = 20),
                   strip.text.y = element_text(size = 20)) +
    ggplot2::guides(color = ggplot2::guide_legend(
      keywidth = 0.0,
      keyheight = 0.3, 
      default.unit = "inch",
      override.aes = list(size = 5))))
  dev.off()
}
#Make a sub function for making individual facets
make_facet_df <- function(gene, hmsc_normal, hmsc_comparison_df, out_dir, treat_type, sirna_plot_order=sirna_order){
  #Filter the normal data frame
  normal_temp <- hmsc_normal %>% dplyr::filter(Gene == gene, Treatment == treat_type, siRNA != "CON")
  normal_temp$siRNA <- factor(normal_temp$siRNA, levels = sirna_order)
  #Now filter the comparison df correctly
  comparison_temp <- hmsc_comparison_df %>% dplyr::filter(Gene == gene) %>% 
    dplyr::filter(str_detect(group1, pattern = paste0("-", treat_type)) & str_detect(group2, pattern = paste0("-", treat_type)) & str_detect(group2, pattern = "CON")) %>%
    mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
    mutate(siRNA = str_replace(group1, pattern=paste0("-",treat_type), replacement="")) %>%
    dplyr::select(siRNA, p.signif)
  #Combine the significance info onto the normal_temp
  normal_temp <- normal_temp %>% inner_join(comparison_temp, by = "siRNA")
  normal_temp$signif <- factor(ifelse(is.na(normal_temp$p.signif), "N", ifelse(normal_temp$p.signif == "*", "T", "F")), levels=c("F","T", "N"))
  #Make plot_df
  plot_df <- normal_temp %>% dplyr::select(siRNA, signif, Value) %>% mutate(subplot = gene)
  return(plot_df)
}

#Call function
make_combined_plot(hmsc_osteo_normal, hmsc_osteo_comparison_df, out_dir, "BMP2")
make_combined_plot(hmsc_adipo_normal, hmsc_adipo_comparison_df, out_dir, "Adipo")

# 6) Make Supplemental Tables of Results ====

### hMSC osteo ###
hmsc_osteo_comparison_filt <- hmsc_osteo_comparison_df %>% dplyr::filter(str_detect(group1, "-BMP2") & str_detect(group2, "-BMP2")) %>%
  dplyr::mutate(siRNA_Gene=paste0(str_remove(group1, "-BMP2"), "_", Gene))
hmsc_osteo_supplement <- hmsc_osteo_normal %>% dplyr::group_by(siRNA, Gene) %>% 
  dplyr::summarize(mean_fold_change=mean(Value, na.rm = TRUE)) %>% 
  dplyr::mutate(siRNA_Gene=paste0(siRNA, "_", Gene)) %>% 
  dplyr::select(siRNA, siRNA_Gene, mean_fold_change) %>%
  inner_join(hmsc_osteo_comparison_filt, by = "siRNA_Gene") %>% 
  dplyr::select(Gene, siRNA, reps, mean_fold_change, p, p.adj) %>% 
  dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj))
write.table(hmsc_osteo_supplement, file = paste0(inp_dir, "hmsc_osteo_markers.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### hMSC adipo ###
hmsc_adipo_comparison_filt <- hmsc_adipo_comparison_df %>% dplyr::filter(str_detect(group1, "-Adipo") & str_detect(group2, "-Adipo")) %>%
  dplyr::mutate(siRNA_Gene=paste0(str_remove(group1, "-Adipo"), "_", Gene))
hmsc_adipo_supplement <- hmsc_adipo_normal %>% dplyr::group_by(siRNA, Gene) %>% 
  dplyr::summarize(mean_fold_change=mean(Value, na.rm = TRUE)) %>% 
  dplyr::mutate(siRNA_Gene=paste0(siRNA, "_", Gene)) %>% 
  dplyr::select(siRNA, siRNA_Gene, mean_fold_change) %>%
  inner_join(hmsc_adipo_comparison_filt, by = "siRNA_Gene") %>% 
  dplyr::select(Gene, siRNA, reps, mean_fold_change, p, p.adj) %>% 
  dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p))
hmsc_adipo_supplement$p.adj = ifelse(is.na(hmsc_adipo_supplement$p.adj), "Not Tested", hmsc_adipo_supplement$p.adj)
write.table(hmsc_adipo_supplement, file = paste0(inp_dir, "hmsc_adipo_markers.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


