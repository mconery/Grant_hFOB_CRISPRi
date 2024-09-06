################################################################################

#make_gene_kd_qpcr_figures.R

#The purpose of this script is to make figures and tables for the siRNA KD 
#validations in the osteoblastogenesis and adipogenesis assays.

################################################################################


# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(tidyverse)
library(ggpubr)

#Set directories and file locations
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/osteoblast_assays/"
hmsc_osteo_file <- paste0(inp_dir, "hMSC-BMP2_osteo_kd_qpcr.tsv")
hmsc_adipo_file <- paste0(inp_dir, "hMSC-BMP2_adipo_kd_qpcr.tsv")
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/osteoblast_assays/"

#Set random seed
set.seed(5)

# 1) Read in files ====

#Read in files
hmsc_osteo_raw <- read.csv(hmsc_osteo_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_adipo_raw <- read.csv(hmsc_adipo_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# 2)  Make functions for generating comparison df ====

#Define function
calc_comparison_df <- function(inp_raw){
  #Get unique list of siRNAs/grnas names
  sirnas <- inp_raw %>% dplyr::filter(siRNA != "CON") %>% dplyr::select(siRNA) %>% unique()
  sirnas <- sirnas[,1]
  #Calculate p-values for all the sirnas and bind into a dataframe
  sirna_p_list <- lapply(sirnas, calc_comparison_sirna, inp_raw=inp_raw)
  sirna_df <- dplyr::bind_rows(sirna_p_list)
  return(sirna_df)
}
#Calc comparisons for a given siRNA
calc_comparison_sirna <- function(sirna, inp_raw){
  #Filter for the control and given sirna only
  temp <- inp_raw %>% dplyr::filter(siRNA %in% c(sirna, "CON", "Control") & Gene == sirna)
  #Execute Tests
  treats <- temp %>% select(Treatment) %>% unique
  treats <- treats[,1]
  #Execute test of comparing siRNA knockdown to control in each treatment state
  treat_p_vals <- t(vapply(treats, FUN = calc_sirna_test_treat, FUN.VALUE = numeric(2), temp=temp, sirna=sirna))
  #Form into dataframe
  sirna_df <- cbind.data.frame(Gene = sirna, Treatment=treats, reps=treat_p_vals[,1], p=treat_p_vals[,2])
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

# 3) Call functions and Correct for Multiple Testing ====

#Call functions
hmsc_osteo_comparison_df <- calc_comparison_df(hmsc_osteo_raw)
hmsc_adipo_comparison_df <- calc_comparison_df(hmsc_adipo_raw)

#Remove results with fewer than three replicates
remove_few_rep_results <- function(comparison_df, rep_cutoff=3){
  temp <- comparison_df %>% dplyr::mutate(p = ifelse(reps < rep_cutoff, NA, p))
  return(temp)
}
hmsc_osteo_comparison_df <- remove_few_rep_results(hmsc_osteo_comparison_df)
hmsc_adipo_comparison_df <- remove_few_rep_results(hmsc_adipo_comparison_df)

#Adjust p-values
adjust_comparison_p <- function(hmsc_comparison_df){
  hmsc_comparison_df <- hmsc_comparison_df %>% dplyr::mutate(p.adj = p.adjust(p)) %>% dplyr::mutate(signif = ifelse(p.adj < 0.05, TRUE, FALSE))
  return(hmsc_comparison_df)
}
hmsc_osteo_comparison_df <- adjust_comparison_p(hmsc_osteo_comparison_df)
hmsc_adipo_comparison_df <- adjust_comparison_p(hmsc_adipo_comparison_df)

#Prep for merging by combining columns as identifiers
hmsc_osteo_comparison_df <- hmsc_osteo_comparison_df %>% dplyr::mutate(Gene_Treatment=paste0(Gene, "_", Treatment)) %>% dplyr::select(Gene_Treatment, p.adj, signif)
hmsc_adipo_comparison_df <- hmsc_adipo_comparison_df %>% dplyr::mutate(Gene_Treatment=paste0(Gene, "_", Treatment)) %>% dplyr::select(Gene_Treatment, p.adj, signif)
hmsc_osteo_raw <- hmsc_osteo_raw %>% dplyr::mutate(Gene_Treatment=paste0(Gene, "_", Treatment))
hmsc_adipo_raw <- hmsc_adipo_raw %>% dplyr::mutate(Gene_Treatment=paste0(Gene, "_", Treatment))
#Merge significance info onto the input dataframes
hmsc_osteo_merge <- inner_join(hmsc_osteo_raw, hmsc_osteo_comparison_df, by = "Gene_Treatment")
hmsc_adipo_merge <- inner_join(hmsc_adipo_raw, hmsc_adipo_comparison_df, by = "Gene_Treatment")

#Calculate expression differences
calculate_expr_diff <- function(hmsc_merge){
  #Split out control sirnas
  con_df <- hmsc_merge %>% dplyr::filter(siRNA == "CON") %>% dplyr::arrange(Donor, Gene, Treatment)
  target_df <- hmsc_merge %>% dplyr::filter(siRNA != "CON") %>% dplyr::arrange(Donor, Gene, Treatment)
  #Creater diff_df
  diff_df <- target_df %>% mutate(pct_diff = (target_df$Value - con_df$Value)/(con_df$Value))
  return(diff_df)
}
hmsc_osteo_diff <- calculate_expr_diff(hmsc_osteo_merge)
hmsc_adipo_diff <- calculate_expr_diff(hmsc_adipo_merge)

# 4) Make Box Plots of Results ====


#Make a function that makes a plot for each dataframe 
make_combined_plot <- function(hmsc_normal, out_dir, treat_type, max_val=3){
  #squish max values for plotting
  hmsc_normal <- hmsc_normal %>% mutate(pct_diff = ifelse(abs(pct_diff) > max_val, sign(pct_diff) * max_val, pct_diff))
  #Get y min and max and plot breaks
  y_lims <- c(min(min(hmsc_normal$pct_diff, na.rm = TRUE),-1), max(max(hmsc_normal$pct_diff, na.rm = TRUE),1))
  y_breaks <- seq(((y_lims[1] %/% 0.5) + 1 * (y_lims[1] %% 0.5)) * 0.5, y_lims[2], 1)
  y_break_labels <- y_breaks * 100 %>% round(digits = 0)
  y_break_labels <- paste0(y_break_labels, "%")
  #Get order of names for plots
  gene_order <- hmsc_normal$Gene %>% unique
  #Factorize variables
  hmsc_normal$subplot = factor(hmsc_normal$Treatment, levels = hmsc_normal$Treatment %>% unique)
  hmsc_normal$signif = factor(ifelse(is.na(hmsc_normal$signif), "N", ifelse(hmsc_normal$signif, "T", "F")), levels = c("F", "T", "N"))
  hmsc_normal$Gene = factor(hmsc_normal$Gene, levels = gene_order)
  #Make combined plot in tiff format
  tiff(paste0(out_dir, paste0("kd_qpcr.", treat_type, ".tif")), width = 12000, height=12000, res=1000)
  print(ggplot(hmsc_normal, aes(x = Gene, y = pct_diff, color=signif)) +
    geom_boxplot(linewidth=1.03, width=0.5) + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
    facet_grid(subplot ~ ., scales = "free", space = "free_x") +
    scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05     ", "Not Tested     ")) + 
    ylab("Percentage Decrease in Expression with siRNA Transduction") + 
    scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_break_labels) + 
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
                   axis.title.x = element_blank(), legend.position = "bottom",
                   axis.title.y = element_text(size=24),
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
  #Save lower res version too
  jpeg(paste0(out_dir, paste0("kd_qpcr.", treat_type, ".jpeg")), width = 3600, height=3600, res=300)
  print(ggplot(hmsc_normal, aes(x = Gene, y = pct_diff, color=signif)) +
          geom_boxplot(linewidth=1.03, width=0.5) + 
          geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
          facet_grid(subplot ~ ., scales = "free", space = "free_x") +
          scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05     ", "Not Tested     ")) + 
          ylab("Percentage Change in Expression with siRNA Transduction") + 
          scale_y_continuous(limits = y_lims, breaks = y_breaks, labels = y_break_labels) + 
          ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
                         axis.title.x = element_blank(), legend.position = "bottom",
                         axis.title.y = element_text(size=24),
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

#Call function
make_combined_plot(hmsc_osteo_diff, out_dir, "Osteo")
make_combined_plot(hmsc_adipo_diff, out_dir, "Adipo")

# 5) Make Supplemental Tables of Results ====

### hMSC osteo ###
hmsc_osteo_supplement <- hmsc_osteo_diff %>% dplyr::group_by(siRNA, Gene, Treatment) %>% 
  dplyr::summarize(median_pct_diff=median(pct_diff, na.rm = TRUE), reps = n(), p.adj = first(p.adj)) %>% 
  dplyr::select(Gene, Treatment, reps, median_pct_diff, p.adj) %>% 
  dplyr::mutate(p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj))
write.table(hmsc_osteo_supplement, file = paste0(inp_dir, "hmsc_osteo_kd_qpcr.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### hMSC adipo ###
hmsc_adipo_supplement <- hmsc_adipo_diff %>% 
  dplyr::mutate(ifelse(is.na(p.adj) == FALSE, p.adj, 10)) %>%
  dplyr::group_by(siRNA, Treatment) %>% 
  dplyr::summarize(median_pct_diff=median(pct_diff, na.rm = TRUE), reps = n(), p.adj = first(p.adj), Gene=first(Gene)) %>% 
  dplyr::mutate(p.adj = as.character(p.adj)) %>% 
  as.data.frame() %>%
  dplyr::mutate(p.adj = ifelse(reps < 3, "Not Tested", p.adj)) %>% 
  dplyr::select(Gene, Treatment, reps, median_pct_diff, p.adj) 
write.table(hmsc_adipo_supplement, file = paste0(inp_dir, "hmsc_adipo_kd_qpcr.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
