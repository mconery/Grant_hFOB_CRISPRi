################################################################################

#make_assay_figures.R

#The purpose of this script is to make figures and tables for the osteoblast
#cell model assay results. This script runs locally

################################################################################


# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(tidyverse)
library(ggpubr)

#Set directories and file locations
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/osteoblast_assays/"
hfob_alp_file <- paste0(inp_dir, "hfob_alp_results.tsv")
hmsc_alp_file <- paste0(inp_dir, "hMSC-BMP2_alp_results.tsv")
hmsc_ars_file <- paste0(inp_dir, "hMSC-BMP2_ars_results.tsv")
hmsc_adipo_file <- paste0(inp_dir, "hMSC-BMP2_adipo_results.tsv")
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/osteoblast_assays/"

#Set random seed
set.seed(5)

# 1) Read in files ====

#Read in files
hfob_alp_raw <- read.csv(hfob_alp_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_alp_raw <- read.csv(hmsc_alp_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_ars_raw <- read.csv(hmsc_ars_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_adipo_raw <- read.csv(hmsc_adipo_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#Process the hfob data set to keep formatting consistent
#Convert all "Control"s to "CON"s
hfob_alp_raw$siRNA = ifelse(hfob_alp_raw$siRNA == "Control", "CON", hfob_alp_raw$siRNA)
#Add Name column
hfob_alp_raw <- hfob_alp_raw %>% mutate(Name = siRNA)
#Remove no siRNA test
hfob_alp_raw <- hfob_alp_raw %>% dplyr::filter(siRNA != "NT" & siRNA != "Not Treated")

# 2) Normalize Data to Control Level of Each Plate ====

#Make a function that normalizes each plate relative to its control for the hfobs and call it
normalize_hfob <- function(inp_raw){
  replicates = unique(inp_raw$Replicate)
  bind_rows(lapply(replicates, normalize_hfob_replicate, inp_raw=inp_raw))
}
normalize_hfob_replicate <- function(replicate, inp_raw){
  temp <- inp_raw %>% dplyr::filter(Replicate == replicate)
  temp$Value <- temp$Value/mean(temp[temp$siRNA == "CON","Value"])
  return(temp)
}
hfob_alp_normal <- normalize_hfob(hfob_alp_raw)

#Make a function that normalizes each plate relative to its control for the hMSC assays and call it
normalize_hmsc <- function(inp_raw, treat_name){
  #Bind on a column representing unique plate
  inp_raw <- inp_raw %>% mutate(Unique_plate = paste(Donor, Plate, Passage, sep = "."))
  unique_plates = unique(inp_raw$Unique_plate)
  bind_rows(lapply(unique_plates, normalize_hmsc_plate, inp_raw=inp_raw, treat_name=treat_name))
}
normalize_hmsc_plate <- function(unique_plate, inp_raw, treat_name){
  temp <- inp_raw %>% dplyr::filter(Unique_plate == unique_plate)
  temp$Value <- temp$Value/temp[temp$siRNA == "CON" & temp$Treatment == treat_name,"Value"]
  return(temp)
}
hmsc_alp_normal <- normalize_hmsc(hmsc_alp_raw, treat_name = "BMP2")
hmsc_ars_normal <- normalize_hmsc(hmsc_ars_raw, treat_name = "BMP2")
hmsc_adipo_normal <- normalize_hmsc(hmsc_adipo_raw, treat_name = "Adipo")

# 3) Make functions for generating comparison df ====

#Define function
calc_comparison_df <- function(inp_raw){
  #Get unique list of siRNAs/plate names
  sirnas <- inp_raw %>% dplyr::select(siRNA) %>% unique()
  sirnas <- sirnas[,1]
  sirnas <- sirnas[sirnas != "CON"]
  #Calculate p-values for all the sirnas and bind into a dataframe
  sirna_p_list <- lapply(sirnas, calc_comparison_sirna, inp_raw=inp_raw)
  sirna_df <- dplyr::bind_rows(sirna_p_list)
  #Run test for controls
  if ("Treatment" %in% colnames(inp_raw)) {
    treats = unique(inp_raw$Treatment)
    if("Plate" %in% colnames(inp_raw)){
      #Need to test multiple plates worth of controls
      temp <- t(vapply(unique(inp_raw$Plate), calc_treat_test_control, FUN.VALUE = numeric(2), inp_raw=inp_raw))
      make_plate_name <- function(name){paste0("CON", name, "-", treats)}
      temp <- cbind.data.frame(t(vapply(unique(inp_raw$Plate), make_plate_name, FUN.VALUE = character(2))), temp)
      colnames(temp) <- c("group1", "group2", "reps", "p")
      sirna_df <- rbind.data.frame(sirna_df, temp)
    }
    else{
      temp <- c(paste("CON", treats, sep = "-"), calc_treat_test_sirna("CON", inp_raw))
      sirna_df <- rbind.data.frame(sirna_df, temp)
    }
  }
  return(sirna_df)
}
#Calc comparisons for a given siRNA
calc_comparison_sirna <- function(sirna, inp_raw){
  #Filter for the control and given sirna only
  temp <- inp_raw %>% dplyr::filter(siRNA %in% c(sirna, "CON", "Control"))
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
  #Check if treatment conditions present
  if ("Treatment" %in% colnames(temp)) {
    treats <- temp %>% select(Treatment) %>% unique
    treats <- treats[,1]
    #Execute test of comparing siRNA knockdown to control in each treatment state
    treat_p_vals <- t(vapply(treats, FUN = calc_sirna_test_treat, FUN.VALUE = numeric(2), temp=temp, sirna=sirna))
    #Execute test of comparing treated to untreated for the sirna and append to dataframe
    treat_untreat_p_val <- calc_treat_test_sirna(sirna, temp)
    #Create output df for the function after checking for whether plates need to be appended to the controls
    if ("Plate" %in% colnames(temp)) {
      plate = temp %>% dplyr::select("Plate") %>% unique %>% as.character
      sirna_df <- cbind.data.frame(group1 = paste(sirna, treats, sep = "-"), group2 = paste0("CON", plate, "-", treats), reps=treat_p_vals[,1], p=treat_p_vals[,2])
    } else {
      sirna_df <- cbind.data.frame(group1 = paste(sirna, treats, sep = "-"), group2 = paste("CON", treats, sep = "-"), reps=treat_p_vals[,1], p=treat_p_vals[,2])
    }
    sirna_df <- rbind.data.frame(sirna_df, c(paste(sirna, treats, sep = "-"), treat_untreat_p_val))
  } else { #handle event where there is no treated vs. untreated case
    temp2 <- calc_sirna_test_no_treat(sirna = sirna, temp)
    sirna_df <- cbind.data.frame(group1 = sirna, group2 = "CON", reps=temp2[1], p=temp2[2])
  }
  return(sirna_df)
}
#Make functions for executing tests at sirna level when treatment present
calc_sirna_test_treat <- function(treat, temp, sirna){
  temp_two <- temp %>% filter(Treatment == treat)
  control_vals <- temp_two %>% filter(siRNA == "CON") %>% select(Value)
  control_vals <- control_vals[,1]
  sirna_vals <- temp_two %>% filter(siRNA == sirna) %>% select(Value)
  sirna_vals <- sirna_vals[,1]
  test_result <- t.test(sirna_vals, control_vals, alternative = "less")
  return(c(length(sirna_vals), test_result$p.value))
}
calc_treat_test_sirna <- function(sirna, temp){
  temp_two <- temp %>% filter(siRNA == sirna)
  control_vals <- temp_two %>% dplyr::filter(Treatment == "CON") %>% dplyr::select(Value)
  control_vals <- control_vals[,1]
  treat_vals <- temp_two %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(Value)
  treat_vals <- treat_vals[,1]
  test_result <- t.test(control_vals, treat_vals, alternative = "less")
  return(c(length(treat_vals), test_result$p.value))
}
#Make function for testing when no treatment present
calc_sirna_test_no_treat <- function(sirna, temp){
  control_vals <- temp %>% as.data.frame %>% dplyr::filter(siRNA == "CON") %>% dplyr::select(Value)
  control_vals <- control_vals[,1]
  sirna_vals <- temp %>% as.data.frame %>% dplyr::filter(siRNA == sirna) %>% dplyr::select(Value)
  sirna_vals <- sirna_vals[,1]
  test_result <- t.test(sirna_vals, control_vals, alternative = "less")
  return(c(length(sirna_vals), test_result$p.value))
}
#Make function that tests treatment relative to control for control siRNA
calc_treat_test_control <- function(plate, inp_raw){
  temp <- inp_raw %>% dplyr::filter(Plate == plate)
  if ("Passage" %in% colnames(temp)){
    temp <- temp %>% dplyr::group_by(siRNA, Treatment, Donor) %>% dplyr::summarize(Value=mean(Value)) %>% as.data.frame
  }
  return(calc_treat_test_sirna("CON", temp))
}

# 4) Call functions and Correct for Multiple Testing ====

#Call functions
hfob_alp_comparison_df <- calc_comparison_df(hfob_alp_normal)
hmsc_alp_comparison_df <- calc_comparison_df(hmsc_alp_normal)
hmsc_ars_comparison_df <- calc_comparison_df(hmsc_ars_normal)
hmsc_adipo_comparison_df <- calc_comparison_df(hmsc_adipo_normal)

#Correct for multiple testing separately for hfob tests
hfob_alp_comparison_df <- hfob_alp_comparison_df %>% mutate(p.adj = p.adjust(p, method = "BH"))
#Correct for multiple testing separately for the two test types in the hMSC alp experiment
hmsc_alp_comparison_df <- hmsc_alp_comparison_df %>% mutate(p.adj = 1)
hmsc_alp_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_alp_comparison_df$group2)] <-
  p.adjust(hmsc_alp_comparison_df$p[str_detect(pattern = "-CON",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_alp_comparison_df$group2)], method = "BH")
hmsc_alp_comparison_df$p.adj[str_detect(pattern = "-BMP2",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_alp_comparison_df$group2)] <-
  p.adjust(hmsc_alp_comparison_df$p[str_detect(pattern = "-BMP2",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_alp_comparison_df$group2)], method = "BH")
hmsc_alp_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_alp_comparison_df$group2)] <-
  p.adjust(hmsc_alp_comparison_df$p[str_detect(pattern = "-CON",hmsc_alp_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_alp_comparison_df$group2)], method = "BH")
#Correct for multiple testing separately for the two test types in the hMSC ars experiment
hmsc_ars_comparison_df <- hmsc_ars_comparison_df %>% mutate(p.adj = 1)
hmsc_ars_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_ars_comparison_df$group2)] <-
  p.adjust(hmsc_ars_comparison_df$p[str_detect(pattern = "-CON",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_ars_comparison_df$group2)], method = "BH")
hmsc_ars_comparison_df$p.adj[str_detect(pattern = "-BMP2",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_ars_comparison_df$group2)] <-
  p.adjust(hmsc_ars_comparison_df$p[str_detect(pattern = "-BMP2",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-BMP2", hmsc_ars_comparison_df$group2)], method = "BH")
hmsc_ars_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_ars_comparison_df$group2)] <-
  p.adjust(hmsc_ars_comparison_df$p[str_detect(pattern = "-CON",hmsc_ars_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_ars_comparison_df$group2)], method = "BH")
#Correct for multiple testing separately for the two test types in the hMSC alp experiment
hmsc_adipo_comparison_df <- hmsc_adipo_comparison_df %>% mutate(p.adj = 1)
hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2)] <-
  p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2)], method = "BH")
hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-Adipo",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2)] <-
  p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-Adipo",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-Adipo", hmsc_adipo_comparison_df$group2)], method = "BH")
hmsc_adipo_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_adipo_comparison_df$group2)] <-
  p.adjust(hmsc_adipo_comparison_df$p[str_detect(pattern = "-CON",hmsc_adipo_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_adipo_comparison_df$group2)], method = "BH")

# 5) Make Box Plots of Results ====

#Get order of names for plots
name_order <- c("CONTROL-1-BMP2", "CONTROL-1-CON", "CONTROL-2-BMP2", "CONTROL-2-CON", sort(unique(as.character(hmsc_alp_normal$Name[str_detect(pattern = "CON-", hmsc_alp_normal$Name) == FALSE]))))
adipo_name_order <- str_replace(name_order, pattern = "-BMP2", "-Adipo")
adipo_name_order <- str_replace(adipo_name_order[str_detect(adipo_name_order, pattern = "CONTROL-2", negate = TRUE)], pattern = "CONTROL-1", replacement = "CONTROL")
sirna_order <- sort(unique(as.character(hmsc_alp_normal$siRNA[str_detect(pattern = "CON", hmsc_alp_normal$Name) == FALSE])))

### hFOB ALP ###
#Combine the measurements together and remove the control siRNA
hfob_alp_combined <- hfob_alp_normal %>% dplyr::group_by(siRNA, Replicate) %>% dplyr::summarize(Value = mean(Value)) %>% filter(siRNA != "CON")
#Add on y-positions to the signficance df
hfob_alp_y_pos_df <- hfob_alp_combined %>% dplyr::group_by(siRNA) %>% dplyr::summarize(max_value = max(Value))
hfob_alp_comparison_df <- hfob_alp_comparison_df %>% mutate(siRNA = group1) %>% inner_join(hfob_alp_y_pos_df, by = "siRNA") %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>% 
  mutate(max_value_padded = max_value + 0.04, signif = p.adj < 0.05)
#Add significance back onto the plotting df
hfob_alp_combined <- hfob_alp_combined %>% inner_join(hfob_alp_comparison_df, by = "siRNA")
#Make plot
hfob_alp_plot <- ggplot(hfob_alp_combined, aes(x=siRNA, y=Value,color = signif)) + geom_boxplot(width=0.3) +
                  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 2, alpha = 0.6) + 
                  stat_pvalue_manual(hfob_alp_comparison_df, x= "siRNA", y = "max_value_padded", label = "p.signif", bracket.shorten = -0.1, label.size = 6, color = "signif") + 
                  scale_color_manual(values = c("gray", "dodgerblue3")) + 
                  scale_x_discrete(labels = sirna_order, breaks = sirna_order) + 
                  ylab("hFOB ALP") + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
                        axis.title=element_text(size=16), 
                        legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 

### hMSC ALP ###
#Make a supplemental figure showing effect of BMP2 stimulation
#Split out the controls onto separate sirnas
hmsc_alp_plot <- hmsc_alp_normal %>% mutate(siRNA = ifelse(siRNA == "CON", paste0("CONTROL-", Plate), siRNA),
                                            Name = ifelse(str_detect(pattern = "CON-", Name), str_replace(Name, pattern = "CON-", replacement = paste0("CONTROL-", Plate, "-")), Name))
#Set axis plotting order
hmsc_alp_plot$Name <- factor(hmsc_alp_plot$Name, levels = name_order)
#Add on y positions for significance df and convert to asterisks
hmsc_alp_y_pos_df <- hmsc_alp_plot %>% dplyr::group_by(siRNA) %>% dplyr::summarize(max_value = max(Value))
hmsc_alp_comparison_df_filt <- hmsc_alp_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-CON") & str_detect(group2, pattern = "-BMP2")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(group1 = str_replace(group1, "^CON", "CONTROL-")) %>% mutate(group2 = str_replace(group2, "^CON", "CONTROL-")) %>%
  mutate(siRNA = str_replace(group1, "-CON", "")) %>%
  inner_join(hmsc_alp_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1) #Prep comparison df for plotting
#Plot the figure
jpeg(paste0(out_dir, "treatment_v_control.hmsc_alp.box.jpeg"), width = 12000, height=6400, res=1000)
ggplot(hmsc_alp_plot, aes(x=Name, y=Value, color=Treatment)) + geom_boxplot() + 
  stat_pvalue_manual(hmsc_alp_comparison_df_filt, label = "p.signif", bracket.shorten = -0.1, label.size = 6) + 
  ylab("hMSC ALP") + 
  scale_color_manual(values = c("dodgerblue3","gray")) + 
  scale_x_discrete(labels = ifelse(str_detect(name_order, pattern = "-CON"), "", str_replace(name_order, pattern = "-BMP2", replacement = ""))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=45, hjust = 0.85, vjust = 0.9), 
        axis.title=element_text(size=30), legend.key = element_rect(fill = "transparent"), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=24), axis.title.x = element_blank()) + 
  ggplot2::theme(legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 24)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))
dev.off()

#Make the subpanel plot for the main figure now
hmsc_alp_main_plot <- hmsc_alp_plot %>% dplyr::filter(!(str_detect(siRNA, pattern = "CONTROL")) & Treatment != "CON")
#Now filter the comparison df correctly
hmsc_alp_comparison_df_main <- hmsc_alp_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-BMP2") & 
                                                                          str_detect(group2, pattern = "-BMP2") & 
                                                                          str_detect(group2, pattern = "CON")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(siRNA = str_replace(group1, pattern="-BMP2", replacement="")) %>%
  inner_join(hmsc_alp_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1, signif = p.adj < 0.05) %>% 
  mutate(max_value_padded = max_value + 0.04)#Prep comparison df for plotting
#Add the significance info back onto the main plotting dataframe
hmsc_alp_main_plot <- hmsc_alp_main_plot %>% inner_join(hmsc_alp_comparison_df_main, by = "siRNA")
#Make plot
hmsc_alp_plot <- ggplot(hmsc_alp_main_plot, aes(x=siRNA, y=Value,color = signif)) + geom_boxplot(width=0.3) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 2, alpha = 0.6) + 
  stat_pvalue_manual(hmsc_alp_comparison_df_main, x= "siRNA", y = "max_value_padded", label = "p.signif", bracket.shorten = -0.1, label.size = 6, color = "signif") + 
  scale_color_manual(values = c("gray", "dodgerblue3")) + 
  scale_x_discrete(labels = sirna_order, breaks = sirna_order) + 
  ylab("hMSC ALP") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
        axis.title=element_text(size=16), 
        legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 


### hMSC Alizarin ###
#Make a supplemental figure showing effect of BMP2 stimulation
#Split out the controls onto separate sirnas
hmsc_ars_plot <- hmsc_ars_normal %>% mutate(siRNA = ifelse(siRNA == "CON", paste0("CONTROL-", Plate), siRNA),
                                            Name = ifelse(str_detect(pattern = "CON-", Name), str_replace(Name, pattern = "CON-", replacement = paste0("CONTROL-", Plate, "-")), Name))
#Set axis plotting order
hmsc_ars_plot$Name <- factor(hmsc_ars_plot$Name, levels = name_order)
#Add on y positions for significance df and convert to asterisks
hmsc_ars_y_pos_df <- hmsc_ars_plot %>% dplyr::group_by(siRNA) %>% dplyr::summarize(max_value = max(Value))
hmsc_ars_comparison_df_filt <- hmsc_ars_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-CON") & str_detect(group2, pattern = "-BMP2")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(group1 = str_replace(group1, "^CON", "CONTROL-")) %>% mutate(group2 = str_replace(group2, "^CON", "CONTROL-")) %>%
  mutate(siRNA = str_replace(group1, "-CON", "")) %>%
  inner_join(hmsc_ars_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1) #Prep comparison df for plotting
#Plot the figure
jpeg(paste0(out_dir, "treatment_v_control.hmsc_ars.box.jpeg"), width = 12000, height=6400, res=1000)
ggplot(hmsc_ars_plot, aes(x=Name, y=Value, color=Treatment)) + geom_boxplot() + 
  stat_pvalue_manual(hmsc_ars_comparison_df_filt, label = "p.signif", bracket.shorten = -0.1, label.size = 6) + 
  ylab("hMSC ARS") + 
  scale_color_manual(values = c("dodgerblue3","gray")) + 
  scale_x_discrete(labels = ifelse(str_detect(name_order, pattern = "-CON"), "", str_replace(name_order, pattern = "-BMP2", replacement = ""))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=45, hjust = 0.85, vjust = 0.9), 
        axis.title=element_text(size=30), legend.key = element_rect(fill = "transparent"), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=24), axis.title.x = element_blank()) + 
  ggplot2::theme(legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 24)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))
dev.off()

#Make the subpanel plot for the main figure now
hmsc_ars_main_plot <- hmsc_ars_plot %>% dplyr::filter(!(str_detect(siRNA, pattern = "CONTROL")) & Treatment != "CON")
#Now filter the comparison df correctly
hmsc_ars_comparison_df_main <- hmsc_ars_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-BMP2") & 
                                                                          str_detect(group2, pattern = "-BMP2") & 
                                                                          str_detect(group2, pattern = "CON")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(siRNA = str_replace(group1, pattern="-BMP2", replacement="")) %>%
  inner_join(hmsc_ars_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1, signif = p.adj < 0.05) %>% 
  mutate(max_value_padded = max_value + 0.04)#Prep comparison df for plotting
#Add the significance info back onto the main plotting dataframe
hmsc_ars_main_plot <- hmsc_ars_main_plot %>% inner_join(hmsc_ars_comparison_df_main, by = "siRNA")
#Make plot
hmsc_ars_plot <- ggplot(hmsc_ars_main_plot, aes(x=siRNA, y=Value,color = signif)) + geom_boxplot(width=0.3) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 2, alpha = 0.6) + 
  stat_pvalue_manual(hmsc_ars_comparison_df_main, x= "siRNA", y = "max_value_padded", label = "p.signif", bracket.shorten = -0.1, label.size = 6, color = "signif") + 
  scale_color_manual(values = c("gray", "dodgerblue3")) + 
  scale_x_discrete(labels = sirna_order, breaks = sirna_order) + 
  ylab("hMSC ARS") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
        axis.title=element_text(size=16), 
        legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 

### hMSC Adipogenesis ###
#Split out the controls onto separate sirnas
hmsc_adipo_plot <- hmsc_adipo_normal %>% mutate(siRNA = ifelse(siRNA == "CON", "CONTROL", siRNA),
                                            Name = ifelse(str_detect(pattern = "CON-", Name), str_replace(Name, pattern = "CON-", replacement = "CONTROL-"), Name))
#Set axis plotting order
hmsc_adipo_plot$Name <- factor(hmsc_adipo_plot$Name, levels = adipo_name_order)
#Add on y positions for significance df and convert to asterisks
hmsc_adipo_y_pos_df <- hmsc_adipo_plot %>% dplyr::group_by(siRNA) %>% dplyr::summarize(max_value = max(Value))
hmsc_adipo_comparison_df_filt <- hmsc_adipo_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-CON") & str_detect(group2, pattern = "-Adipo")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(group1 = str_replace(group1, "^CON", "CONTROL-")) %>% mutate(group2 = str_replace(group2, "^CON", "CONTROL-")) %>%
  mutate(siRNA = str_replace(group1, "-CON", "")) %>%
  inner_join(hmsc_adipo_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1) #Prep comparison df for plotting
#Plot the figure
jpeg(paste0(out_dir, "treatment_v_control.hmsc_adipo.box.jpeg"), width = 12000, height=6400, res=1000)
ggplot(hmsc_adipo_plot, aes(x=Name, y=Value, color=Treatment)) + geom_boxplot() + 
  stat_pvalue_manual(hmsc_adipo_comparison_df_filt, label = "p.signif", bracket.shorten = -0.1, label.size = 6) + 
  ylab("hMSC Adipogenesis") + 
  scale_color_manual(values = c("dodgerblue3","gray")) + 
  scale_x_discrete(labels = ifelse(str_detect(adipo_name_order, pattern = "-CON"), "", str_replace(adipo_name_order, pattern = "-Adipo", replacement = ""))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=45, hjust = 0.6, vjust = 0.7), 
        axis.title=element_text(size=30), legend.key = element_rect(fill = "transparent"), 
        legend.position = "bottom", axis.text.y = element_text(color="black", size=24), axis.title.x = element_blank()) + 
  ggplot2::theme(legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 24)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))
dev.off()

#Make the subpanel plot for the main figure now
hmsc_adipo_main_plot <- hmsc_adipo_plot %>% dplyr::filter(!(str_detect(siRNA, pattern = "CONTROL")) & Treatment != "CON")
#Now filter the comparison df correctly
hmsc_adipo_comparison_df_main <- hmsc_adipo_comparison_df %>% dplyr::filter(str_detect(group1, pattern = "-Adipo") & 
                                                                          str_detect(group2, pattern = "-Adipo") & 
                                                                          str_detect(group2, pattern = "CON")) %>%
  mutate(p.signif = ifelse(p.adj < 0.05, "*", "")) %>%
  mutate(siRNA = str_replace(group1, pattern="-Adipo", replacement="")) %>%
  inner_join(hmsc_adipo_y_pos_df, by = "siRNA") %>% 
  mutate(y.position = max_value + 0.1, signif = p.adj < 0.05) %>% 
  mutate(max_value_padded = max_value + 0.04)#Prep comparison df for plotting
#Add the significance info back onto the main plotting dataframe
hmsc_adipo_main_plot <- hmsc_adipo_main_plot %>% inner_join(hmsc_adipo_comparison_df_main, by = "siRNA")
#Set factor levels
hmsc_adipo_main_plot$siRNA <- factor(hmsc_adipo_main_plot$siRNA, levels = sirna_order)
hmsc_adipo_comparison_df_main$siRNA <- factor(hmsc_adipo_comparison_df_main$siRNA, levels = sirna_order)
#Make plot
hmsc_adipo_plot <- ggplot(hmsc_adipo_main_plot, aes(x=siRNA, y=Value,color = signif)) + geom_boxplot(width=0.3) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 2, alpha = 0.6) + 
  stat_pvalue_manual(hmsc_adipo_comparison_df_main, x= "siRNA", y = "max_value_padded", label = "p.signif", bracket.shorten = -0.1, label.size = 6, color = "signif") + 
  scale_color_manual(values = c("gray", "dodgerblue3")) + 
  scale_x_discrete(labels = sirna_order, breaks = sirna_order) + 
  ylab("hMSC Adipogenesis") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
        axis.title=element_text(size=16), 
        legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 

#Make stacked combined plot data frames
combined_hfob_alp_plot <- hfob_alp_combined %>% dplyr::select(siRNA, signif, reps, Value) %>% mutate(reps, subplot = "hFOB ALP")
combined_hmsc_alp_plot <- hmsc_alp_main_plot %>% dplyr::select(siRNA, signif, reps, Value) %>% mutate(subplot = "hMSC ALP")
combined_hmsc_ars_plot <- hmsc_ars_main_plot %>% dplyr::select(siRNA, signif, reps, Value) %>% mutate(subplot = "hMSC ARS")
combined_hmsc_adipo_plot <- hmsc_adipo_main_plot %>% dplyr::select(siRNA, signif, reps, Value) %>% mutate(subplot = "hMSC Adipo")
combined_plot_df <- rbind.data.frame(combined_hfob_alp_plot, combined_hmsc_alp_plot, combined_hmsc_ars_plot, combined_hmsc_adipo_plot)
combined_plot_df$subplot <- factor(combined_plot_df$subplot, levels = c("hFOB ALP", "hMSC ALP", "hMSC ARS", "hMSC Adipo"))
#And now make significance data frames
combined_hfob_alp_comparison <- hfob_alp_comparison_df %>% dplyr::select(group1, group2, siRNA, max_value_padded, p.signif, signif) %>% mutate(subplot = "hFOB ALP")
combined_hmsc_alp_comparison <- hmsc_alp_comparison_df_main %>% dplyr::select(group1, group2, siRNA, max_value_padded, p.signif, signif) %>% mutate(subplot = "hMSC ALP")
combined_hmsc_ars_comparison <- hmsc_ars_comparison_df_main %>% dplyr::select(group1, group2, siRNA, max_value_padded, p.signif, signif) %>% mutate(subplot = "hMSC ARS")
combined_hmsc_adipo_comparison <- hmsc_adipo_comparison_df_main %>% dplyr::select(group1, group2, siRNA, max_value_padded, p.signif, signif) %>% mutate(subplot = "hMSC Adipo")
combined_comparison_df <- rbind.data.frame(combined_hfob_alp_comparison, combined_hmsc_alp_comparison, combined_hmsc_ars_comparison, combined_hmsc_adipo_comparison)
combined_comparison_df$subplot <- factor(combined_comparison_df$subplot, levels = c("hFOB ALP", "hMSC ALP", "hMSC ARS", "hMSC Adipo"))

#Make combined plot
jpeg(paste0(out_dir, "main_assay_figure.box.jpeg"), width = 12000, height=12000, res=1000)
ggplot(combined_plot_df, aes(x = siRNA, y = Value, color=signif)) +
  geom_boxplot(linewidth=1.03, width=0.5) + 
#  stat_pvalue_manual(combined_comparison_df, x= "siRNA", y = "max_value_padded", label = "p.signif", label.size = 6, color = "signif", show.legend=FALSE) + 
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
  facet_grid(subplot ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c("gray", "dodgerblue3"), labels=c("Not Significant", "Adj. P-Value < 0.05")) + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
        axis.title = element_blank(), legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
        legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
        legend.title = ggplot2::element_blank(),
        legend.text = element_text(size = 20), legend.key = element_blank(),
        axis.text.y = element_text(size = 20),
        strip.text.y = element_text(size = 20)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))
dev.off()

#Make combined plot in tiff format
tiff(paste0(out_dir, "main_assay_figure.box.tiff"), width = 12000, height=12000, res=1000)
ggplot(combined_plot_df, aes(x = siRNA, y = Value, color=signif)) +
  geom_boxplot(linewidth=1.03, width=0.5) + 
  #  stat_pvalue_manual(combined_comparison_df, x= "siRNA", y = "max_value_padded", label = "p.signif", label.size = 6, color = "signif", show.legend=FALSE) + 
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
  facet_grid(subplot ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c("gray", "dodgerblue3"), labels=c("Not Significant", "Adj. P-Value < 0.05")) + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
                 axis.title = element_blank(), legend.position = "bottom",
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
                 legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 20), legend.key = element_blank(),
                 axis.text.y = element_text(size = 20),
                 strip.text.y = element_text(size = 20)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))
dev.off()

# 6) Make Supplemental Tables of Results ====

### Start with hfob alp ###
hfob_alp_supplement <- combined_hfob_alp_plot %>% dplyr::group_by(siRNA) %>% 
  dplyr::summarize(mean_fold_change=mean(Value)) %>% 
  inner_join(hfob_alp_comparison_df, by = "siRNA") %>% 
  dplyr::select(siRNA, reps, mean_fold_change, p, p.adj)
write.table(hfob_alp_supplement, file = paste0(inp_dir, "hfob_alp.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### hMSC alp ###
hmsc_alp_supplement <- combined_hmsc_alp_plot %>% dplyr::group_by(siRNA) %>% 
  dplyr::summarize(mean_fold_change=mean(Value)) %>% 
  inner_join(hmsc_alp_comparison_df_main, by = "siRNA") %>% 
  dplyr::select(siRNA, reps, mean_fold_change, p, p.adj)
write.table(hmsc_alp_supplement, file = paste0(inp_dir, "hmsc_alp.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### hMSC alizarin ###
hmsc_ars_supplement <- combined_hmsc_ars_plot %>% dplyr::group_by(siRNA) %>% 
  dplyr::summarize(mean_fold_change=mean(Value)) %>% 
  inner_join(hmsc_ars_comparison_df_main, by = "siRNA") %>% 
  dplyr::select(siRNA, reps, mean_fold_change, p, p.adj)
write.table(hmsc_ars_supplement, file = paste0(inp_dir, "hmsc_ars.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### hMSC adipo###
hmsc_adipo_supplement <- combined_hmsc_adipo_plot %>% dplyr::group_by(siRNA) %>% 
  dplyr::summarize(mean_fold_change=mean(Value)) %>% 
  inner_join(hmsc_adipo_comparison_df_main, by = "siRNA") %>% 
  dplyr::select(siRNA, reps, mean_fold_change, p, p.adj)
write.table(hmsc_adipo_supplement, file = paste0(inp_dir, "hmsc_adipo.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


