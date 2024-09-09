################################################################################

#make_plate-normalized_assay_figures.R

#The purpose of this script is to make figures and tables for the osteoblast
#cell model assay results. This script uses the plate-normalized intensity 
#measurements to compare siRNAs. This script runs locally.

################################################################################


# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(tidyverse)
library(ggpubr)

#Set directory locations
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/osteoblast_assays/"
out_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/osteoblast_assays/"

#Set assay result locations
hfob_alp_file <- paste0(inp_dir, "hfob_alp_results.tsv")
hmsc_alp_file <- paste0(inp_dir, "hMSC-BMP2_alp_results.tsv")
hmsc_ars_file <- paste0(inp_dir, "hMSC-BMP2_ars_results.tsv")
hmsc_adipo_file <- paste0(inp_dir, "hMSC-BMP2_adipo_results.tsv")

#Set assay cell count locations
hfob_alp_count_file <- paste0(inp_dir, "hFOB_alp_cell_counts.tsv")
hmsc_alp_count_file <- paste0(inp_dir, "hMSC-BMP2_alp_cell_counts.tsv")
hmsc_ars_count_file <- paste0(inp_dir, "hMSC-BMP2_ars_cell_counts.tsv")
hmsc_adipo_count_file <- paste0(inp_dir, "hMSC-ADIPO_adipo_cell_counts.tsv")

#Set random seed
set.seed(5)

# 1) Read in assay files ====

#Read in assay files
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

#Combine different measurments for hFOB ALP assay wells together by mean
hfob_alp_raw <- hfob_alp_raw %>% group_by(siRNA, Replicate, Name) %>% dplyr::summarise(Value = mean(Value)) %>% as.data.frame()

#Remove EPDR1
hfob_alp_raw <- hfob_alp_raw %>% dplyr::filter(siRNA != "EPDR1")
hmsc_alp_raw <- hmsc_alp_raw %>% dplyr::filter(siRNA != "EPDR1")
hmsc_ars_raw <- hmsc_ars_raw %>% dplyr::filter(siRNA != "EPDR1")
hmsc_adipo_raw <- hmsc_adipo_raw %>% dplyr::filter(siRNA != "EPDR1")

# 2) Read in Cell Count files ====

#Read in cell count files
hfob_alp_count_raw <- read.csv(hfob_alp_count_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_alp_count_raw <- read.csv(hmsc_alp_count_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_ars_count_raw <- read.csv(hmsc_ars_count_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
hmsc_adipo_count_raw <- read.csv(hmsc_adipo_count_file, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#Process the hfob data set to keep formatting consistent
#Convert all "Control"s to "CON"s
hfob_alp_count_raw$siRNA = ifelse(hfob_alp_count_raw$siRNA == "Control", "CON", hfob_alp_count_raw$siRNA)
#Remove no siRNA test from hFOB data
hfob_alp_count_raw <- hfob_alp_count_raw %>% dplyr::filter(siRNA != "NT" & siRNA != "Not Treated")

#Duplicate over the ARS control to BMP2
temp <- hmsc_ars_count_raw
temp$Treatment <- "BMP2"
temp$Name <- paste(temp$siRNA, temp$Treatment, sep = "-")
hmsc_ars_count_raw <- rbind.data.frame(hmsc_ars_count_raw, temp)

# 3) Match Up Replicates and Create DAPI-Normalized Intensity Measures ====

### Start with hFOB ALP ###
#Add rownames onto the dataframes that can be used for filtering
rownames(hfob_alp_raw) <- with(hfob_alp_raw, paste(siRNA, Replicate, sep="."))
rownames(hfob_alp_count_raw) <- with(hfob_alp_count_raw, paste(siRNA, Replicate, sep="."))
hfob_alp_intersect = intersect(rownames(hfob_alp_raw), rownames(hfob_alp_count_raw))
#Filter down the dataframes
hfob_alp_filt <- hfob_alp_raw[hfob_alp_intersect,]
hfob_alp_count_filt <- hfob_alp_count_raw[hfob_alp_intersect,]
#Make the dataframe of intensity normalized by cell count
temp <- intersect(colnames(hfob_alp_filt), colnames(hfob_alp_count_filt))
temp <- temp[temp != "Value"]
hfob_alp_intensity_dapi <- cbind.data.frame(hfob_alp_filt[,temp], Value = hfob_alp_filt$Value/hfob_alp_count_filt$Value)

### Next hMSC ALP ###
#Add rownames onto the dataframes that can be used for filtering
rownames(hmsc_alp_raw) <- with(hmsc_alp_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
rownames(hmsc_alp_count_raw) <- with(hmsc_alp_count_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
hmsc_alp_intersect = intersect(rownames(hmsc_alp_raw), rownames(hmsc_alp_count_raw))
#Filter down the dataframes
hmsc_alp_filt <- hmsc_alp_raw[hmsc_alp_intersect,]
hmsc_alp_count_filt <- hmsc_alp_count_raw[hmsc_alp_intersect,]
#Make the dataframe of intensity normalized by cell count
temp <- intersect(colnames(hmsc_alp_filt), colnames(hmsc_alp_count_filt))
temp <- temp[temp != "Value"]
hmsc_alp_intensity_dapi <- cbind.data.frame(hmsc_alp_filt[,temp], Value = hmsc_alp_filt$Value/hmsc_alp_count_filt$Value)

### Next hMSC ARS ###
#Add rownames onto the dataframes that can be used for filtering
rownames(hmsc_ars_raw) <- with(hmsc_ars_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
rownames(hmsc_ars_count_raw) <- with(hmsc_ars_count_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
hmsc_ars_intersect = intersect(rownames(hmsc_ars_raw), rownames(hmsc_ars_count_raw))
#Filter down the dataframes
hmsc_ars_filt <- hmsc_ars_raw[hmsc_ars_intersect,]
hmsc_ars_count_filt <- hmsc_ars_count_raw[hmsc_ars_intersect,]
#Make the dataframe of intensity normalized by cell count
temp <- intersect(colnames(hmsc_ars_filt), colnames(hmsc_ars_count_filt))
temp <- temp[temp != "Value"]
hmsc_ars_intensity_dapi <- cbind.data.frame(hmsc_ars_filt[,temp], Value = hmsc_ars_filt$Value/hmsc_ars_count_filt$Value)

### Last hMSC Adipo ###
#Add rownames onto the dataframes that can be used for filtering
rownames(hmsc_adipo_raw) <- with(hmsc_adipo_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
rownames(hmsc_adipo_count_raw) <- with(hmsc_adipo_count_raw, paste(Name, siRNA, Treatment, Plate, Donor, Passage, sep="."))
hmsc_adipo_intersect = intersect(rownames(hmsc_adipo_raw), rownames(hmsc_adipo_count_raw))
#Filter down the dataframes
hmsc_adipo_filt <- hmsc_adipo_raw[hmsc_adipo_intersect,]
hmsc_adipo_count_filt <- hmsc_adipo_count_raw[hmsc_adipo_intersect,]
#Make the dataframe of intensity normalized by cell count
temp <- intersect(colnames(hmsc_adipo_filt), colnames(hmsc_adipo_count_filt))
temp <- temp[temp != "Value"]
hmsc_adipo_intensity_dapi <- cbind.data.frame(hmsc_adipo_filt[,temp], Value = hmsc_adipo_filt$Value/hmsc_adipo_count_filt$Value)

# 4) Normalize Data to Control Level of Each Plate ====

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
hfob_alp_normal <- normalize_hfob(hfob_alp_filt)
hfob_alp_count_normal <- normalize_hfob(hfob_alp_count_filt)
hfob_alp_intensity_dapi_normal <- normalize_hfob(hfob_alp_intensity_dapi)

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
#Call the function on the hmsc alp results
hmsc_alp_normal <- normalize_hmsc(hmsc_alp_filt, treat_name = "BMP2")
hmsc_alp_count_normal <- normalize_hmsc(hmsc_alp_count_filt, treat_name = "BMP2")
hmsc_alp_intensity_dapi_normal <- normalize_hmsc(hmsc_alp_intensity_dapi, treat_name = "BMP2")
#Call the function on the hmsc ars results
hmsc_ars_normal <- normalize_hmsc(hmsc_ars_filt, treat_name = "BMP2")
hmsc_ars_count_normal <- normalize_hmsc(hmsc_ars_count_filt, treat_name = "BMP2")
hmsc_ars_intensity_dapi_normal <- normalize_hmsc(hmsc_ars_intensity_dapi, treat_name = "BMP2")
#Call the function on the hmsc adipogenesis results
hmsc_adipo_normal <- normalize_hmsc(hmsc_adipo_filt, treat_name = "Adipo")
hmsc_adipo_count_normal <- normalize_hmsc(hmsc_adipo_count_filt, treat_name = "Adipo")
hmsc_adipo_intensity_dapi_normal <- normalize_hmsc(hmsc_adipo_intensity_dapi, treat_name = "Adipo")

# 5) Make functions for generating comparison df ====

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
  #Filter for donor controls with matched sirna treatment
  sirna_donors <- temp_two %>% dplyr::filter(siRNA == sirna) %>% dplyr::select(Donor)
  sirna_donors <- sirna_donors[,1]
  temp_two <- temp_two %>% dplyr::filter(Donor %in% sirna_donors) %>% dplyr::arrange(Donor)
  control_vals <- temp_two %>% filter(siRNA == "CON") %>% select(Value)
  control_vals <- control_vals[,1]
  sirna_vals <- temp_two %>% filter(siRNA == sirna) %>% select(Value)
  sirna_vals <- sirna_vals[,1]
  test_result <- t.test(sirna_vals, control_vals, alternative = "less", paired = TRUE)
  return(c(length(sirna_vals), test_result$p.value))
}
calc_treat_test_sirna <- function(sirna, temp){
  temp_two <- temp %>% filter(siRNA == sirna)
  #Filter for donor controls with matched sirna treatment
  sirna_donors <- temp_two %>% dplyr::filter(siRNA == sirna) %>% dplyr::select(Donor)
  sirna_donors <- sirna_donors[,1]
  temp_two <- temp_two %>% dplyr::filter(Donor %in% sirna_donors) %>% dplyr::arrange(Donor)
  control_vals <- temp_two %>% dplyr::filter(Treatment == "CON") %>% dplyr::select(Value)
  control_vals <- control_vals[,1]
  treat_vals <- temp_two %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(Value)
  treat_vals <- treat_vals[,1]
  test_result <- t.test(control_vals, treat_vals, alternative = "less", paired = TRUE)
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

# 6) Call functions and Correct for Multiple Testing ====

#Call functions for raw intensity data
hfob_alp_comparison_df <- calc_comparison_df(hfob_alp_normal)
hmsc_alp_comparison_df <- calc_comparison_df(hmsc_alp_normal)
hmsc_ars_comparison_df <- calc_comparison_df(hmsc_ars_normal)
hmsc_adipo_comparison_df <- calc_comparison_df(hmsc_adipo_normal)

#Call functions for raw intensity data
hfob_alp_count_comparison_df <- calc_comparison_df(hfob_alp_count_normal)
hmsc_alp_count_comparison_df <- calc_comparison_df(hmsc_alp_count_normal)
hmsc_ars_count_comparison_df <- calc_comparison_df(hmsc_ars_count_normal)
hmsc_adipo_count_comparison_df <- calc_comparison_df(hmsc_adipo_count_normal)

#Call functions for intensity data normalized by DAPI cell counts
hfob_alp_intensity_dapi_comparison_df <- calc_comparison_df(hfob_alp_intensity_dapi_normal)
hmsc_alp_intensity_dapi_comparison_df <- calc_comparison_df(hmsc_alp_intensity_dapi_normal)
hmsc_ars_intensity_dapi_comparison_df <- calc_comparison_df(hmsc_ars_intensity_dapi_normal)
hmsc_adipo_intensity_dapi_comparison_df <- calc_comparison_df(hmsc_adipo_intensity_dapi_normal)

#Remove results with fewer than three replicates
remove_few_rep_results <- function(comparison_df, rep_cutoff=3){
  temp <- comparison_df %>% dplyr::mutate(p = ifelse(reps < rep_cutoff, NA, p))
  return(temp)
}
hfob_alp_comparison_df <- remove_few_rep_results(hfob_alp_comparison_df)
hmsc_alp_comparison_df <- remove_few_rep_results(hmsc_alp_comparison_df)
hmsc_ars_comparison_df <- remove_few_rep_results(hmsc_ars_comparison_df)
hmsc_adipo_comparison_df <- remove_few_rep_results(hmsc_adipo_comparison_df)
hfob_alp_count_comparison_df <- remove_few_rep_results(hfob_alp_count_comparison_df)
hmsc_alp_count_comparison_df <- remove_few_rep_results(hmsc_alp_count_comparison_df)
hmsc_ars_count_comparison_df <- remove_few_rep_results(hmsc_ars_count_comparison_df)
hmsc_adipo_count_comparison_df <- remove_few_rep_results(hmsc_adipo_count_comparison_df)
hfob_alp_intensity_dapi_comparison_df <- remove_few_rep_results(hfob_alp_intensity_dapi_comparison_df)
hmsc_alp_intensity_dapi_comparison_df <- remove_few_rep_results(hmsc_alp_intensity_dapi_comparison_df)
hmsc_ars_intensity_dapi_comparison_df <- remove_few_rep_results(hmsc_ars_intensity_dapi_comparison_df)
hmsc_adipo_intensity_dapi_comparison_df <- remove_few_rep_results(hmsc_adipo_intensity_dapi_comparison_df)

#Correct for multiple testing separately for hfob tests
hfob_alp_comparison_df <- hfob_alp_comparison_df %>% mutate(p.adj = p.adjust(p, method = "BH"))
hfob_alp_count_comparison_df <- hfob_alp_count_comparison_df %>% mutate(p.adj = p.adjust(p, method = "BH"))
hfob_alp_intensity_dapi_comparison_df <- hfob_alp_intensity_dapi_comparison_df %>% mutate(p.adj = p.adjust(p, method = "BH"))

#Make a function to correct for multiple testing separately for the different test types in the hMSC experiments
correct_hmsc_bh <- function(hmsc_comparison_df, treat = "BMP2"){
  treat = paste0("-", treat)
  hmsc_comparison_df <- hmsc_comparison_df %>% mutate(p.adj = 1)
  hmsc_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_comparison_df$group1) & str_detect(pattern = treat, hmsc_comparison_df$group2)] <-
    p.adjust(hmsc_comparison_df$p[str_detect(pattern = "-CON",hmsc_comparison_df$group1) & str_detect(pattern = treat, hmsc_comparison_df$group2)], method = "BH")
  hmsc_comparison_df$p.adj[str_detect(pattern = treat,hmsc_comparison_df$group1) & str_detect(pattern = treat, hmsc_comparison_df$group2)] <-
    p.adjust(hmsc_comparison_df$p[str_detect(pattern = treat,hmsc_comparison_df$group1) & str_detect(pattern = treat, hmsc_comparison_df$group2)], method = "BH")
  hmsc_comparison_df$p.adj[str_detect(pattern = "-CON",hmsc_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_comparison_df$group2)] <-
    p.adjust(hmsc_comparison_df$p[str_detect(pattern = "-CON",hmsc_comparison_df$group1) & str_detect(pattern = "-CON", hmsc_comparison_df$group2)], method = "BH")
  return(hmsc_comparison_df)
}
#Call the function on the hMSC ALP Results
hmsc_alp_comparison_df <- correct_hmsc_bh(hmsc_alp_comparison_df)
hmsc_alp_count_comparison_df <- correct_hmsc_bh(hmsc_alp_count_comparison_df)
hmsc_alp_intensity_dapi_comparison_df <- correct_hmsc_bh(hmsc_alp_intensity_dapi_comparison_df)
#Call the function on the hMSC ARS Results
hmsc_ars_comparison_df <- correct_hmsc_bh(hmsc_ars_comparison_df)
hmsc_ars_count_comparison_df <- correct_hmsc_bh(hmsc_ars_count_comparison_df)
hmsc_ars_intensity_dapi_comparison_df <- correct_hmsc_bh(hmsc_ars_intensity_dapi_comparison_df)
#Call the function on the hMSC Adipo Results
hmsc_adipo_comparison_df <- correct_hmsc_bh(hmsc_adipo_comparison_df, treat = "Adipo")
hmsc_adipo_count_comparison_df <- correct_hmsc_bh(hmsc_adipo_count_comparison_df, treat = "Adipo")
hmsc_adipo_intensity_dapi_comparison_df <- correct_hmsc_bh(hmsc_adipo_intensity_dapi_comparison_df, treat = "Adipo")

# 7) Append on significance status to Data Frames ====

#Create a function that appends on significance 
append_signif <- function(plot_df, comparison_df){
  #Check to see if we have hMSC or hFOB Data
  if (str_detect(comparison_df$group1[1],"-") & str_detect(comparison_df$group2[1],"-")) {
    comparison_df <- comparison_df %>% dplyr::filter(!(str_detect(group1, "-CON")) & !(str_detect(group2, "-CON"))) %>%
      dplyr::mutate(group1 = sub("-.*", "", group1)) %>% dplyr::mutate(group2 = sub("-.*", "", group2))
  }
  #Get significant sirnas
  signif_sirnas <- comparison_df[comparison_df$p.adj <= 0.05,"group1"]
  #Get untested sirnas
  untested_sirnas <- c(comparison_df[is.na(comparison_df$p),"group1"], "CON", "CON1", "CON2")
  #Remove control siRNAs from results
  plot_df <- plot_df %>% dplyr::filter(siRNA != "CON")
  #Append on significance info to the plot_df
  plot_df <- plot_df %>% dplyr::mutate(signif = ifelse(siRNA %in% untested_sirnas, "Cat 3", ifelse(siRNA %in% signif_sirnas, "Cat 2", "Cat 1")))
  return(plot_df)
}

#Call the function on all the dataframes
hfob_alp_normal <- append_signif(hfob_alp_normal, hfob_alp_comparison_df)
hfob_alp_count_normal <- append_signif(hfob_alp_count_normal, hfob_alp_count_comparison_df)
hfob_alp_intensity_dapi_normal <- append_signif(hfob_alp_intensity_dapi_normal, hfob_alp_intensity_dapi_comparison_df)
hmsc_alp_normal <- append_signif(hmsc_alp_normal, hmsc_alp_comparison_df)
hmsc_alp_count_normal <- append_signif(hmsc_alp_count_normal, hmsc_alp_count_comparison_df)
hmsc_alp_intensity_dapi_normal <- append_signif(hmsc_alp_intensity_dapi_normal, hmsc_alp_intensity_dapi_comparison_df)
hmsc_ars_normal <- append_signif(hmsc_ars_normal, hmsc_ars_comparison_df)
hmsc_ars_count_normal <- append_signif(hmsc_ars_count_normal, hmsc_ars_count_comparison_df)
hmsc_ars_intensity_dapi_normal <- append_signif(hmsc_ars_intensity_dapi_normal, hmsc_ars_intensity_dapi_comparison_df)
hmsc_adipo_normal <- append_signif(hmsc_adipo_normal, hmsc_adipo_comparison_df)
hmsc_adipo_count_normal <- append_signif(hmsc_adipo_count_normal, hmsc_adipo_count_comparison_df)
hmsc_adipo_intensity_dapi_normal <- append_signif(hmsc_adipo_intensity_dapi_normal, hmsc_adipo_intensity_dapi_comparison_df)

# 8) Make Box Plots of Results ====

#Combine dataframes across the 4 experiments by type
raw_intensity_df <- rbind.data.frame(
  hfob_alp_normal %>% dplyr::select(siRNA, Replicate, signif, Value) %>% dplyr::mutate(subplot = "hFOB ALP"),
  hmsc_alp_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ALP"),
  hmsc_ars_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA,Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ARS"),
  hmsc_adipo_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC Adipo")
)
count_df <- rbind.data.frame(
  hfob_alp_count_normal %>% dplyr::select(siRNA, Replicate, signif, Value) %>% dplyr::mutate(subplot = "hFOB ALP"),
  hmsc_alp_count_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ALP"),
  hmsc_ars_count_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ARS"),
  hmsc_adipo_count_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC Adipo")
)
intensity_dapi_df <- rbind.data.frame(
  hfob_alp_intensity_dapi_normal %>% dplyr::select(siRNA, Replicate, signif, Value) %>% dplyr::mutate(subplot = "hFOB ALP"),
  hmsc_alp_intensity_dapi_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ALP"),
  hmsc_ars_intensity_dapi_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC ARS"),
  hmsc_adipo_intensity_dapi_normal %>% dplyr::filter(Treatment != "CON") %>% dplyr::select(siRNA, Donor, signif, Value) %>% dplyr::rename(Replicate = Donor) %>% dplyr::mutate(subplot = "hMSC Adipo")
)

#Reset factor levels
reset_factor_levels <- function(combined_df){
  #Reset siRNA factor levels
  combined_df$siRNA <- as.factor(combined_df$siRNA)
  #Set signif factor levels
  combined_df$signif <- factor(combined_df$signif, levels = c("Cat 1", "Cat 2", "Cat 3"))
  #Set subplot levels
  combined_df$subplot <- factor(combined_df$subplot, levels = c("hFOB ALP", "hMSC ALP", "hMSC ARS", "hMSC Adipo"))
  #Reset Plate Levels
  combined_df$Plate <- factor(paste0("Plate ", combined_df$Plate), levels = c("Plate 1", "Plate 2"))
  #Return the dataframe
  return(combined_df)
}
raw_intensity_df <- reset_factor_levels(raw_intensity_df)
count_df <- reset_factor_levels(count_df)
intensity_dapi_df <- reset_factor_levels(intensity_dapi_df)

#Make raw intensity plot
intensity_plot <- ggplot(raw_intensity_df, aes(x = siRNA, y = Value, colour = signif)) +
  geom_boxplot(linewidth=1.03, width=0.5) + 
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
  facet_grid(subplot ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05      ", "Not Tested     "),
                     guide = guide_legend(override.aes = list(color = "white"))) + 
  ylab("Fold Change Stain Intensity") + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                 axis.title.x = element_blank(), legend.position = "bottom",
                 axis.title.y = element_text(size = 16), 
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
                 legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 12, colour = "white"), 
                 axis.text.y = element_text(size = 12),
                 strip.text = element_text(size = 12))

#Make DAPI cell-count plot
count_plot <- ggplot(count_df, aes(x = siRNA, y = Value, colour = signif)) +
  geom_boxplot(linewidth=1.03, width=0.5) + 
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
  facet_grid(subplot ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05      ", "Not Tested     ")) + 
  ylab("Fold Change DAPI Cell Count") + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                 axis.title.x = element_blank(), legend.position = "bottom",
                 axis.title.y = element_text(size = 16),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
                 legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 12), 
                 axis.text.y = element_text(size = 12),
                 strip.text = element_text(size = 12)) +
  ggplot2::guides(color = ggplot2::guide_legend(
    keywidth = 0.0,
    keyheight = 0.3, 
    default.unit = "inch",
    override.aes = list(size = 5)))

#Make Plot of intensity normalized by cell-count
intensity_dapi_plot <- ggplot(intensity_dapi_df, aes(x = siRNA, y = Value, colour = signif)) +
  geom_boxplot(linewidth=1.03, width=0.5) + 
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", linewidth = 1.4, alpha = 0.6) +
  facet_grid(subplot ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c("gray", "dodgerblue3", "#48494B"), labels=c("Not Significant     ", "Adj. P-Value < 0.05      ", "Not Tested     "),
                     guide = guide_legend(override.aes = list(color = "white"))) + 
  ylab("Fold Change Intensity / DAPI Cell Count") + 
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                 axis.title.x = element_blank(), legend.position = "bottom",
                 axis.title.y = element_text(size = 16),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"), panel.grid.major.x = element_blank(),
                 legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
                 legend.title = ggplot2::element_blank(),
                 legend.text = element_text(size = 12, colour = "white"), 
                 axis.text.y = element_text(size = 12),
                 strip.text = element_text(size = 12))

#Make combined plot
tiff(paste0(out_dir, "main_assay_figure.main.box.tif"), width = 18000, height=6000, res=1000)
ggarrange(intensity_plot, count_plot, intensity_dapi_plot, ncol = 3, nrow = 1, common.legend = FALSE, align = "v", heights = c(4.25, 4.25), 
          labels = c("A", "B", "C"), font.label = list(face = "plain", size = 20))
dev.off()

# 9) Make Supplemental Tables of Results ====

### Make a Single Table for all hFOB ALP Results ###
rbind.data.frame(
  #Intensity Values
  raw_intensity_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hFOB ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    inner_join(hfob_alp_comparison_df, join_by(siRNA ==group1)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Stain Intensity"),
  #DAPI Cell Counts
  count_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hFOB ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    inner_join(hfob_alp_count_comparison_df, join_by(siRNA ==group1)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Cell Count"),
  #DAPI Cell Counts
  intensity_dapi_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hFOB ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    inner_join(hfob_alp_intensity_dapi_comparison_df, join_by(siRNA ==group1)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Intensity / Cell Count")
) %>%
write.table(file = paste0(inp_dir, "hfob_alp.plate-normalized.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### Make a Single Table for all hMSC ALP Results ###
rbind.data.frame(
  #Intensity Values
  raw_intensity_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_alp_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Stain Intensity"),
  #DAPI Cell Counts
  count_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_alp_count_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Cell Count"),
  #DAPI Cell Counts
  intensity_dapi_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ALP") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_alp_intensity_dapi_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Intensity / Cell Count")
) %>%
  write.table(file = paste0(inp_dir, "hmsc_alp.plate-normalized.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### Make a Single Table for all hMSC ARS Results ###
rbind.data.frame(
  #Intensity Values
  raw_intensity_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ARS") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_ars_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Stain Intensity"),
  #DAPI Cell Counts
  count_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ARS") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_ars_count_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Cell Count"),
  #DAPI Cell Counts
  intensity_dapi_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC ARS") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "BMP2"), temp2 = "CON-BMP2") %>%
    inner_join(hmsc_ars_intensity_dapi_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Intensity / Cell Count")
) %>%
  write.table(file = paste0(inp_dir, "hmsc_ars.plate-normalized.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

### Make a Single Table for all hMSC Adipogenesis Results ###
rbind.data.frame(
  #Intensity Values
  raw_intensity_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC Adipo") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "Adipo"), temp2 = "CON-Adipo") %>%
    inner_join(hmsc_adipo_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Stain Intensity"),
  #DAPI Cell Counts
  count_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC Adipo") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "Adipo"), temp2 = "CON-Adipo") %>%
    inner_join(hmsc_adipo_count_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Cell Count"),
  #DAPI Cell Counts
  intensity_dapi_df %>% dplyr::select(!Plate) %>%
    dplyr::filter(subplot == "hMSC Adipo") %>% 
    dplyr::group_by(siRNA) %>%
    dplyr::summarize(mean_fold_change=mean(Value)) %>%
    dplyr::mutate(temp = paste0(siRNA, "-", "Adipo"), temp2 = "CON-Adipo") %>%
    inner_join(hmsc_adipo_intensity_dapi_comparison_df, join_by(temp ==group1, temp2==group2)) %>% 
    dplyr::select(siRNA, reps, mean_fold_change, p, p.adj) %>%
    dplyr::mutate(p = ifelse(is.na(p), "Not Tested", p), p.adj = ifelse(is.na(p.adj), "Not Tested", p.adj)) %>%
    dplyr::mutate(type = "Intensity / Cell Count")
) %>%
  write.table(file = paste0(inp_dir, "hmsc_adipo.plate-normalized.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
