################################################################################

#analyze_differential_expression.R

#The purpose of this script is to analyze the test results for each of the four
#testing approaches we tried and determine DE genes.The script is set to run 
#with the following modules:

#  R/4.2.3

################################################################################

# 0) Call libraries and set directories, file locations, and universal variables ====

#Call libraries
library(stringr)
library(pbapply)
library(ggplot2)

#Set directories and file locations
inp_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/crispri_screen/HFOB_screen_bone/differential_expression/"

#Set test and cell types
test_types <- c("in_out", "ntc")
cell_types <- c("single_sgrna", "all")

# 1) Create analysis functions ====

#Define a function to calculate percentiles
calc_percentile_on_df <- function(inp_df){
  inp_df <- inp_df[order(inp_df$p_val, decreasing = TRUE),]
  rownames(inp_df) <- NULL
  inp_df$p_val_percentile <- as.integer(rownames(inp_df))/nrow(inp_df)
  return(inp_df)
}

#Create a function for calculating expected p-values (Both dataframes must be sorted per above function)
calc_expected_p_val <- function(inp_df_row, null_df){
  inp_p_val_percentile <- as.numeric(inp_df_row["p_val_percentile"])
  return(min(null_df$p_val[null_df$p_val_percentile <= inp_p_val_percentile]))
}

#Create a function for calculating empirical p-values (Both dataframes must be sorted per above function)
calc_empirical_p_val <- function(test_df_row, ntc_df){
  test_p_val <- as.numeric(test_df_row["p_val"])
  return((sum(ntc_df$p_val <= test_p_val) + 1)/(nrow(ntc_df) + 1))
}

#Define main analysis function
analyze_de_results <- function(inp_dir, test_type, cell_type, p_val_thresh=0.10){
  
  #Create output dir
  inp_dir <- ifelse(substr(inp_dir, nchar(inp_dir), nchar(inp_dir)) == "/", inp_dir, paste0(inp_dir, "/"))
  out_dir <- paste0(inp_dir, test_type, "_", cell_type, "/")
  
  #Read in RDS objects
  null_results <- readRDS(paste0(inp_dir, "null_", test_type, "_", cell_type, "_df.rds"))
  pos_results <- readRDS(paste0(inp_dir, "pos_control_", test_type, "_", cell_type, "_df.rds"))
  ntc_results <- readRDS(paste0(inp_dir, "non_targeting_", test_type, "_", cell_type, "_df.rds"))
  test_results <- readRDS(paste0(inp_dir, "targeting_", test_type, "_", cell_type, "_df.rds"))
  
  #Recast matrix test results to dataframes
  null_df <- as.data.frame(null_results)
  pos_df <- as.data.frame(pos_results)
  ntc_df <- as.data.frame(ntc_results)
  test_df <- as.data.frame(test_results)
  
  #Recast p-value columns as numerics
  null_df$p_val <- as.numeric(null_df$p_val)
  pos_df$p_val <- as.numeric(pos_df$p_val)
  ntc_df$p_val <- as.numeric(ntc_df$p_val)
  test_df$p_val <- as.numeric(test_df$p_val)
  #Calculate p_value quantiles
  ntc_df <- calc_percentile_on_df(ntc_df)
  null_df <- calc_percentile_on_df(null_df)
  test_df <- calc_percentile_on_df(test_df)
  #Calculate expected p-values
  ntc_df$expected_p_val <- pbapply(ntc_df, FUN = calc_expected_p_val, MARGIN = 1, null_df=null_df)
  test_df$expected_p_val <- pbapply(test_df, FUN = calc_expected_p_val, MARGIN = 1, null_df=null_df)
  
  #Make plot df
  plot_df <- rbind.data.frame(
    cbind.data.frame(category=rep("test guide", nrow(test_df)), observed_p_val=test_df$p_val, expected_p_val=test_df$expected_p_val),
    cbind.data.frame(category=rep("ntc guide", nrow(ntc_df)), observed_p_val=ntc_df$p_val, expected_p_val=ntc_df$expected_p_val)
    )
  plot_df$neg_log_observed <- -log10(plot_df$observed_p_val)
  plot_df$neg_log_expected <- -log10(plot_df$expected_p_val)
  
  #Make plot
  temp_plot <- ggplot(data = plot_df, aes(x=neg_log_expected, y=neg_log_observed, color=category)) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) + 
    xlab("-log10(expected P-value)") + xlim(c(0, max(plot_df$neg_log_expected))) + 
    ylab("-log10(observed P-value)") + ylim(c(0, max(plot_df$neg_log_observed))) + 
    theme_minimal() +
    scale_color_manual(values = c("#D4D4D4", "#F87850")) + 
    theme(legend.title=element_blank(), axis.title=element_text(size=16), legend.position = "none", 
          legend.text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), legend.box="vertical") + 
  jpeg(paste0(out_dir, "qq.scatter.jpeg"), width = 6, height = 11.5, units = 'in', res = 500)
  print(temp_plot)
  dev.off()
  
  #Create empirical FDR for test_df and positive control df
  test_df$empirical_p_val <- pbapply(test_df, FUN = calc_empirical_p_val, MARGIN = 1, ntc_df=ntc_df)
  pos_df$empirical_p_val <- pbapply(pos_df, FUN = calc_empirical_p_val, MARGIN = 1, ntc_df=ntc_df)
  #Calculated adjusted empirical p-value
  test_df$empirical_p_val_adj <- p.adjust(test_df$empirical_p_val, method = "BH")
  pos_df$empirical_p_val_adj <- p.adjust(pos_df$empirical_p_val, method = "BH")
  
  #Create test export dataframe and export
  export_df <- test_df[,c("sgrna", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "empirical_p_val", "empirical_p_val_adj")]
  export_df <- export_df[order(export_df$p_val),]
  export_df <- export_df[which(export_df$empirical_p_val_adj <= p_val_thresh),]
  write.table(export_df, file = paste0(out_dir,test_type,"_",cell_type,"_significant_results.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
  #Do the same for pos dataframe
  export_df <- pos_df[,c("sgrna", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "empirical_p_val", "empirical_p_val_adj")]
  export_df <- export_df[order(export_df$p_val),]
  export_df <- export_df[which(export_df$empirical_p_val_adj <= p_val_thresh),]
  write.table(export_df, file = paste0(out_dir,test_type,"_",cell_type,"_significant_pos_control_results.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

} 

# 2) Create KS-Test Function for Testing True non-Targeting Ness of NTCs ====

#Create execution function for test
execute_KS <- function(ntc_sgrna, ntc_df, col_name="avg_log2FC", test_alternative = "two.sided"){
  target_vals <- as.numeric(ntc_df[which(ntc_df$sgrna == ntc_sgrna), col_name])
  compare_vals <- as.numeric(ntc_df[which(ntc_df$sgrna != ntc_sgrna), col_name])
  temp <- suppressWarnings(ks.test(target_vals, compare_vals, alternative = test_alternative))
  return(temp$p.value)
}

#Define main wrapper function
check_KS <- function(inp_dir, test_type, cell_type, col_name="p_val") {
  
  #Create output dir
  inp_dir <- ifelse(substr(inp_dir, nchar(inp_dir), nchar(inp_dir)) == "/", inp_dir, paste0(inp_dir, "/"))
  out_dir <- paste0(inp_dir, test_type, "_", cell_type, "/")
  
  #Read in ntc RDS object
  ntc_results <- readRDS(paste0(inp_dir, "non_targeting_", test_type, "_", cell_type, "_df.rds"))
  #Recast matrix test results to dataframes
  ntc_df <- as.data.frame(ntc_results)
  #Recast p-value columns as numerics
  ntc_df$p_val <- as.numeric(ntc_df$p_val)
  
  #Get list of unique ntc sgrnas
  ntc_sgrnas <- unique(ntc_df$sgrna)
  #Execute tests
  temp <- vapply(ntc_sgrnas, FUN = execute_KS, FUN.VALUE = 0, ntc_df=ntc_df, col_name=col_name, test_alternative="two.sided")
  
  #Return list
  return(temp)
  
}

# 3) Create a Wilcoxon Function for Testing True non-Targeting Ness of NTCs ====


#Create execution function for test
execute_wilcox <- function(ntc_sgrna, ntc_df, col_name="p_val", test_alternative = "less"){
  target_vals <- as.numeric(ntc_df[which(ntc_df$sgrna == ntc_sgrna), col_name])
  compare_vals <- as.numeric(ntc_df[which(ntc_df$sgrna != ntc_sgrna), col_name])
  temp <- suppressWarnings(wilcox.test(target_vals, compare_vals, alternative = test_alternative))
  return(temp$p.value)
}

#Define main wrapper function
check_wilcox <- function(inp_dir, test_type, cell_type, col_name="p_val") {
  
  #Create output dir
  inp_dir <- ifelse(substr(inp_dir, nchar(inp_dir), nchar(inp_dir)) == "/", inp_dir, paste0(inp_dir, "/"))
  out_dir <- paste0(inp_dir, test_type, "_", cell_type, "/")
  
  #Read in ntc RDS object
  ntc_results <- readRDS(paste0(inp_dir, "non_targeting_", test_type, "_", cell_type, "_df.rds"))
  #Recast matrix test results to dataframes
  ntc_df <- as.data.frame(ntc_results)
  #Recast p-value columns as numerics
  ntc_df$p_val <- as.numeric(ntc_df$p_val)
  
  #Get list of unique ntc sgrnas
  ntc_sgrnas <- unique(ntc_df$sgrna)
  #Execute tests
  temp <- vapply(ntc_sgrnas, FUN = execute_wilcox, FUN.VALUE = 0, ntc_df=ntc_df, col_name=col_name, test_alternative="less")
  
  #Return list
  return(temp)
  
}


# 4) Call Analysis Functions ====

#Loop over cell types and test types and call function
for (test_type in test_types) {
  for (cell_type in cell_types) {
    analyze_de_results(inp_dir, test_type, cell_type)
  }
}
