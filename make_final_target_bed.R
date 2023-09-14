# 1) Load libraries ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(eulerr)
library(stringi)
library(stringr)

# 2) Read in data and Process ====

#Set directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/varianttogene/bone_cells/"

#Read in final targets file
final_targets_raw <- read.table(paste0(inp_dir, "final_signal_targets.txt"), header = TRUE, sep = "\t")

#Make list of proxy snps
signal_proxies <- final_targets_raw$proxies
signal_proxies <- strsplit(signal_proxies, ", ")
#Convert to a vector format
proxies_vector <- vector()
for (i in 1:length(signal_proxies)) {
  proxies_vector <- append(proxies_vector, signal_proxies[[i]])
}
#Remove parentheses
proxies_vector <- str_replace(proxies_vector, pattern = "[)]", replacement = "")
proxies_vector <- str_replace(proxies_vector, pattern = "[(]", replacement = ":")
#Split string
proxies_df <- as.data.frame(str_split(proxies_vector, ":", simplify = TRUE))
#Reorganize columns
proxies_df <- proxies_df[,c(2,3,1)]
#Convert bp to numeric
colnames(proxies_df) <- c("CHR", "BP", "SNP")
proxies_df$BP <- as.numeric(proxies_df$BP)

#Strip the rsid from final_targets_raw and append to a separate column
final_targets_raw <- cbind.data.frame(final_targets_raw, 
                                      just_proxies=stri_replace_all(final_targets_raw$proxies, replacement = "", regex = "[(]chr\\d+:\\d+[)]"))

final_targets_raw$just_proxies <- stri_replace_all(final_targets_raw$just_proxies, replacement = ",", regex = ", ")

#Append mid positions to final_targets_raw
append_mid <- function(targets, proxies_df){
  targets_vec <- str_split(targets, ",")[[1]]
  temp_df <- proxies_df[which(proxies_df$SNP %in% targets_vec),]
  min_pos <- min(temp_df$BP)
  max_pos <- max(temp_df$BP)
  mid_pos <- round((max_pos - min_pos)/2, 0) + min_pos
  dist_to_far <- max(temp_df$BP - mid_pos)
  return(c(mid_pos, dist_to_far))
}
final_targets_raw <- cbind.data.frame(final_targets_raw, 
                                      t(vapply(final_targets_raw$just_proxies, append_mid, proxies_df=proxies_df, FUN.VALUE = numeric(2))))

#Make output df
output_df <- cbind.data.frame(CHR=paste0("chr",final_targets_raw$CHR), BP_START=final_targets_raw["1"], BP_END=final_targets_raw["1"] + 1,
                              SNP=final_targets_raw$just_proxies, MAX_DIST_TO_MID=final_targets_raw["2"])

colnames(output_df) <- c("CHR", "BP_START", "BP_END", "SNP", "MAX_DIST_TO_MID")

#Combine the two paired rows
#rs17130569,rs7554551
#rs73232633,rs73402679
combine_rows <- function(row1, row2){
  chromo=row1$CHR
  BP_START=round(as.numeric(row1$BP_START)/2+as.numeric(row2$BP_START)/2,0)
  BP_END=BP_START+1
  SNP=paste(row1$SNP,row2$SNP, sep = ",")
  MAX_DIST_TO_MID=max(row1$BP_START-BP_START, row2$BP_START-BP_START)
  return(c(chromo, BP_START, BP_END, SNP, MAX_DIST_TO_MID))
}
merge_rows <- function(SNP1, SNP2, output_df){
  row1 <- output_df[which(output_df$SNP == SNP1),]
  row2 <- output_df[which(output_df$SNP == SNP2),]
  temp_df <- output_df[which(!(output_df$SNP %in% c(SNP1, SNP2))),]
  temp_df <- rbind.data.frame(temp_df, combine_rows(row1, row2))
  temp_df$BP_START <- as.numeric(temp_df$BP_START)
  temp_df$BP_END <- as.numeric(temp_df$BP_END)
  temp_df$MAX_DIST_TO_MID <- as.numeric(temp_df$MAX_DIST_TO_MID)
  return(temp_df)
}
output_df <- merge_rows("rs7554551", "rs17130569", output_df)
output_df <- merge_rows("rs73232633", "rs73402679", output_df)

#Sort the bed file
output_df <- output_df[order(output_df$CHR, output_df$BP_START),]
#Write to file
write.table(output_df, file = paste0(inp_dir, "osteoblast_target_regions.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




