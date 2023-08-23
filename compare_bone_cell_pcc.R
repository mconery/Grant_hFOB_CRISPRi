# 1) Load libraries ====

#Load libraries
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(eulerr)

# 2) Read in data and Minimally Process ====

#Set directory
inp_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/varianttogene/bone_cells/"
sub_dirs <- c("hMSC_Osteo", "hMSC_Control", "hFOBsDiff", "hFOBsundiff")

#Read in files
hMSC_Osteo_promoter_raw <- read.table(paste0(inp_dir, sub_dirs[1],"/ld_proxies.", sub_dirs[1], ".snpOpenPromoter.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hFOBsDiff_promoter_raw <- read.table(paste0(inp_dir, sub_dirs[3],"/ld_proxies.", sub_dirs[3], ".snpOpenPromoter.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hFOBsundiff_promoter_raw <- read.table(paste0(inp_dir, sub_dirs[4],"/ld_proxies.", sub_dirs[4], ".snpOpenPromoter.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hMSC_Osteo_interaction_raw <- read.table(paste0(inp_dir, sub_dirs[1],"/ld_proxies.", sub_dirs[1], ".snpOE_interaction.both.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hFOBsDiff_interaction_raw <- read.table(paste0(inp_dir, sub_dirs[3],"/ld_proxies.", sub_dirs[3], ".snpOE_interaction.both.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
hFOBsundiff_interaction_raw <- read.table(paste0(inp_dir, sub_dirs[4],"/ld_proxies.", sub_dirs[4], ".snpOE_interaction.both.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#Make a function to remove TCONS and then remove TCONS
remove_TCONS <- function(implication_raw){
  temp <- implication_raw[which(substr(implication_raw$IMPLICATED_GENE, 1, 5) != "TCONS"),]
  return(temp)
}
hMSC_Osteo_promoter_scrubbed <- remove_TCONS(hMSC_Osteo_promoter_raw)
hFOBsDiff_promoter_scrubbed <- remove_TCONS(hFOBsDiff_promoter_raw)
hFOBsundiff_promoter_scrubbed <- remove_TCONS(hFOBsundiff_promoter_raw)
hMSC_Osteo_interaction_scrubbed <- remove_TCONS(hMSC_Osteo_interaction_raw)
hFOBsDiff_interaction_scrubbed <- remove_TCONS(hFOBsDiff_interaction_raw)
hFOBsundiff_interaction_scrubbed <- remove_TCONS(hFOBsundiff_interaction_raw)

#Make a table of sentinel/gene relationships
#First get vectors from each table
hMSC_Osteo_promoter_relations <- unique(paste(hMSC_Osteo_promoter_scrubbed$SENTINEL, hMSC_Osteo_promoter_scrubbed$IMPLICATED_GENE, sep = "|"))
hFOBsDiff_promoter_relations <- unique(paste(hFOBsDiff_promoter_scrubbed$SENTINEL, hFOBsDiff_promoter_scrubbed$IMPLICATED_GENE, sep = "|"))
hFOBsundiff_promoter_relations <- unique(paste(hFOBsundiff_promoter_scrubbed$SENTINEL, hFOBsundiff_promoter_scrubbed$IMPLICATED_GENE, sep = "|"))
hMSC_Osteo_interaction_relations <- unique(paste(hMSC_Osteo_interaction_scrubbed$SENTINEL, hMSC_Osteo_interaction_scrubbed$IMPLICATED_GENE, sep = "|"))
hFOBsDiff_interaction_relations <- unique(paste(hFOBsDiff_interaction_scrubbed$SENTINEL, hFOBsDiff_interaction_scrubbed$IMPLICATED_GENE, sep = "|"))
hFOBsundiff_interaction_relations <- unique(paste(hFOBsundiff_interaction_scrubbed$SENTINEL, hFOBsundiff_interaction_scrubbed$IMPLICATED_GENE, sep = "|"))

#Make a dataframe of results
relation_df <- rbind.data.frame(cbind.data.frame(type=rep("hMSC_Osteo.promoter", length(hMSC_Osteo_promoter_relations)),relation=hMSC_Osteo_promoter_relations),
                                cbind.data.frame(type=rep("hFOBsDiff.promoter", length(hFOBsDiff_promoter_relations)),relation=hFOBsDiff_promoter_relations),
                                cbind.data.frame(type=rep("hFOBsundiff.promoter", length(hFOBsundiff_promoter_relations)),relation=hFOBsundiff_promoter_relations),
                                cbind.data.frame(type=rep("hMSC_Osteo.interaction", length(hMSC_Osteo_interaction_relations)),relation=hMSC_Osteo_interaction_relations),
                                cbind.data.frame(type=rep("hFOBsDiff.interaction", length(hFOBsDiff_interaction_relations)),relation=hFOBsDiff_interaction_relations),
                                cbind.data.frame(type=rep("hFOBsundiff.interaction", length(hFOBsundiff_interaction_relations)),relation=hFOBsundiff_interaction_relations))

#Make usable table
relation_table <- data.frame(matrix(data=0, nrow = length(unique(relation_df$relation)), ncol = length(unique(relation_df$type))))
rownames(relation_table) <- unique(relation_df$relation)
colnames(relation_table) <- unique(relation_df$type)
#Fill the lazy way
for (i in 1:nrow(relation_table)) {
  for (j in 1:ncol(relation_table)) {
    relation_table[i,j] <- sum(relation_df$type == colnames(relation_table)[j] & relation_df$relation == rownames(relation_table)[i])
  }
}
#Create Either implications
relation_table <- cbind.data.frame(relation_table, 
                                   hMSC_Osteo.either=ifelse(relation_table$hMSC_Osteo.promoter==1 | relation_table$hMSC_Osteo.interaction==1, 1, 0),
                                   hFOBsDiff.either=ifelse(relation_table$hFOBsDiff.promoter==1 | relation_table$hFOBsDiff.interaction==1, 1, 0),
                                   hFOBsundiff.either=ifelse(relation_table$hFOBsundiff.promoter==1 | relation_table$hFOBsundiff.interaction==1, 1, 0))

#Make a master implication file with proxies and genes
imp_table <- cbind.data.frame(rbind.data.frame(hMSC_Osteo_promoter_scrubbed[,1:9], hMSC_Osteo_interaction_scrubbed[,1:9], 
                                  hFOBsDiff_promoter_scrubbed[,1:9], hFOBsDiff_interaction_scrubbed[,1:9], 
                                  hFOBsundiff_promoter_scrubbed[,1:9], hFOBsundiff_interaction_scrubbed[,1:9]),
                 cell_type=c(rep("hMSC_Osteo", nrow(hMSC_Osteo_promoter_scrubbed) + nrow(hMSC_Osteo_interaction_scrubbed)),
                             rep("hFOBsDiff", nrow(hFOBsDiff_promoter_scrubbed) + nrow(hFOBsDiff_interaction_scrubbed)),
                             rep("hFOBsundiff", nrow(hFOBsundiff_promoter_scrubbed) + nrow(hFOBsundiff_interaction_scrubbed))),
                 imp_type=c(rep("promoter", nrow(hMSC_Osteo_promoter_scrubbed)),
                            rep("interaction", nrow(hMSC_Osteo_interaction_scrubbed)),
                            rep("promoter", nrow(hFOBsDiff_promoter_scrubbed)),
                            rep("interaction", nrow(hFOBsDiff_interaction_scrubbed)),
                            rep("promoter", nrow(hFOBsundiff_promoter_scrubbed)),
                            rep("interaction", nrow(hFOBsundiff_interaction_scrubbed))))

################################################## Make either Venn ###############################################################
#Make Venn Diagram
VennDiag <- euler((relation_table[,c(7,8,9)]))
jpeg(paste(inp_dir, "implication_combos.either.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

################################################## Make interaction Venn #################################################################
#Make Venn Diagram
VennDiag <- euler((relation_table[,c(4,5,6)]))
jpeg(paste(inp_dir, "implication_combos.interaction.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

################################################## Make promoter Venn ####################################################################
#Make Venn Diagram
VennDiag <- euler((relation_table[,c(1,2,3)]))
jpeg(paste(inp_dir, "implication_combos.promoter.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

##########################################################################################################################################

#Read in expression data
tpm_dir <- "C:/Users/mitch/Documents/UPenn/Voight_Lab/Trait_Colocalization/data/grant_osteoblasts/"
hMSC_tpm_raw <- read.table(paste0(tpm_dir, "osteoblast_tpm.txt"), header=TRUE)
hFOB_tpm_raw <- read.table(paste0(tpm_dir, "hFOB_tpm.txt"), header=TRUE)
#Tack on gene/sentinel variant info
temp <- as.data.frame(strsplit(rownames(relation_table), "[|]"))
relation_table <- cbind.data.frame(relation_table, t(as.data.frame(strsplit(rownames(relation_table), "[|]"))))
colnames(relation_table)[length(colnames(relation_table))-1] <- "RSID"
colnames(relation_table)[length(colnames(relation_table))] <- "IMPLICATED_GENE"

#Clean-up junk
remove(hMSC_Osteo_promoter_raw, hFOBsDiff_promoter_raw, hFOBsundiff_promoter_raw, VennDiag)
remove(hMSC_Osteo_interaction_raw, hFOBsDiff_interaction_raw, hFOBsundiff_interaction_raw, temp)

#Make dataframe of average tpm of genes across samples
hMSC_Osteo_average_tpm <- aggregate(hMSC_tpm_raw$tpm, by = list(hMSC_tpm_raw$gene_name), FUN=mean)
colnames(hMSC_Osteo_average_tpm) <- c("gene_name", "tpm")
rownames(hMSC_Osteo_average_tpm) <- hMSC_Osteo_average_tpm$gene_name
hFOBsDiff_tpm_raw <- hFOB_tpm_raw[which(substr(hFOB_tpm_raw$sample, 1, 8) == "hFOBdiff"),]
hFOBsundiff_tpm_raw <- hFOB_tpm_raw[which(substr(hFOB_tpm_raw$sample, 1, 10) == "hFOBundiff"),]
hFOBsDiff_average_tpm <- aggregate(hFOBsDiff_tpm_raw$tpm, by = list(hFOBsDiff_tpm_raw$gene_name), FUN=mean)
colnames(hFOBsDiff_average_tpm) <- c("gene_name", "tpm")
rownames(hFOBsDiff_average_tpm) <- hFOBsDiff_average_tpm$gene_name
hFOBsundiff_average_tpm <- aggregate(hFOBsundiff_tpm_raw$tpm, by = list(hFOBsundiff_tpm_raw$gene_name), FUN=mean)
colnames(hFOBsundiff_average_tpm) <- c("gene_name", "tpm")
rownames(hFOBsundiff_average_tpm) <- hFOBsundiff_average_tpm$gene_name

#Make function to extract expression data
get_expr <- function(imp_gene, average_tpm_df){
  temp <- average_tpm_df[which(average_tpm_df$gene_name == imp_gene),]
  if(nrow(temp) > 0){
    return(mean(temp$tpm))
  } else {
    return(0)
  }
}
#Tack on expression info to the relation table
relation_table_expr <- cbind.data.frame(relation_table, 
                                        hMSC_Osteo_expr=vapply(relation_table$IMPLICATED_GENE, FUN = get_expr, average_tpm_df = hMSC_Osteo_average_tpm, FUN.VALUE =  1),
                                        hFOBsDiff_expr=vapply(relation_table$IMPLICATED_GENE, FUN = get_expr, average_tpm_df = hFOBsDiff_average_tpm, FUN.VALUE =  1),
                                        hFOBsundiff_expr=vapply(relation_table$IMPLICATED_GENE, FUN = get_expr, average_tpm_df = hFOBsundiff_average_tpm, FUN.VALUE =  1))
#Make implication columns accounting for expression
tpm_cutoff  <- 1 #tpm
relation_table_expr <- cbind.data.frame(relation_table_expr,
                                        hMSC_Osteo.promoter.expr=ifelse(relation_table_expr$hMSC_Osteo.promoter == 1 & relation_table_expr$hMSC_Osteo_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsDiff.promoter.expr=ifelse(relation_table_expr$hFOBsDiff.promoter == 1 & relation_table_expr$hFOBsDiff_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsundiff.promoter.expr=ifelse(relation_table_expr$hFOBsundiff.promoter == 1 & relation_table_expr$hFOBsundiff_expr > tpm_cutoff, TRUE, FALSE),
                                        hMSC_Osteo.interaction.expr=ifelse(relation_table_expr$hMSC_Osteo.interaction == 1 & relation_table_expr$hMSC_Osteo_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsDiff.interaction.expr=ifelse(relation_table_expr$hFOBsDiff.interaction == 1 & relation_table_expr$hFOBsDiff_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsundiff.interaction.expr=ifelse(relation_table_expr$hFOBsundiff.interaction == 1 & relation_table_expr$hFOBsundiff_expr > tpm_cutoff, TRUE, FALSE),
                                        hMSC_Osteo.either.expr=ifelse(relation_table_expr$hMSC_Osteo.either == 1 & relation_table_expr$hMSC_Osteo_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsDiff.either.expr=ifelse(relation_table_expr$hFOBsDiff.either == 1 & relation_table_expr$hFOBsDiff_expr > tpm_cutoff, TRUE, FALSE),
                                        hFOBsundiff.either.expr=ifelse(relation_table_expr$hFOBsundiff.either == 1 & relation_table_expr$hFOBsundiff_expr > tpm_cutoff, TRUE, FALSE))

#Make function to extract expression data
check_imp_expr <- function(imp_row){
  imp_gene <- as.character(imp_row["IMPLICATED_GENE"])
  cell_type <- as.character(imp_row["cell_type"])
  temp_text <- paste0(cell_type, "_average_tpm")
  average_tpm_df <- eval(parse(text = temp_text))
  temp <- average_tpm_df[which(average_tpm_df$gene_name == imp_gene),]
  if(nrow(temp) > 0){
    return(mean(temp$tpm))
  } else {
    return(0)
  }
}
#Tack on expression data to the imp_table
imp_table <- cbind.data.frame(imp_table, 
                              expressed=ifelse(apply(imp_table, MARGIN = 1, FUN = check_imp_expr) > tpm_cutoff, TRUE, FALSE))



################################################## Make either expression Venn ###############################################################
#Make Venn Diagram
VennDiag <- euler((relation_table_expr[,c(21,22,23)]))
jpeg(paste(inp_dir, "implication_combos.either.expr.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

################################################## Make interaction expression Venn #################################################################
#Make Venn Diagram
VennDiag <- euler((relation_table_expr[,c(18,19,20)]))
jpeg(paste(inp_dir, "implication_combos.interaction.expr.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

################################################## Make promoter expression Venn ####################################################################
#Make Venn Diagram
VennDiag <- euler((relation_table_expr[,c(15,16,17)]))
jpeg(paste(inp_dir, "implication_combos.promoter.expr.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

##############################################################################################################################

#Tack on either sum columns
relation_table_expr <- cbind.data.frame(relation_table_expr, cells_imp=rowSums(relation_table_expr[,c(7,8,9)]),
                                        cells_imp.expr=rowSums(relation_table_expr[,c(21,22,23)]))
#Clean-up
remove(relation_df, relation_table, hMSC_Osteo_average_tpm, hFOBsDiff_average_tpm, hFOBsundiff_average_tpm)
remove(hFOB_tpm_raw, hFOBsDiff_tpm_raw, hFOBsundiff_tpm_raw, hMSC_tpm_raw)
remove(i, j, hFOBsDiff_interaction_scrubbed, hFOBsDiff_promoter_scrubbed, hFOBsundiff_interaction_scrubbed, hFOBsundiff_promoter_scrubbed, hMSC_Osteo_interaction_scrubbed, hMSC_Osteo_promoter_scrubbed)

#Read in Morris et al signal list and append information
morris_bmd_df <- read.table(paste0(inp_dir, "morris_bmd_COJO_signals.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#Create a function that takes an rsid and finds out how many/which genes are linked to it
find_imps <- function(rsid, relation_table_expr){
  temp <- relation_table_expr[which(relation_table_expr$RSID == rsid),]
  count_any = sum(ifelse(rowSums(temp[,c("hMSC_Osteo.either.expr", "hFOBsDiff.either.expr", "hFOBsundiff.either.expr")])>0,1,0))
  count_hMSC_Osteo = sum(temp$hMSC_Osteo.either.expr)
  count_hFOBsDiff = sum(temp$hFOBsDiff.either.expr)
  count_hFOBsundiff = sum(temp$hFOBsundiff.either.expr) 
  genes_hMSC_Osteo = ifelse(count_hMSC_Osteo > 0, paste(temp[which(temp$hMSC_Osteo.either.expr == TRUE),"IMPLICATED_GENE"], collapse = ","), "")
  genes_hFOBsDiff = ifelse(count_hFOBsDiff > 0, paste(temp[which(temp$hFOBsDiff.either.expr == TRUE),"IMPLICATED_GENE"], collapse = ","), "")
  genes_hFOBsundiff = ifelse(count_hFOBsundiff > 0, paste(temp[which(temp$hFOBsundiff.either.expr == TRUE),"IMPLICATED_GENE"], collapse = ","), "")
  genes_any = ifelse(count_any > 0, paste(temp$IMPLICATED_GENE[ifelse(rowSums(temp[,c("hMSC_Osteo.either.expr", "hFOBsDiff.either.expr", "hFOBsundiff.either.expr")])>0,TRUE,FALSE)], collapse = ","), "")
  return(cbind(count_hMSC_Osteo, count_hFOBsDiff, count_hFOBsundiff, count_any, genes_hMSC_Osteo, genes_hFOBsDiff, genes_hFOBsundiff, genes_any))
}
morris_bmd_imp <- cbind.data.frame(morris_bmd_df,
                                   t(sapply(morris_bmd_df$RSID, find_imps, relation_table_expr=relation_table_expr)))
colnames(morris_bmd_imp)[25:32] <- c("count_hMSC_Osteo", "count_hFOBsDiff", "count_hFOBsundiff", "count_any", 
                                     "genes_hMSC_Osteo", "genes_hFOBsDiff", "genes_hFOBsundiff", "genes_any")
morris_bmd_imp$count_hMSC_Osteo <- as.numeric(morris_bmd_imp$count_hMSC_Osteo)
morris_bmd_imp$count_hFOBsDiff <- as.numeric(morris_bmd_imp$count_hFOBsDiff)
morris_bmd_imp$count_hFOBsundiff <- as.numeric(morris_bmd_imp$count_hFOBsundiff)
morris_bmd_imp$count_any <- as.numeric(morris_bmd_imp$count_any)

#Make histograms of the number of signals implicated
make_histogram <- function(field_names, labels, ind_sig_df, output_file, max_val=FALSE){
  b <- min(ind_sig_df[,field_names], na.rm = TRUE)-1 # Set the minimum for the breakpoints
  e <- max(ind_sig_df[,field_names], na.rm = TRUE) # Set the maximum for the breakpoints
  if (max_val!=FALSE) {
    e <- max(e, max_val)
  }
  ax <- pretty(b:e, n = e-b) # Make a neat vector for the breakpoints
  colors <- rainbow(length(field_names), alpha = 1) #Get colors for the histogram
  jpeg(file=output_file)
  if (length(field_names) > 1) {
    hg <- hist(ind_sig_df[which(is.na(ind_sig_df[,field_names[1]]) == FALSE),field_names[1]], breaks = ax, plot = FALSE)
    plot(hg, col = colors[1], xlab = labels[2], main = labels[1], )
    for (i in 2:length(field_names)) {
      hg <- hist(ind_sig_df[which(is.na(ind_sig_df[,field_names[i]]) == FALSE),field_names[i]], breaks = ax, plot = FALSE)
      plot(hg, col = colors[i], add = TRUE)
    }
    legend("topright", legend = field_names, fill=colors)
  } else if (length(field_names) == 1) {
    hg <- hist(ind_sig_df[which(is.na(ind_sig_df[,field_names]) == FALSE),field_names], breaks = ax, plot = FALSE)
    plot(hg, col = colors[1], xlab = labels[2], main = labels[1])
  }
  dev.off()
}
make_histogram(c("count_any", "count_hMSC_Osteo", "count_hFOBsDiff", "count_hFOBsundiff"), 
               labels=c("Genes per Signal", "Genes"), ind_sig_df = morris_bmd_imp, output_file = paste0(inp_dir, "genes_per_signal.jpg"))

#Make venn diagrams
temp <- ifelse(morris_bmd_imp[,c(25,26,27)] > 0, 1, 0)
VennDiag <- euler((temp))
jpeg(paste(inp_dir, "genes_per_signal_venn.jpg", sep = ""), width = 5, height = 5, units = 'in', res = 200)
print(plot(VennDiag, quantities = TRUE, fill = c("red", "#03AC13", "#99EDC3"), alpha = c(0.9, 0.90, 0.2), 
           labels = list(labels=c("hMSC_Osteo", "hFOBsDiff", "hFOBsundiff"), cex = rep(1.1,3))))
dev.off()

#Count number of genes implicated
sum(morris_bmd_imp$count_any > 0) #305 signals with 1 or more expressed implicated genes
sum(morris_bmd_imp$count_any == 1) #157 implicate only 1 expressed gene and 148 implicate multiple expressed genes

#############################################################################################################################

#Clean-up
remove(VennDiag, morris_bmd_df)
remove(hFOBsDiff_interaction_relations, hFOBsDiff_promoter_relations, hMSC_Osteo_interaction_relations, hMSC_Osteo_promoter_relations, hFOBsundiff_interaction_relations, hFOBsundiff_promoter_relations)
remove(temp)

#Add the proxy info on to the relation_table_expr
get_proxies <- function(relation_row, imp_table){
  rsid <- as.character(relation_row["RSID"])
  imp_gene <- as.character(relation_row["IMPLICATED_GENE"])
  #Pull individual temp files for the all and interactions
  temp <- imp_table[which(imp_table$SENTINEL == rsid & imp_table$IMPLICATED_GENE == imp_gene),]
  temp_interaction <- imp_table[which(imp_table$SENTINEL == rsid & imp_table$IMPLICATED_GENE == imp_gene & imp_table$imp_type == "interaction"),]
  temp_expr <- imp_table[which(imp_table$SENTINEL == rsid & imp_table$IMPLICATED_GENE == imp_gene & imp_table$expressed == TRUE),]
  temp_interaction_expr <- imp_table[which(imp_table$SENTINEL == rsid & imp_table$IMPLICATED_GENE == imp_gene & imp_table$imp_type == "interaction" & imp_table$expressed == TRUE),]
  proxies <- unique(temp$PROXY)
  int_proxies <- unique(temp_interaction$PROXY)
  locations <- unique(temp$PROXY_POSITION)
  int_locations <- unique(temp_interaction$PROXY_POSITION)
  proxy_out <- paste0(proxies, "(", locations, ")", collapse = ", ")
  int_proxy_out <- paste0(int_proxies, "(", int_locations, ")", collapse = ", ")
  proxies_expr <- unique(temp_expr$PROXY)
  int_proxies_expr <- unique(temp_interaction_expr$PROXY)
  count_proxies <- length(proxies)
  count_int_proxies <- length(int_proxies)
  count_proxies_expr <- length(proxies_expr)
  count_int_proxies_expr <- length(int_proxies_expr)
  locations_expr <- unique(temp_expr$PROXY_POSITION)
  int_locations_expr <- unique(temp_interaction_expr$PROXY_POSITION)
  proxy_out_expr <- paste0(proxies_expr, "(", locations_expr, ")", collapse = ", ")
  int_proxy_out_expr <- paste0(int_proxies_expr, "(", int_locations_expr, ")", collapse = ", ")
  return(c(count_proxies, proxy_out, count_int_proxies, int_proxy_out, count_proxies_expr, proxy_out_expr, count_int_proxies_expr, int_proxy_out_expr))
}
#Tack on the proxy SNP info
relation_table_out <- cbind.data.frame(relation_table_expr, t(apply(relation_table_expr, MARGIN = 1, FUN = get_proxies, imp_table=imp_table)))
colnames(relation_table_out)[(ncol(relation_table_out) - 7):ncol(relation_table_out)] <- 
  c("count_all_proxies", "all_proxies", "count_interaction_proxies", "interaction_proxies", 
    "count_all_proxies_expr", "all_proxies_expr", "count_interaction_proxies_expr", "interaction_proxies_expr")
relation_table_out$count_all_proxies <- as.numeric(relation_table_out$count_all_proxies)
relation_table_out$count_interaction_proxies <- as.numeric(relation_table_out$count_interaction_proxies)
relation_table_out$count_all_proxies_expr <- as.numeric(relation_table_out$count_all_proxies_expr)
relation_table_out$count_interaction_proxies_expr <- as.numeric(relation_table_out$count_interaction_proxies_expr)

#Write full relation table to file
write.table(relation_table_out, paste0(inp_dir, "all_implications.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

#Make histogram of all implication proxies
make_histogram(c("count_all_proxies"), 
               labels=c("Proxies per Signal/Gene", "Proxies"), ind_sig_df = relation_table_out, output_file = paste0(inp_dir, "all_proxies_per_signal_gene.jpg"))
make_histogram(c("count_all_proxies_expr"), 
               labels=c("Proxies per Signal/Gene", "Proxies"), ind_sig_df = relation_table_out[which(relation_table_out$count_all_proxies_expr>0),], output_file = paste0(inp_dir, "all_proxies_per_signal_expressed_gene.jpg"))

#Append proxy info onto morris_bmd_imp
append_proxy <- function(morris_row, imp_table){
  rsid <- as.character(morris_row["RSID"])
  temp <- imp_table[which(imp_table$SENTINEL == rsid),]
  imp_types <- unique(temp$imp_type)
  all_imps_ints <- ifelse(length(imp_types) == 1, ifelse(imp_types == "interaction", TRUE, FALSE), FALSE) #Indicator variable for whether all implications are via interactions
  proxies <- unique(temp$PROXY)
  count_proxies <- length(proxies)
  locations <- unique(temp$PROXY_POSITION)
  proxy_out <- paste0(proxies, "(", locations, ")", collapse = ", ")
  return(c(all_imps_ints, count_proxies, proxy_out))
}
morris_bmd_proxy <- cbind.data.frame(morris_bmd_imp, 
                                     t(apply(morris_bmd_imp, MARGIN = 1, FUN = append_proxy, imp_table=imp_table[which(imp_table$expressed == TRUE),])))
colnames(morris_bmd_proxy)[(ncol(morris_bmd_proxy)-2):ncol(morris_bmd_proxy)] <- c("all_implication_via_interaction", 
                                                                                   "count_proxies_for_expressed_genes",
                                                                                   "proxies")
morris_bmd_proxy$count_proxies_for_expressed_genes <- as.numeric(morris_bmd_proxy$count_proxies_for_expressed_genes)

#Write morris_bmd_imp to file
write.table(morris_bmd_proxy, paste0(inp_dir, "morris_signals_w_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Make histogram of proxies per signal
make_histogram(c("count_proxies_for_expressed_genes"), 
               labels=c("Proxies per Signal", "Proxies"), ind_sig_df = morris_bmd_proxy[which(morris_bmd_proxy$count_proxies_for_expressed_genes > 0),], output_file = paste0(inp_dir, "all_proxies_per_signal.jpg"))

#Calculate proportion of implication signals with a single proxy
nrow(morris_bmd_proxy[which(morris_bmd_proxy$count_proxies_for_expressed_genes == 1 ),])/nrow(morris_bmd_proxy[which(morris_bmd_proxy$count_proxies_for_expressed_genes > 0),])

############################################################################################################################

#Filter for ideal targets
ideal_target_table <- morris_bmd_proxy[which(morris_bmd_proxy$count_proxies_for_expressed_genes == 1 & 
                                               morris_bmd_proxy$all_implication_via_interaction == TRUE),]
#Write to file
write.table(ideal_target_table, paste0(inp_dir, "single_proxy_expressed_gene_signals.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Calculate number of non-hMSC ideal targets
nrow(ideal_target_table[which(ideal_target_table$count_hMSC_Osteo == 0),])
