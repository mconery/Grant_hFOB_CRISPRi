################################################################################

#make_s-ldsc_enrichment_plots.R

#The purpose of this script is to make enrichment plots showing the level of
#heritability enrichment in the activating chromatin regions across each of the
#tissue types included in our S-LDSC Analysis

#This code runs with the R/4.2 module

################################################################################


# 0) Call libraries and set variables and directories ====

#Call libraries
library(tidyverse)
library(scales)
library(ggtext)
library(glue)

#Set directory
work_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/roadmap_epigenomics_chip/s-ldsc_data/"
plot_dir <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/s-ldsc_results/"
#Set files
cell_types_file=paste0(work_dir, "complete_cell_types.tsv")
bmd_enrich_file=paste0(work_dir, "summary_BMD.txt")
bf_enrich_file=paste0(work_dir, "summary_BF.txt")
pdb_enrich_file=paste0(work_dir, "summary_PDB.txt")

# 1) Read in files (and reformat PDB file) ====

#Read in files
cell_types_raw <- read.table(cell_types_file, header = TRUE, sep = "\t")
bmd_enrich_raw <- read.table(bmd_enrich_file, header = TRUE, sep = "\t")
bf_enrich_raw <- read.table(bf_enrich_file, header = TRUE, sep = "\t")
pdb_enrich_raw <- read.table(pdb_enrich_file, header = TRUE, sep = "\t")

#Modify format of cell-types file to work with function
colnames(cell_types_raw) <- str_replace(colnames(cell_types_raw), "bedfile", "Celltype")

#Map new cell types in new enrichment files to original cell types
map_cell_types <- function(enrich_raw){
  enrich_raw$Celltype <- enrich_raw$Celltype %>%
    str_replace(pattern = "adult_hMSC_Osteoblasts_atac", replacement = "hMSC_BMP2") %>%
    str_replace(pattern = "chondrocytes_atac", replacement = "Chondrocytes.Grant.OCRs") %>% 
    str_replace(pattern = "hFOB_Diff_atac", replacement = "hFOB_Diff") %>% 
    str_replace(pattern = "hFOB_Perm_atac", replacement = "hFOB_Perm") %>% 
    str_replace(pattern = "hFOBs_chip", replacement = "hFOB.H3K27ac.Cottone") %>% 
    str_replace(pattern = "monocytes_atac", replacement = "Monocytes.Grant.OCRs") %>% 
    str_replace(pattern = "osteoclasts_bae_2022_atac", replacement = "OSC.RANKL.Diff.Atac.Bae") %>% 
    str_replace(pattern = "osteoclasts_Bae_2022_chip", replacement = "OSC.RANKL.Diff.H3K27ac.Bae") %>% 
    str_replace(pattern = "osteoclasts_Grant_D4_atac", replacement = "OSC.Day4.Diff.Atac.Grant") %>% 
    str_replace(pattern = "pediatric_hMSCs_Control_atac", replacement = "Control_3D") %>% 
    str_replace(pattern = "pediatric_hMSCs_Osteoblasts_3D_atac", replacement = "BMP2_3D") %>% 
    str_replace(pattern = "pediatric_hMSCs_Osteoblasts_6D_atac", replacement = "BMP2_6D")
    return(enrich_raw)
}
pdb_enrich_raw <- map_cell_types(pdb_enrich_raw)
bf_enrich_raw <- map_cell_types(bf_enrich_raw)
bmd_enrich_raw <- map_cell_types(bmd_enrich_raw)

# 2) Make Enrichment Plots ====

#Make a thin version of the cell-types dataframe
cell_types_thin <- cell_types_raw %>% select(Celltype, Cell.type.group, E.Mnemonic, MARK)

#Create a function to make plots for each enrichment
make_enrichment_plot <- function(enrich_raw, cell_types_thin, plot_dir, file_prefix){
  #Append on cell types to heritability enrichment data and adjust p-values
  enrich_append <- enrich_raw %>% inner_join(cell_types_thin, by="Celltype") %>%
    mutate(bonf_adj_pval = p.adjust(Enrichment_p, method="bonferroni"), 
           neg_log_adj_pval=(-log10(p.adjust(Enrichment_p, method="bonferroni"))))
  write.table(enrich_append, paste0(plot_dir, file_prefix, ".with_tissues_bonf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  #Check if E.Mnemonic column was already present
  if ("E.Mnemonic.x" %in% colnames(enrich_append)) {
    colnames(enrich_append) <- str_replace_all(colnames(enrich_append), "[.]x", "")
  }
  #Make cell-type grouped box-plot
  enrich_append$Cell.type.group <- as.factor(enrich_append$Cell.type.group)
  type_violin <- ggplot(enrich_append, aes(x=Cell.type.group, y=neg_log_adj_pval, color=Cell.type.group)) + geom_boxplot() + 
    ylab("-log(Bonferroni P-Value)") + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 2) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(color="black", size=12, angle=90, hjust = 1, vjust = 0.5), axis.title=element_text(size=12), 
          legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 
  jpeg(paste0(plot_dir,file_prefix, ".grouped_cell_type.violin.jpeg"), width = 7200, height=4800, res = 1000)
    print(type_violin)
  dev.off()
  #Prep data frame for making cell-type separated bar plot
  enrich_append <- enrich_append %>% arrange(Cell.type.group, E.Mnemonic, MARK)
  enrich_append$Cell.type.group <- as.factor(enrich_append$Cell.type.group)
  enrich_append$MARK <- as.factor(enrich_append$MARK)
  enrich_append$E.Mnemonic <- as.factor(enrich_append$E.Mnemonic)
  enrich_append$Celltype <- factor(enrich_append$Celltype, levels = enrich_append$Celltype)
  #Calculate break points for bar plot axis labels and append them to the enrich_append dataframe
  temp <- enrich_append %>% group_by(Cell.type.group) %>% dplyr::summarize(count=n()) %>% mutate(axis_break=floor(cumsum(count)-count/2))
  temp <- temp %>% mutate(discrete_break=enrich_append$Celltype[axis_break], label_color = hue_pal()(nrow(temp)))
  temp2 <- temp$label_color
  names(temp2) <- temp$discrete_break
  enrich_append <- enrich_append %>% mutate(Celltype_num = seq(1, nrow(enrich_append), 1)) %>%
    mutate(axis_label=ifelse(Celltype %in% temp$discrete_break, as.character(Cell.type.group), as.character(Celltype_num)),
    color = ifelse(Celltype %in% temp$discrete_break, as.character(temp2[as.character(Celltype)]), "#FFFFFF00")) %>% #Use 8-digit hex code for transparent text
    mutate(name = glue("<i style='color:{color}'>{axis_label}"))
  #Reorder the factors for the axis labels (name field)
  enrich_append$name <- factor(enrich_append$name, levels = enrich_append$name)
  #Make bar plot
  separate_bar <- ggplot(enrich_append, aes(x=name, y=neg_log_adj_pval, fill=Cell.type.group)) + geom_bar(stat="identity") + 
    ylab("-log(Bonferroni P-Value)") + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 2) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_markdown(size = 14, angle=90, hjust = 1, vjust = 0.75), 
          axis.title=element_text(size=14), 
          legend.position = "none", axis.text.y = element_text(color="black", size=14), axis.title.x = element_blank()) 
  jpeg(paste0(plot_dir,file_prefix, ".split_cell_type.bar.jpeg"), width = 10800, height=4000, res=1000)
  print(separate_bar)
  dev.off()
}
#Call function
make_enrichment_plot(bmd_enrich_raw, cell_types_thin=cell_types_thin, plot_dir=plot_dir, "bmd")
make_enrichment_plot(bf_enrich_raw, cell_types_thin=cell_types_thin, plot_dir=plot_dir, "fracture")
make_enrichment_plot(pdb_enrich_raw, cell_types_thin=cell_types_thin, plot_dir=plot_dir, "pagets_disease")
