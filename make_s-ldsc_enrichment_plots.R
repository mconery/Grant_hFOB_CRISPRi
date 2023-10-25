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
bmd_enrich_file=paste0(work_dir, "ldsc_BMD.txt")
bf_enrich_file=paste0(work_dir, "ldsc_BF.txt")

# 1) Read in files ====

cell_types_raw <- read.table(cell_types_file, header = TRUE, sep = "\t")
bmd_enrich_raw <- read.table(bmd_enrich_file, header = TRUE, sep = "\t")
bf_enrich_raw <- read.table(bf_enrich_file, header = TRUE, sep = "\t")

# 2) Make Enrichment Plots ====

#Make a thin version of the cell-types dataframe
cell_types_thin <- cell_types_raw %>% select(bedfile, Cell.type.group)

#Create a function to make plots for each enrichment
make_enrichment_plot <- function(enrich_raw, cell_types_thin, plot_dir, file_prefix){
  #Append on cell types to heritability enrichment data and adjust p-values
  enrich_append <- enrich_raw %>% inner_join(cell_types_thin, by="bedfile") %>%
    mutate(neg_log_adj_pval=(-log10(p.adjust(Enrichment_p, method="bonferroni"))))
  #Make cell-type grouped box-plot
  enrich_append$Cell.type.group <- as.factor(enrich_append$Cell.type.group)
  type_violin <- ggplot(enrich_append, aes(x=Cell.type.group, y=neg_log_adj_pval, color=Cell.type.group)) + geom_boxplot() + 
    ylab("-log(Bonferroni P-Value)") + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 2) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(color="black", size=16, angle=90, hjust = 1, vjust = 0.5), axis.title=element_text(size=16), 
          legend.position = "none", axis.text.y = element_text(color="black", size=16), axis.title.x = element_blank()) 
  jpeg(paste0(plot_dir,file_prefix, ".grouped_cell_type.violin.jpeg"), width = 720, height=480)
    print(type_violin)
  dev.off()
  #Prep data frame for making cell-type separated bar plot
  enrich_append <- enrich_append %>% arrange(Cell.type.group, E.Mnemonic, MARK)
  enrich_append$Cell.type.group <- as.factor(enrich_append$Cell.type.group)
  enrich_append$MARK <- as.factor(enrich_append$MARK)
  enrich_append$E.Mnemonic <- as.factor(enrich_append$E.Mnemonic)
  enrich_append$bedfile <- factor(enrich_append$bedfile, levels = enrich_append$bedfile)
  #Calculate break points for bar plot axis labels and append them to the enrich_append dataframe
  temp <- enrich_append %>% group_by(Cell.type.group) %>% summarize(count=n()) %>% mutate(axis_break=floor(cumsum(count)-count/2))
  temp <- temp %>% mutate(discrete_break=enrich_append$bedfile[axis_break], label_color = hue_pal()(nrow(temp)))
  temp2 <- temp$label_color
  names(temp2) <- temp$discrete_break
  enrich_append <- enrich_append %>% mutate(bedfile_num = seq(1, nrow(enrich_append), 1)) %>%
    mutate(axis_label=ifelse(bedfile %in% temp$discrete_break, as.character(Cell.type.group), as.character(bedfile_num)),
    color = ifelse(bedfile %in% temp$discrete_break, as.character(temp2[as.character(bedfile)]), "#FFFFFF00")) %>% #Use 8-digit hex code for transparent text
    mutate(name = glue("<i style='color:{color}'>{axis_label}"))
  #Reorder the factors for the axis labels (name field)
  enrich_append$name <- factor(enrich_append$name, levels = enrich_append$name)
  #Make bar plot
  separate_bar <- ggplot(enrich_append, aes(x=name, y=neg_log_adj_pval, fill=Cell.type.group)) + geom_bar(stat="identity") + 
    ylab("-log(Bonferroni P-Value)") + 
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", linewidth = 2) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_markdown(size = 20, angle=90, hjust = 1, vjust = 0.75), 
          axis.title=element_text(size=20), 
          legend.position = "none", axis.text.y = element_text(color="black", size=16), axis.title.x = element_blank()) 
  jpeg(paste0(plot_dir,file_prefix, ".split_cell_type.bar.jpeg"), width = 1080, height=400)
  print(separate_bar)
  dev.off()
}
#Call function
make_enrichment_plot(bmd_enrich_raw, cell_types_thin=cell_types_thin, plot_dir=plot_dir, "bmd")
make_enrichment_plot(bf_enrich_raw, cell_types_thin=cell_types_thin, plot_dir=plot_dir, "fracture")
