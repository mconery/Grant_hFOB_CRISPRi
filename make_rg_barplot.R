################################################################################

#make_rg_barplot.R

#The purpose of this script is to make a bar-plot of the genetic correlation 
#results 

#This code runs with the R/4.2 module

################################################################################


# 0) Call libraries and set variables and directories ====

library(ggplot2)
library(tidyverse)
library(tools)

#Set file locations
global_correlations <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/genetic_correlations/globalRg_820pairs.txt"
plot_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/genetic_correlations/BMD_genetic_correlations.bar.jpeg"
no_second_bmd_plot_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/genetic_correlations/BMD_genetic_correlations.no_panukbb_bmd.bar.jpeg"

# 1) Make Bar-Plot ====

#Read in file
global_raw <- read.table(global_correlations, header = TRUE) 

#Filter for the BMD results and Benjamini-Hochberg adjust the BMD genetic correlation p-values
global_bmd <- global_raw %>% filter(p1=="BMD" | p2 == "BMD") %>%
  mutate(non_bmd_trait=ifelse(p1=="BMD", p2, p1)) %>%
  filter(non_bmd_trait != "PDB") %>% #Remove non-UKBB Paget's Disease
  mutate(p_bh=p.adjust(p, method="BH"), non_bmd_trait = ifelse(non_bmd_trait == "phecode-251.1-both_sexes", "Hypoglycemia", non_bmd_trait)) %>%
  select(non_bmd_trait, rg, p_bh) %>% mutate(v_adj=ifelse(rg>0, 0.4, 1.1))


#Adjust Trait names
global_bmd <- global_bmd %>% mutate(non_bmd_trait_clean=toTitleCase(str_replace_all(non_bmd_trait, "_", " "))) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Alp", "ALP")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Alt", "ALT")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"BF", "Bone Fracture")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Bmi", "BMI")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Ldl", "LDL")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Hdl", "HDL")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Egfr Creat", "eGFR Creatinine")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Ggt", "GGT")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Hba1c", "HbA1c")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Vitamin d", "Vitamin D")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Smoking Ever Never", "Smoking Ever/Never")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Lbs", "(lbs)")) %>% 
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Bone Mineral Density", "BMD (Pan-UKBB)"))

# Create a bar plot
p <- ggplot(global_bmd, aes(x = non_bmd_trait_clean, y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Genetic Correlation") +
  xlab("Phenotype") +
  theme_minimal()

# Add asterisks for significant values (p < 0.05)
p <- p + geom_text(data = subset(global_bmd, p_bh < 0.05), aes(label = "*", vjust = v_adj), size = 10)

#Set min and max values
p <- p + scale_y_continuous(limits = c(-0.9, 0.9))
#Add an axis at x=0
p <- p + geom_hline(yintercept = 0, color="black", linewidth=1)

# Customize the appearance of the asterisks
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
               axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
               axis.title=element_text(size=12), 
               legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 
jpeg(plot_loc, width = 10800, height=4000, res=1000)
print(p)
dev.off()

#Also make a version with Pan-UKBB eBMD removed
global_bmd_filt <- global_bmd %>% filter(non_bmd_trait != "bone_mineral_density")
# Create a bar plot
p <- ggplot(global_bmd_filt, aes(x = non_bmd_trait_clean, y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Genetic Correlation") +
  xlab("Phenotype") +
  theme_minimal()
# Add asterisks for significant values (p < 0.05)
p <- p + geom_text(data = subset(global_bmd_filt, p_bh < 0.05), aes(label = "*", vjust = v_adj), size = 10)
#Set min and max values
p <- p + scale_y_continuous(limits = c(-0.9, 0.9))
#Add an axis at x=0
p <- p + geom_hline(yintercept = 0, color="black", linewidth=1)
# Customize the appearance of the asterisks
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
               axis.text.x = element_text(size = 12, angle=45, hjust = 1, vjust = 1), 
               axis.title=element_text(size=12), 
               legend.position = "none", axis.text.y = element_text(color="black", size=12), axis.title.x = element_blank()) 
jpeg(no_second_bmd_plot_loc, width = 10800, height=4000, res=1000)
print(p)
dev.off()