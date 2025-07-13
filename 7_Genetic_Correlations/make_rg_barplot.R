################################################################################

#make_rg_barplot.R

#The purpose of this script is to make a bar-plot of the genetic correlation 
#results and any related supplementary tables.

#This code runs with R/4.2 on my local.

################################################################################


# 0) Call libraries and set variables and directories ====

library(ggplot2)
library(tidyverse)
library(tools)

#Set file locations
data_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/genetic_correlations/"
global_correlations <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/data/genetic_correlations/gR_alltraits.txt"
plot_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/genetic_correlations/BMD_genetic_correlations.bar.jpeg"
no_second_bmd_plot_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/genetic_correlations/BMD_genetic_correlations.no_panukbb_bmd.bar.jpeg"
impede_loc <- "C:/Users/mitch/Documents/UPenn/Grant_Lab/hFOB_CRISPRi_Screen/figures/genetic_correlations/impedance_traits_genetic_correlations.heatmap.jpeg"


# 1) Make Main Figure Bar-Plot ====

#Read in file
global_raw <- read.table(global_correlations, header = TRUE) 

#Filter for the BMD results and Benjamini-Hochberg adjust the BMD genetic correlation p-values
global_bmd <- global_raw %>% filter(p1=="BMD" | p2 == "BMD") %>%
  mutate(non_bmd_trait=ifelse(p1=="BMD", p2, p1)) %>%
  filter(non_bmd_trait != "PDB") %>% #Remove non-UKBB Paget's Disease
  filter(!(is.na(rg))) %>% #NA values result from not enough signal in the data to estimate the genetic variance, so we'll remove these signals
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
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Bone Mineral Density", "BMD (Pan-UKBB)")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean, pattern = "Whole Body", replacement = "Whole-Body"))

# Create a bar plot
p <- ggplot(global_bmd, aes(x = non_bmd_trait_clean, y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Genetic Correlation with BMD") +
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
#Add annotation
p <- p + annotate(
  "text",
  x = Inf, y = Inf,
  label = expression("* BH-Adj. " * italic(P) * " < 0.05"),
  hjust = 1.1, vjust = 1.5,
  size = 6,
  color = "black"
)
jpeg(plot_loc, width = 10800, height=4000, res=1000)
print(p)
dev.off()

#Also make a version with Pan-UKBB eBMD removed
global_bmd_filt <- global_bmd %>% filter(non_bmd_trait != "bone_mineral_density")
# Create a bar plot
p <- ggplot(global_bmd_filt, aes(x = non_bmd_trait_clean, y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Genetic Correlation with BMD") +
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
#Add annotation
p <- p + annotate(
  "text",
  x = Inf, y = Inf,
  label = expression("* BH-Adj. " * italic(P) * " < 0.05"),
  hjust = 1.1, vjust = 1.5,
  size = 6,
  color = "black"
)
jpeg(no_second_bmd_plot_loc, width = 10800, height=4000, res=1000)
print(p)
dev.off()

#Create tiff image
p <- ggplot(global_bmd_filt, aes(x = non_bmd_trait_clean, y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Genetic Correlation with BMD") +
  xlab("Phenotype") +
  theme_minimal()
# Add asterisks for significant values (p < 0.05)
p <- p + geom_text(data = subset(global_bmd_filt, p_bh < 0.05), aes(label = "*", vjust = v_adj), size = 6)
#Set min and max values
p <- p + scale_y_continuous(limits = c(min(1.1*min(global_bmd_filt$rg),-0.5), max(0.5,1.1*max(global_bmd_filt$rg))))
#Add an axis at x=0
p <- p + geom_hline(yintercept = 0, color="black", linewidth=0.3)
# Customize the appearance of the asterisks
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_rect(fill="white", color = "white"), legend.title = element_blank(),
               axis.text.x = element_text(size = 7, angle=45, hjust = 1, vjust = 1, color = "black", family = "Arial"),
               plot.background = element_rect(color = "white"),
               axis.line = element_line(color="black", linewidth = 0.3),
               axis.title=element_text(size=7), 
               legend.position = "none", 
               axis.text.y = element_text(color="black", size=7, family = "Arial"), axis.title.x = element_blank()) 
#Add annotation
p <- p + annotate(
  "text",
  x = Inf, y = Inf,
  label = expression("* BH-Adj. " * italic(P) * " < 0.05"),
  hjust = 1.1, vjust = 1.5,
  size = 3,
  color = "black"
)
ggsave(str_replace(no_second_bmd_plot_loc, pattern = ".jpeg", replacement = ".tif"), 
       width = 6, height=2.4, device = "tif", units = "in", dpi=1000, plot = p)

#Write table to file for use in supplement
global_raw %>% filter(p1=="BMD" | p2 == "BMD") %>%
  mutate(non_bmd_trait=ifelse(p1=="BMD", p2, p1)) %>%
  filter(non_bmd_trait != "PDB") %>% #Remove non-UKBB Paget's Disease
  mutate(p_bh=p.adjust(p, method="BH"), non_bmd_trait = ifelse(non_bmd_trait == "phecode-251.1-both_sexes", "Hypoglycemia", non_bmd_trait)) %>%
  dplyr::select(non_bmd_trait, rg, se, p, p_bh) %>% mutate(v_adj=ifelse(rg>0, 0.4, 1.1)) %>% 
  mutate(non_bmd_trait_clean=toTitleCase(str_replace_all(non_bmd_trait, "_", " "))) %>% 
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
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean,"Bone Mineral Density", "BMD (Pan-UKBB)")) %>%
  mutate(non_bmd_trait_clean=str_replace(non_bmd_trait_clean, pattern = "Whole Body", replacement = "Whole-Body")) %>% 
  dplyr::select(non_bmd_trait_clean, rg, se, p, p_bh) %>% 
  write.table(file=paste0(data_loc, "bmd_genetic_correlations.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# 2) Make a Matrix of Genetic-correlations Between the Impedance-Derived Traits ====

#Filter for impedance/weight traits
impede_weight_traits <- c("basal_metabolic_rate", "body_fat_percentage", "impedance_whole_body", 
                          "trunk_fat_percentage", "weight", "whole_body_fat-free_mass", "whole_body_fat_mass", "whole_body_water_mass", "bmi")
global_impede_weight <- global_raw %>% 
  dplyr::filter(p1 %in% impede_weight_traits & p2 %in% impede_weight_traits) %>%
  mutate(p_bh=p.adjust(p, method="BH")) %>% 
  mutate(p1=toTitleCase(str_replace_all(p1, "_", " "))) %>%
  mutate(p2=toTitleCase(str_replace_all(p2, "_", " "))) %>%
  mutate(p1=str_replace(p1,"Bmi", "BMI")) %>% 
  mutate(p2=str_replace(p2,"Bmi", "BMI")) %>% 
  mutate(p1=str_replace(p1, pattern = "Whole Body", replacement = "Whole-Body")) %>% 
  mutate(p2=str_replace(p2, pattern = "Whole Body", replacement = "Whole-Body"))

#Output plot
jpeg(impede_loc, width = 10800, height=10800, res=1000)
ggplot(global_impede_weight, aes(x = p2, y = p1)) +
  scale_y_discrete(position = "right") + 
  geom_tile(aes(fill = rg)) +
  geom_text(aes(label = round(rg, 2)), size = 8) +
  scale_fill_gradient(low = "blue", high = "red") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=-45, hjust = 0, vjust = 1), 
        axis.title=element_blank(), 
        legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title.x = element_blank()) 
dev.off()
tiff(str_replace(impede_loc, pattern = ".jpeg", replacement = ".tiff"), width = 10800, height=10800, res=1000)
ggplot(global_impede_weight, aes(x = p2, y = p1)) +
  scale_y_discrete(position = "right") + 
  geom_tile(aes(fill = rg)) +
  geom_text(aes(label = round(rg, 2)), size = 8) +
  scale_fill_gradient(low = "blue", high = "red") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 24, angle=-45, hjust = 0, vjust = 1), 
        axis.title=element_blank(), 
        legend.position = "none", axis.text.y = element_text(color="black", size=24), axis.title.x = element_blank()) 
dev.off()

#Write table to file
global_impede_weight %>% dplyr::select(p1, p2, rg, se, p, p_bh) %>% 
  write.table(file=paste0(data_loc, "impedance_weight_genetic_correlations.supplement.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

