# Grant hFOB CRISPRi Screen and Associated Analyses
This Github is a repository for all the code used to complete the analyses
presented in "GWAS-informed data integration and non-coding CRISPRi screen 
reveals etiological insights for bone mineral density" (Conery ###et al.###). 
It is organized into the following eight sections:

### 1) 1_ATAC-seq_Processing ###
This folder contains the ENCODE configuration files necessary to call 
peaks for the new and re-processed ATAC/ChiP-seq data sets used in the S-LDSC 
analysis and CRISPRi target selection. It also contains the scripts used 
to download the reprocessed FASTQ files and the files necessary to check the 
quality of the osteoclast datasets.

### 2) 2_S-LDSC_Heritability_Enrichment ###
This folder contains the scripts used for executing the S-LDSC analysis which
includes the script used to munge the summary statistics used in the S-LDSC
and cross-trait LDSC analyses. The following figures were generated by the 
code in this section:

Fig. 1 - **make_s-ldsc_enrichment_plots.R**  
Supplementary Fig. 1 - **make_s-ldsc_enrichment_plots.R** 

### 3) 3_hFOB_Capture-C_and_RNA-seq ###
This folder contains the code used to generate the new hFOB RNA-seq results 
and the Chicago config file used to create the Capture-C dataset.

### 4) 4_Target_Selection_and_sgRNA_Design ###
This folder conains the scripts used to select the screen targets from the 
hMSC-Osteoblast and hFOB RNA-seq, ATAC-seq, and Capture-C results. It also 
contains the scripts and commands used to design the sgRNAs. The following 
figures were generated by the code in this section:

Fig. 2a - **compare_bone_cell_pcc.R**
Fig. 2d - **run_pygenome.sh**
Fig. 2e - **run_pygenome.sh**
Supplementary Fig. 11 - **run_pygenome.sh**
Supplementary Fig. 12 - **run_pygenome.sh**
Supplementary Fig. 13 - **run_pygenome.sh**

### 5) 5_Single-Cell_QC_and_Sceptre ###
The folder contains the scripts used to quality filter the scRNA-seq screen 
readouts and the SCEPTRE analysis used to test for significant perturbations. 
The following figures were generated by the code in this section:

Fig. 2b - **sceptre_differential_expression.R**
Fig. 2c - **sceptre_differential_expression.R**
Supplementary Fig. 3 - **crispri_screen_qc.py**
Supplementary Fig. 4 - **crispri_screen_qc.py**
Supplementary Fig. 5 - **crispri_screen_qc.py**
Supplementary Fig. 6 - **sceptre_differential_expression.R**

### 6) 6_Osteoblast_Assays ###
This folder contains the script that tests for significant effects of siRNA
knockdown in the osteoblast and adipogenesis assays. The following figure was
generated by the code in this section:

Fig. 3 - **make_assay_figures.R**

### 7) 7_Genetic_Correlations ###
This folder contains the codes for executing the cross-trait LDSC analysis for
genetic correlations. The following figures were generated by the code in this
section:

Fig. 4 - **make_rg_barplot.R**
Supplementary Fig. 14 - **make_rg_barplot.R**

### 8) 8_Multi-Trait_Fine-Mapping ###
This folder contains all the scripts needed to execute the cross-trait fine-
mapping experiment and analyze the full results. The following figures are 
generated by the code in this section

\Fig. 5 - 4_analyze_BMD_signals/**cluster_signals_and_check_overlap.R**
\Supplementary Fig. 15 - 4_analyze_BMD_signals/**cluster_signals_and_check_overlap.R**
\Supplementary Fig. 16 - 4_analyze_BMD_signals/**cluster_signals_and_check_overlap.R**
\Supplementary Fig. 17 - 4_analyze_BMD_signals/**cluster_signals_and_check_overlap.R**
\Supplementary Fig. 18 - 4_analyze_BMD_signals/**cluster_signals_and_check_overlap.R**


\*Some of the figures were manually edited after being created in the scripts above
to include or adjust labels
\\\*\*Other figures not covered by the above scripts were either created from 
photographs or manually drawn
