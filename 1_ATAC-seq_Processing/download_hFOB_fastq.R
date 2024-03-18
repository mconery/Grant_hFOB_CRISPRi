################################################################################

#download_hFOB_fastq.R

#The purpose of this script is to download the H3K27ac ChIP-seq hFOB FASTQ 
#files. These will be used to generate peak calls and ultimately hFOB-relevant
#fine-mapped BMD variants. The script runs with the following modules:

#  R/4.2.3

################################################################################

# 0) Load libraries ====

library(GEOfastq)

# 1) Set file locations and key variables ====

download_dir = "/mnt/isilon/sfgi/conerym/data/ChIP-seq/hFOBs/cottone_2022/"
gse_name <- 'GSE152942'

# 2) Download fastq files ====

#Set working directory
setwd(download_dir)
#Get text list needed for GEOfastq corresponding to accession number
gse_text <- crawl_gse(gse_name)
#Get names of all sub accessions for the main accession
gsm_names <- extract_gsms(gse_text)
srp_meta <- crawl_gsms(gsm_names)
#Filter for runs that are either WT input or H3K27ac for runs 1-3
srp_meta <- srp_meta[c(grep(pattern = "hFOB_WT[1-3]_Input", x = srp_meta$title),
                       grep(pattern = "hFOB_WT[1-3]_H3K27ac", x = srp_meta$title)),]
#Download files
res <- get_fastqs(srp_meta, download_dir)