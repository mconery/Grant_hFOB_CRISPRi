################################################################################

#coloc_BMD_pairwise.R

#This script is designed to execute coloc on the output of a multi-trait fine-
#mapping experiment of BMD-associated loci. The script executes pairwise 
#colocalizations between BMD and all the other mapped traits.

################################################################################

# 0) Load Required Packages ----
suppressPackageStartupMessages({
  library(coloc)
  library(stringr)
  library(data.table)
  library(dplyr)
})

# 1) Parse Command Line Arguments ----
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Usage: Rscript coloc_pairwise.R <locus_prefix> <rds_dir> <bmd_trait> <out_dir> [--p1 1e-4] [--p2 1e-4] [--p12 1e-5] [--p4_cut 0.5]")
}

# Required arguments
locus_prefix <- args[1]    # e.g. "chr6.12345678.12355678"
rds_dir <- args[2]         # Directory containing .rds files
bmd_trait <- args[3]       # Name of BMD trait (e.g. "eBMD")
out_dir <- args[4]         # Output directory

#Test locations
locus_prefix <- "chr12.131250001.132000000"
rds_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/susie_results"
bmd_trait <- "BMD"
out_dir <- "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/susie_coloc_results"
# Optional arguments with defaults
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5
p4_cut <- 0.5

# Parse optional parameters
if(length(args) > 4) {
  opt_args <- args[5:length(args)]
  for(i in seq_along(opt_args)) {
    if(opt_args[i] == "--p1") p1 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p2") p2 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p12") p12 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p4_cut") p4_cut <- as.numeric(opt_args[i+1])
  }
}

# 2) File Validation ----
check_file <- function(path) {
  if(!file.exists(path)) {
    stop(paste("File does not exist:", path))
  }
}

# Get all RDS files for this locus
pattern <- paste0(locus_prefix, "\\.susie\\.rds$")
rds_files <- list.files(rds_dir, pattern = pattern, full.names = TRUE)

if(length(rds_files) < 2) {
  stop(paste("Insufficient RDS files found for locus:", locus_prefix))
}

# Verify BMD file exists
bmd_file <- grep(paste0(bmd_trait, "\\."), rds_files, value = TRUE)
if(length(bmd_file) != 1) {
  stop(paste("Could not find unique BMD file for", bmd_trait))
}

# 3) Load Data ----
message("\nLoading SuSiE results...")
susie_objects <- lapply(rds_files, function(f) {
  obj <- readRDS(f)
  if(!inherits(obj, "susie")) {
    stop(paste("Invalid SuSiE object in file:", f))
  }
  return(obj)
})
names(susie_objects) <- str_remove(str_remove(basename(rds_files), "\\.susie\\.rds$"), paste0("\\.",locus_prefix))

# 4) Prepare for Coloc ----
message("\nPreparing coloc analysis...")
bmd_obj <- susie_objects[[grep(bmd_trait, names(susie_objects))]]
other_traits <- setdiff(names(susie_objects), bmd_trait)

# 5) Run Pairwise Coloc ----
results <- list()

for(trait in other_traits) {
  message(paste("\nColocalizing", bmd_trait, "vs", trait))
  
  trait_obj <- susie_objects[[trait]]
  
  # Perform coloc analysis
  coloc_res <- coloc.susie(
    dataset1 = bmd_obj,
    dataset2 = trait_obj,
    p1 = p1,
    p2 = p2,
    p12 = p12
  )
  
  results[[trait]] <- coloc_res
}

# [Remaining code for output formatting and saving would go here]

# 6) Extract BMD information into a Table ====

#Extract BMD Signals and send to a table
bmd_sets <- bmd_obj$sets$cs_index
extract_cs <- function(index, bmd_object=bmd_obj){bmd_object$sets$cs[[paste0("L",index)]]}
extract_field <- function(variant_indices, field ,bmd_object=bmd_obj){bmd_object[[field]][variant_indices]}
bmd_table <- cbind.data.frame(signal_id=paste(locus_prefix, bmd_sets, sep = "."), 
                              locus=rep(locus_prefix, length(bmd_sets)),
                              signal=bmd_sets,
                              set_size=vapply(lapply(bmd_sets, FUN=extract_cs), FUN = length, FUN.VALUE = numeric(1)),
                              snp_ids=vapply(lapply(lapply(bmd_sets, FUN=extract_cs), FUN = names), FUN = paste0, FUN.VALUE = character(1), collapse = ","),
                              rsids=vapply(lapply(lapply(bmd_sets, FUN=extract_cs), extract_field, field="rsid"), FUN = paste0, FUN.VALUE = character(1), collapse = ","),
                              purity=bmd_obj$sets$purity[paste0("L", bmd_sets),"min.abs.corr"],
                              max_neglogp=vapply(lapply(lapply(bmd_sets, FUN=extract_cs), extract_field, field="neglogp_gwas"), FUN=max, FUN.VALUE = numeric(1)),
                              max_pip=vapply(lapply(lapply(bmd_sets, FUN=extract_cs), extract_field, field="pip"), FUN=max, FUN.VALUE = numeric(1)))

# 7) Extract other trait info into the table ====

#Create a function to extract all the desired results
extract_coloc <- function(non_bmd_trait, results_list=results){
  result_df = results_list[[non_bmd_trait]]$summary
  if (!is.null(result_df)) {
    good_coloc_rows = result_df %>% dplyr::filter(PP.H4.abf > p4_cut)
    max_neglogp <- vapply(lapply(lapply(good_coloc_rows$idx2, FUN = extract_cs, bmd_object=susie_objects[[non_bmd_trait]]), extract_field, bmd_object=susie_objects[[non_bmd_trait]], field="neglogp_gwas"), FUN = max, FUN.VALUE = numeric(1))
    max_pip <- vapply(lapply(lapply(good_coloc_rows$idx2, FUN = extract_cs, bmd_object=susie_objects[[non_bmd_trait]]), extract_field, bmd_object=susie_objects[[non_bmd_trait]], field="pip"), FUN = max, FUN.VALUE = numeric(1))
    return_df <- cbind.data.frame(signal=good_coloc_rows$idx1,
                                  other_traits=rep(non_bmd_trait, nrow(good_coloc_rows)), 
                                  PP4=good_coloc_rows$PP.H4.abf,
                                  other_traits_max_pip=max_pip,
                                  other_traits_max_neglogp=max_neglogp)
    return(return_df)
  } else{
    return_df <- cbind.data.frame(signal=NA,
                                  other_traits=NA, 
                                  PP4=NA,
                                  other_traits_max_pip=NA,
                                  other_traits_max_neglogp=NA)
    return(return_df)
  }
}
#Extract the results
summarized_results <- na.omit(dplyr::bind_rows(lapply(other_traits, extract_coloc))) %>% group_by(signal) %>%
  summarize(num_other_traits=n(), other_traits=paste0(unique(other_traits), collapse = ","), 
            PP4=paste0(PP4, collapse = ","), other_traits_max_pip=paste0(other_traits_max_pip, collapse = ","), other_traits_max_neglogp=paste0(other_traits_max_neglogp, collapse = ",")) %>%
  as.data.frame

#Merge the dataframes
export_table <- bmd_table %>% left_join(summarized_results, by = "signal")

#Write to file
write.table(export_table, file = file.path(out_dir, paste0(locus_prefix, ".susie-coloc.tsv")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
