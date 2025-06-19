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
  library(jsonlite)
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
trait_file <- "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/trait_sample_sizes_cc.tsv"
locus_file <- "/mnt/isilon/sfgi/conerym/analyses/grant/multi-trait_fine-mapping/bmd_and_related/loci_files/traits_per_loci.json"

# Parse optional parameters
if(length(args) > 4) {
  opt_args <- args[5:length(args)]
  for(i in seq_along(opt_args)) {
    if(opt_args[i] == "--p1") p1 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p2") p2 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p12") p12 <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--p4_cut") p4_cut <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--trait_file") trait_file <- as.numeric(opt_args[i+1])
    if(opt_args[i] == "--locus_file") locus_file <- as.numeric(opt_args[i+1])
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

#Read in the json file
locus_json <- fromJSON(locus_file)
#Extract the traits for the locus
locus_traits <- locus_json[[locus_prefix]]
#Get traits corresponding to rds files
rds_traits <- str_remove(str_remove(basename(rds_files), "\\.susie\\.rds$"), paste0("\\.",locus_prefix))

if(length(rds_files) < 2) {
  stop(paste("ERROR: Insufficient RDS files found for locus:", locus_prefix))
} else if (any(sort(rds_traits) != sort(locus_traits))){
  paste("ERROR: At least one trait did not complete SuSiE fine-mapping for locus:", locus_prefix)
  quit(save="no")
}

# Verify BMD file exists
bmd_file <- grep(paste0(bmd_trait, "\\."), rds_files, value = TRUE)
if(length(bmd_file) != 1) {
  stop(paste("ERROR: Could not find unique BMD file for", bmd_trait))
}

# 3) Load Data ----
message("\nLoading SuSiE results...")
susie_objects <- lapply(rds_files, function(f) {
  obj <- readRDS(f)
  if(!inherits(obj, "susie")) {
    paste("WARNING: NULL object in file:", f)
    return()
  } else{
    return(obj)
  }
})
names(susie_objects) <- str_remove(str_remove(basename(rds_files), "\\.susie\\.rds$"), paste0("\\.",locus_prefix))

#Remove null objects
susie_objects <- susie_objects[!unlist(lapply(susie_objects, is.null))]
#Check if all mappings failed
if (length(susie_objects) == 0 || !("BMD" %in% names(susie_objects))){
  print(paste0("WARNING: ", locus_prefix, " was not mapped for any traits"))
  quit(save = "no")
}

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

# 7) Extract summary of signals into a table of export ====

#Create a function to extract table values
extract_signs <- function(index_pair, check_obj, bmd_object=bmd_obj, field="mu"){
  cs <- names(extract_cs(index_pair[1], bmd_object))
  max_pip_snp <- names(which.max(bmd_object$pip[cs])) #Get max pip snp
  bmd_effect <- check_obj[[field]][index_pair[1],max_pip_snp]
  check_effect <- check_obj[[field]][index_pair[2],max_pip_snp]
  return(ifelse(bmd_effect * check_effect > 0, "+", "-"))
}

#Create a function to extract all the desired results
extract_coloc <- function(non_bmd_trait, results_list=results, susie_objs=susie_objects){
  result_df = results_list[[non_bmd_trait]]$summary
  if (!is.null(result_df)) {
    good_coloc_rows = result_df %>% dplyr::filter(PP.H4.abf > p4_cut)
    max_neglogp <- vapply(lapply(lapply(good_coloc_rows$idx2, FUN = extract_cs, bmd_object=susie_objs[[non_bmd_trait]]), extract_field, bmd_object=susie_objs[[non_bmd_trait]], field="neglogp_gwas"), FUN = max, FUN.VALUE = numeric(1))
    max_pip <- vapply(lapply(lapply(good_coloc_rows$idx2, FUN = extract_cs, bmd_object=susie_objs[[non_bmd_trait]]), extract_field, bmd_object=susie_objs[[non_bmd_trait]], field="pip"), FUN = max, FUN.VALUE = numeric(1))
    signs <- apply(good_coloc_rows[,c("idx1", "idx2")], MARGIN = 1, FUN = extract_signs, check_obj=susie_objects[[non_bmd_trait]])
    return_df <- cbind.data.frame(signal=good_coloc_rows$idx1,
                                  other_traits=rep(non_bmd_trait, nrow(good_coloc_rows)), 
                                  PP4=good_coloc_rows$PP.H4.abf,
                                  other_traits_max_pip=max_pip,
                                  other_traits_max_neglogp=max_neglogp,
                                  other_traits_max_pip_sign=signs)
    return(return_df)
  } else{
    return_df <- cbind.data.frame(signal=NA,
                                  other_traits=NA, 
                                  PP4=NA,
                                  other_traits_max_pip=NA,
                                  other_traits_max_neglogp=NA,
                                  other_traits_max_pip_sign=NA)
    return(return_df)
  }
}
#Extract the results
summarized_results <- na.omit(dplyr::bind_rows(lapply(other_traits, extract_coloc))) %>% group_by(signal) %>%
  summarize(num_other_traits=n(), other_traits=paste0(unique(other_traits), collapse = ","), 
            PP4=paste0(PP4, collapse = ","), other_traits_max_pip=paste0(other_traits_max_pip, collapse = ","), 
            other_traits_max_neglogp=paste0(other_traits_max_neglogp, collapse = ","),
            other_traits_max_pip_signs=paste0(other_traits_max_pip_sign, collapse = ",")) %>%
  as.data.frame

#Merge the dataframes
export_table <- bmd_table %>% left_join(summarized_results, by = "signal")  %>%
  dplyr::mutate(num_other_traits=ifelse(is.na(num_other_traits), 0, num_other_traits))

#Write to file
write.table(export_table, file = file.path(out_dir, paste0(locus_prefix, ".susie-coloc.tsv")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 8) Extract Signal Activities into a Table ====

#Read in trait file
trait_raw <- read.csv(trait_file, sep = "\t", header = TRUE)
#Extract list of traits
all_traits <- trait_raw$Trait

#Create a function set to extract signed PP4s of top signal
extract_top_pp4 <- function(non_bmd_trait, bmd_cs, results_list=results, susie_objs=susie_objects, bmd_name="BMD"){
  non_bmd_result <- results_list[[non_bmd_trait]]
  if (any(is.na(non_bmd_result))) {
    return(0)
  } else {
    cs_rows <- non_bmd_result$summary %>% filter(idx1 == bmd_cs)
    max_pp4_row <- cs_rows[which.max(PP.H4.abf),]
    sign <- extract_table_values(unlist(max_pp4_row[,c("idx1", "idx2")]), check_obj=susie_objs[[non_bmd_trait]])
    return(ifelse(sign == "+", max_pp4_row$PP.H4.abf, -max_pp4_row$PP.H4.abf))
  }
}
extract_signal_signed_pp4s <- function(bmd_cs, locus_trait_list=locus_traits[locus_traits != "BMD"], all_trait_list=all_traits){
  locus_traits_pp4s <- vapply(locus_trait_list, FUN = extract_top_pp4, FUN.VALUE = numeric(1), bmd_cs=bmd_cs)
  temp <- rep(0, length(all_trait_list))
  names(temp) <- all_trait_list
  temp[names(locus_traits_pp4s)] <- locus_traits_pp4s
  return(temp)
}
#Apply the functions to all the bmd signals
signal_signed_pp4s <- lapply(bmd_obj$sets$cs_index, FUN=extract_signal_signed_pp4s)

#Create the output dataframe for the activities
activity_summary <- bind_rows(signal_signed_pp4s) %>% as.data.frame()
row.names(activity_summary) <- paste(locus_prefix, bmd_obj$sets$cs_index, sep = ".")
#Write to file
write.table(activity_summary, file = file.path(out_dir, paste0(locus_prefix, ".signed_pp4s.tsv")), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
