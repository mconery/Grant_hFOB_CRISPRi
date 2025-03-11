################################################################################

#finemap_SuSiE.R

#This script is designed to execute SuSiE fine-mapping given a fixed region 
#when called by finemap_SuSiE.py after temp file generation. The code saves its 
#output as a .rds file that can later be used for plotting or for downstream 
#analysis.

################################################################################

# 0) Load Needed Libraries ====

library(coloc)
library(stringr)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following: 
  # 1) summary stats file filtered for locus snps (Assumes trait first part of name; filtered to mapping snps)
  # 2) ld ref matrix (assumes snp names present in file; filtered to mapping snps)
  # 3) sample size (maps binary traits as continous)
  # 4) out file name
  # 5) random seed [OPTIONAL; default is 5]
  # 6) confidence level [OPTIONAL; default is 5]

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 6) {
  #Print confirmation statement
  print("Necessary and Optional Parameters Input")
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  ld_ref_matrix <- args[2]
  sample_size_file <- args[3]
  out_file <- args[4]
  random_seed <- as.integer(args[5])
  confidence <- as.numeric(args[6])
} else if(length(args) == 4) {
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  ld_ref_matrix <- args[2]
  sample_size_file <- args[3]
  out_file <- args[4]
  random_seed <- 5
  confidence <- 0.95
} else {
  print("ERROR: Invalid inputs given.")
  print("Usage: %> Rscript finemap_SuSiE.R summary_stats_file ld_matrix_file sample_size_file out_file [random_seed] [confidence_level]");
  quit(save="no");
}

# 2) Prep for executing SuSiE ====

#Get trait name and file prefix
trait_name = str_split(basename(summary_stats_file), "[.]")[[1]][1] 
file_prefix = str_remove(basename(out_file), ".susie.rds") 

#Read in summary_statistics and ld matrix
summary_stats_filt <- read.table(summary_stats_file, header = TRUE)
ld_matrix_filt <- read.table(ld_ref_matrix, header = TRUE, row.names = 1)

#Reset colnames of ld matrix to rownames due to read in issue
colnames(ld_matrix_filt) <- rownames(ld_matrix_filt)

### Identify trait type and process accordingly ###
#Read in sample size file
sample_size_raw <- read.csv(sample_size_file, header = TRUE, sep = "\t", row.names = 1)
#Check if trait has any cases and controls entries
#Then add relevant info to summary statistics
if(is.na(sample_size_raw[trait_name,"Cases"]) && is.na(sample_size_raw[trait_name,"Controls"])){
  trait_type = "quant"
  sample_size = as.numeric(sample_size_raw[trait_name,"Sample.Size"])
  summary_stats_filt <- cbind.data.frame(summary_stats_filt, 
                                         num_samples= sample_size)
} else if(is.na(sample_size_raw[trait_name,"Cases"]) || is.na(sample_size_raw[trait_name,"Controls"])) {
  print(paste0("ERROR: Indeterminate Trait Type. Check Sample size file for ", trait_name, "."))
  quit(save="no");
} else {
  trait_type = "cc"
  cases = as.numeric(sample_size_raw[trait_name,"Cases"])
  controls = as.numeric(sample_size_raw[trait_name,"Controls"])
  sample_size = as.numeric(sample_size_raw[trait_name,"Sample.Size"])
  summary_stats_filt <- cbind.data.frame(summary_stats_filt, 
                                         num_samples= sample_size,
                                         prop_cases=(cases/(cases + controls)))
}

#Check allele frequencies. Increase any 0% up to to a minimum value
#Get minimum values of the allele frequency columns
min_FRQ <- min(min(summary_stats_filt[which(summary_stats_filt$MAF > 0), "MAF"]), 1/(2*max(summary_stats_filt$num_samples)))
#Over write values
summary_stats_filt$MAF <- ifelse(summary_stats_filt$MAF < min_FRQ, min_FRQ, summary_stats_filt$MAF)
#Also overwrite max values
summary_stats_filt$MAF <- ifelse(summary_stats_filt$MAF > 1-min_FRQ, 1-min_FRQ, summary_stats_filt$MAF)

#Set variables for SuSiE
n <- max(summary_stats_filt$num_samples)
R <- ld_matrix_filt %>% as.matrix()
beta <- summary_stats_filt$BETA
varbeta <- (summary_stats_filt$SE)^2
MAF <- summary_stats_filt$MAF

#Check for cc and make lists to create coloc datasets
print("Creating coloc datasets")
if(trait_type == "cc") {
  coloc_data <- list(snp=summary_stats_filt$SNP, position=summary_stats_filt$BP,
                     type="cc", N=n, MAF=MAF, s=summary_stats_filt$prop_cases, LD=R, beta=beta, varbeta=varbeta) 
} else {
  coloc_data <- list(snp=summary_stats_filt$SNP, position=summary_stats_filt$BP,
                     type="quant", N=n, MAF=MAF, LD=R, beta=beta, varbeta=varbeta)
}
#Check dataset
if (is.null(check_dataset(coloc_data))) {
  print("Dataset check passed.")
} else{
  print("Dataset check failed.")
  check_dataset(coloc_data)
}

# 3) Set error checking function ====

check <- function(expression){
  
  withCallingHandlers(expression,
                      
                      warning = function(w){
                      },
                      error = function(e){
                        message(paste0("ERROR: ", file_prefix, " is unmappable."))
                        #Write an empty output file and quit
                        fitted_rss1 = NULL
                        saveRDS(fitted_rss1, file = out_file)
                        quit(save="no");
                      },
                      finally = {
                      })
}

# 4) Execute SuSiE ====

#Execute SuSiE
set.seed(random_seed) #Nothing up until here has been random
fitted_rss1 = check(runsusie(coloc_data, n=n, min_abs_corr=0, coverage=confidence, estimate_residual_variance = FALSE, L = 5, maxit=10000, repeat_until_convergence=FALSE))

#Verify convergence
if (fitted_rss1$converged == FALSE) {
  print(paste0("ERROR: ", file_prefix, " did not converge in < 10,000 iterations. Outputting empty rds file."))
  #Write an empty output file and quit
  fitted_rss1 = NULL
  saveRDS(fitted_rss1, file = out_file)
  quit(save="no");
}

## Calculate residuals ##
#Get credible sets
cred_sets=names(fitted_rss1$sets$cs)
cred_sets=vapply(cred_sets, str_replace, character(1), pattern="L", replacement="")
cred_sets=vapply(cred_sets, as.numeric, 0)
cred_sets=sort(cred_sets)

#Calculate residuals
for (i in cred_sets) {
  cred_set <- names(fitted_rss1$sets$cs[[paste0("L", i)]])
  #Calculate absolute residuals
  a = fitted_rss1$mu[-i,]
  b = fitted_rss1$alpha[-i,]
  z = beta / sqrt(varbeta)
  if (is.null(ncol(a)) && is.null(ncol(b))) {
    prediction = colSums(t(a) * t(b))
  } else{
    prediction = colSums(a * b)
  }
  test = matrix(sqrt(n - 1) * z) - ((n - 1) * R %*% prediction)
  test2 = (test)/sqrt(fitted_rss1$sigma2 * (n-1) * diag(R))
  neglogp_absolute = (-log(2)-pt(-abs(test2), df = n - 1, lower.tail = TRUE, log.p = TRUE))/log(10)
  #Check if i == 1 to process info to data frames correctly
  if(i == 1){
    neglogp_absolute_df = as.data.frame(neglogp_absolute)
  } else{
    neglogp_absolute_df = cbind.data.frame(neglogp_absolute_df, neglogp_absolute)
  }
}

#Check if the residual dfs were actually populated (They won't be if no credible sets were identifiable)
if (exists("neglogp_absolute_df")) {
  #Set column names
  colnames(neglogp_absolute_df) <- cred_sets
} else {
  print(paste0("WARNING: No credible sets detected for ", file_prefix))
  neglogp_absolute_df <- NULL
}

#Append residuals to SuSiE object
fitted_rss1$neglogp_absolute <- neglogp_absolute_df

#Append chr, bp, and neglogp_gwas to SuSiE object
fitted_rss1$chr <- median(summary_stats_filt$CHR)
fitted_rss1$bp <- summary_stats_filt$BP
fitted_rss1$rsid <- summary_stats_filt$RSID
fitted_rss1$neglogp_gwas <- -log10(summary_stats_filt$P)

#Save susie_object
saveRDS(fitted_rss1, file = out_file)
