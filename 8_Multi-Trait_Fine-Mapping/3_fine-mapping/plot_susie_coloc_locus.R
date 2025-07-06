################################################################################
# Plot SuSiE–coloc results for a single locus
#
# USAGE:
#   Rscript plot_susie_coloc_locus.R locus_prefix \
#           [--rds_dir ./] [--out_dir ./] [--bmd_trait BMD] \
#           [--p1 1e-4] [--p2 1e-4] [--p12 1e-5] [--pp4_cut 0.5]
#
# EXAMPLE:
#   Rscript plot_susie_coloc_locus.R chr6.12345678.12355678 \
#           --rds_dir ./susie_rds/ --out_dir ./plots/ --bmd_trait eBMD
#
# The script creates a stacked set of scatter plots (PDF & PNG) where:
#   •  x–axis  = genomic position (bp)
#   •  y–axis  = –log10(P-value) from GWAS
#   •  colour  = signal membership / colocalisation status
#
# Colour scheme
#   light‐grey  : variant not in any credible set
#   dark‐grey   : variant in a non-BMD credible set that does *not* colocalise
#   rainbow(n)  : variants in BMD signals and the non-BMD signals that
#                 colocalise with them (one unique colour per BMD signal)
################################################################################

suppressPackageStartupMessages({
  library(coloc)       # coloc.susie
  library(ggplot2)
  library(cowplot)     # plot_grid
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tools)
})

################################################################################
# 1. COMMAND-LINE ARGUMENTS
################################################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_susie_coloc_locus.R locus_prefix [--options]")
}
## defaults
opt <- list(
  locus_prefix = args[1],
  rds_dir   = "./",
  out_dir   = "./",
  bmd_trait = "BMD",
  p1        = 1e-4,
  p2        = 1e-4,
  p12       = 1e-5,
  pp4_cut   = 0.50
)

if (length(args) > 1) {
  flag_pos <- which(grepl("^--", args))
  for (f in flag_pos) {
    key   <- gsub("^--", "", args[f])
    value <- args[f + 1]
    if (key %in% names(opt)) {
      if (key %in% c("p1","p2","p12","pp4_cut")) {
        opt[[key]] <- as.numeric(value)
      } else {
        opt[[key]] <- value
      }
    }
  }
}

locus_prefix <- opt$locus_prefix

################################################################################
# 2. LOAD SuSiE OBJECTS
################################################################################
rds_pattern <- paste0("*.", locus_prefix, ".*\\.susie\\.rds$")
rds_files   <- list.files(opt$rds_dir, pattern = rds_pattern, full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No *.susie.rds files found for locus ", locus_prefix)
}

susie_objects <- lapply(rds_files, readRDS)
names(susie_objects) <- str_remove(
  str_remove(basename(rds_files), "\\.susie\\.rds$"),
  paste0("\\.", locus_prefix)
)

## split BMD vs other traits
bmd_name <- names(susie_objects)[grepl(opt$bmd_trait, names(susie_objects), ignore.case = TRUE)]
if (length(bmd_name) != 1) {
  stop("Could not uniquely identify BMD trait – use --bmd_trait to set it.")
}
bmd_obj      <- susie_objects[[bmd_name]]
other_traits <- setdiff(names(susie_objects), bmd_name)

################################################################################
# 3. MAP COLOCALISING SIGNALS (re-run coloc quickly – lightweight)
################################################################################
message("Running pair-wise coloc to map signal pairs ...")
bmd_cs   <- bmd_obj$sets$cs_index
n_signals <- length(bmd_cs)
signal_colours <- setNames(rainbow(n_signals), bmd_cs)   # colour per BMD signal

## list to hold mappings: one row per colocating pair
pair_map <- list()

for (tr in other_traits) {
  coloc_res <- coloc.susie(
    dataset1 = bmd_obj,
    dataset2 = susie_objects[[tr]],
    p1 = opt$p1, p2 = opt$p2, p12 = opt$p12
  )
  
  summ <- coloc_res$summary %>%
    filter(PP.H4.abf > opt$pp4_cut)
  if (nrow(summ) > 0) {
    pair_map[[tr]] <- summ %>% select(idx1, idx2) %>% 
      mutate(trait = tr)
  }
}

pair_df <- bind_rows(pair_map)

# Only keep traits that have at least one signal colocalizing with BMD
traits_with_coloc <- unique(pair_df$trait)
other_traits_filtered <- intersect(other_traits, traits_with_coloc)


################################################################################
# 4. FUNCTION: ASSIGN COLOURS PER VARIANT
################################################################################
get_variant_df <- function(obj, trait_name) {
  df <- data.frame(
    snp  = names(obj$pip),
    pos  = obj$bp,
    y    = obj$neglogp_gwas,
    cs   = NA_character_,                 # credible set ID (e.g. "L3")
    col  = "lightgrey"                    # default colour
  )
  
  ## variants in any CS
  if (!is.null(obj$sets$cs) && length(obj$sets$cs) > 0) {
    for (cs_name in names(obj$sets$cs)) {
      v <- obj$sets$cs[[cs_name]]
      df$cs[v] <- cs_name
    }
  }
  
  ## colour rules
  for (i in seq_len(nrow(df))) {
    cs <- df$cs[i]
    if (is.na(cs)) next                       # stays light-grey
    if (trait_name == bmd_name) {
      ## BMD: map cs index → colour
      idx <- as.integer(str_remove(cs, "L"))
      df$col[i] <- signal_colours[as.character(idx)]
    } else {
      ## Non-BMD: does this CS colocalise?
      idx2_match <- pair_df %>% 
        filter(trait == trait_name,
               idx2 == as.integer(str_remove(cs,"L"))) %>%
        pull(idx1)
      if (length(idx2_match) > 0) {
        df$col[i] <- signal_colours[as.character(idx2_match[1])]
      } else {
        df$col[i] <- "darkgrey"
      }
    }
  }
  df$trait <- trait_name
  return(df)
}

################################################################################
# 5. BUILD DATA FOR ALL TRAITS
################################################################################
plot_dfs <- lapply(c(bmd_name, other_traits_filtered), function(tr) {
  get_variant_df(susie_objects[[tr]], tr)
})
names(plot_dfs) <- c(bmd_name, other_traits_filtered)

#Create function to clean up trait names
correct_trait_names <- function(trait_names){
  trait_names %>% str_replace_all(pattern = "_", replacement = " ") %>%
    str_replace_all(pattern = "[.]", replacement = "-") %>%
    toTitleCase %>%
    str_replace(pattern = "Alp", "ALP") %>%
    str_replace(pattern = "Alt", "ALT") %>%
    str_replace(pattern = "BF", "Bone Fracture") %>%
    str_replace(pattern = "Bmi", "BMI") %>%
    str_replace(pattern = "Ldl", "LDL") %>%
    str_replace(pattern = "Hdl", "HDL") %>%
    str_replace(pattern = "Egfr Creat", "eGFR Creatinine") %>%
    str_replace(pattern = "Ggt", "GGT") %>%
    str_replace(pattern = "Hba1c", "HbA1c") %>%
    str_replace(pattern = "Vitamin d", "Vitamin D") %>%
    str_replace(pattern = "Smoking Ever Never", "Smoking Ever/Never") %>%
    str_replace(pattern = "Lbs", "(lbs)") %>%
    str_replace(pattern = "Bone Mineral Density", "BMD (Pan-UKBB)") %>%
    str_replace(pattern = "Whole Body", replacement = "Whole-Body") %>%
    return
}
#Create cleaned trait names for plot
clean_trait_names <- correct_trait_names(names(plot_dfs))


################################################################################
# 6. GENERATE STACKED PLOTS
################################################################################
message("Generating plots ...")

make_panel <- function(df, trait_label) {
  ggplot(df, aes(x = pos, y = y, colour = col)) +
    geom_point(size = 1.3) +
    scale_colour_identity() +
    labs(title = trait_label, x = NULL, y = expression(-log[10](P))) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0, face = "bold"),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      panel.grid.minor = element_blank()
    )
}

panels <- mapply(plot_dfs, FUN = make_panel, 
                 trait_label = clean_trait_names,
                 SIMPLIFY = FALSE)

## bottom panel keeps x-axis labels
panels[[length(panels)]] <- panels[[length(panels)]] +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_line())

plot_stack <- plot_grid(plotlist = panels, ncol = 1, align = "v")

################################################################################
# 7. SAVE OUTPUT
################################################################################
out_base <- file.path(opt$out_dir, paste0(locus_prefix, ".susie_coloc_plot"))
ggsave(paste0(out_base, ".pdf"), plot_stack, width = 10, height = 2.5*length(panels))
ggsave(paste0(out_base, ".png"), plot_stack, width = 10, height = 2.5*length(panels), dpi = 300)

message("Plot written to: ", out_base, ".pdf / .png")
