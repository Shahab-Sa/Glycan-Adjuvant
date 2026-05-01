#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT:          BCR Hierarchical Clonotyping
# PURPOSE:         Executes scoper's hierarchicalClones to cluster single-cell B cell 
#                  repertoires into clonal lineages based on junctional similarity.
#
# AUTHOR:          Shahab Saghaei (Research Fellow & Data Manager)
# AFFILIATIONS:    Harvard Medical School | Broad Institute | Ragon Institute
# CONTACT:         https://www.linkedin.com/in/shahabsa/
#
# PUBLICATION:     "A glycan-based adjuvant expands the breadth and duration of 
#                  protection of mRNA-based vaccines." (Nature Immunology)
# CONTAINER IMAGE: Immcantation Framework (docker.io/immcantation/suite:4.5.0)
# LICENSE:         BSD 3-Clause
# ==============================================================================

suppressPackageStartupMessages({
  library(scoper)
  library(dplyr)
})

# --- 1. SETUP & ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[1]
threshold   <- as.numeric(args[2])
norm_method <- args[3]
output_file <- args[4]
n_cores     <- as.numeric(args[5])

db <- read.delim(input_file, sep = "\t", header = TRUE)

# --- 2. CLONOTYPING ---
# Perform hierarchical clustering on the heavy chain CDR3 regions
results <- hierarchicalClones(
    db,
    cell_id = "cell_id",
    locus = "locus",
    threshold = threshold,
    normalize = norm_method,
    only_heavy = TRUE,
    split_light = TRUE,
    summarize_clones = FALSE,
    nproc = n_cores # Use the variable passed from Python
)

# --- 3. REPORTING & EXPORT ---
cat("Clonotyping complete:", output_file, "\n")

cat("\n          START> hierarchicalClones\n")
cat(paste("      THRESHOLD>", threshold, "\n"))
cat(paste("    NORM_METHOD>", norm_method, "\n")) 
cat(paste("        THREADS>", n_cores, "\n"))
cat(paste("          INPUT>", input_file, "\n"))
cat(paste("         OUTPUT>", output_file, "\n"))
cat(paste("        RECORDS>", nrow(db), "\n"))
cat(paste("     NUM_CLONES>", length(unique(results$clone_id)), "\n"))
cat("            END> hierarchicalClones\n\n")

write.table(results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
