#!/bin/bash

# ==============================================================================
# SCRIPT:          Master BCR Repertoire Pipeline Execution
# PURPOSE:         Orchestrates the sequential execution of the BCR processing pipeline.
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

# --- STRICT ERROR HANDLING ---
set -euo pipefail

# --- CLI Splash Screen Colors ---
CYAN='\033[0;36m'
GREEN='\033[0;32m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Catch errors and print exactly where they happened
trap 'echo -e "\n${RED}❌ CRITICAL ERROR: Pipeline failed at line $LINENO. Halting execution.${NC}"' ERR

# --- Print the Splash Screen ---
echo -e "${CYAN}${BOLD}"
cat << 'EOF'
===============================================================================
             🧬 BCR REPERTOIRE ANALYSIS PIPELINE (v1.0) 🧬
===============================================================================
 Author:          Shahab Saghaei (Research Fellow & Data Manager)
 Affiliations:    Harvard Medical School | Broad Institute | Ragon Institute

 Publication:     "A glycan-based adjuvant expands the breadth and duration of 
                  protection of mRNA-based vaccines." (Nature Immunology)

 Container Image: Immcantation Framework (docker.io/immcantation/suite:4.5.0)
===============================================================================
EOF
echo -e "${NC}"
sleep 1 # Pauses for 1 second so the user can read the header

# --- Prepare Directories ---
mkdir -p ../output

echo -e "Initializing pipeline execution sequence...\n"

# --- Step 0: Reconstruct Input Structure from GEO ---
echo -e "---> [Step 0/3] Reconstructing input structure from GEO data..."

# Define the source and target directories
GEO_DIR="../GSE315320_RAW"
TARGET_DIR="../input"

# Ensure the target input directory exists
mkdir -p "$TARGET_DIR"

# Loop through the GEO folder to find 'filtered_contig.fasta.gz' files
# These identify our primary BCR datasets
for fasta_path in "$GEO_DIR"/*_filtered_contig.fasta.gz; do
    if [ -f "$fasta_path" ]; then
        # 1. Extract the filename
        filename=$(basename "$fasta_path")
        
        # 2. Determine the Library ID (e.g., WL-156-BCR) by removing 
        # the GSM prefix and the file suffix.
        lib_id=$(echo "$filename" | sed 's/GSM[0-9]*_//; s/_filtered_contig.fasta.gz//')
        
        echo "      Processing: $lib_id"
        
        # 3. Create the nested directory expected by the pipeline
        lib_dir="$TARGET_DIR/$lib_id"
        mkdir -p "$lib_dir"
        
        # 4. Find and copy the matching annotation file
        gsm_prefix=$(echo "$filename" | cut -d'_' -f1)
        annot_source="$GEO_DIR/${gsm_prefix}_${lib_id}_filtered_contig_annotations.csv.gz"
        
        # 5. Copy and rename files into the nested structure
        cp "$fasta_path" "$lib_dir/filtered_contig.fasta.gz"
        cp "$annot_source" "$lib_dir/filtered_contig_annotations.csv.gz"
        
        # 6. Decompress the files in the new location (overwriting if exists)
        gunzip -f "$lib_dir/filtered_contig.fasta.gz"
        gunzip -f "$lib_dir/filtered_contig_annotations.csv.gz"
    fi
done

echo -e "     ✅ Input structure successfully reconstructed in '$TARGET_DIR'.\n"

# --- Execute Pipeline Steps ---
echo -e "---> [Step 1/3] Running Clonotyping Pipeline..."
python3 step1_clonotyping.py

echo -e "\n---> [Step 2/3] Formatting CDR3 Lengths for Prism..."
python3 step2_cdr3_length.py

echo -e "\n---> [Step 3/3] Running Diversity & V-Gene Analysis..."
python3 step3_diversity.py

echo -e "\n==========================================================================="
echo -e "${GREEN}🎉 PIPELINE COMPLETE! All results are saved in the 'output/' folder.${NC}"
echo "==========================================================================="
