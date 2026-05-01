#!/usr/bin/env python3
"""
BCR Diversity Analysis Wrapper
Iterates over paired BCR datasets and executes the R-based diversity pipeline.

Author:          Shahab Saghaei (Research Fellow & Data Manager)
Affiliations:    Harvard Medical School | Broad Institute | Ragon Institute
Contact:         https://www.linkedin.com/in/shahabsa/

Publication:     "A glycan-based adjuvant expands the breadth and duration of 
                 protection of mRNA-based vaccines." (Nature Immunology)
Container Image: Immcantation Framework (docker.io/immcantation/suite:4.5.0)
License:         BSD 3-Clause
"""

import os
import sys
import subprocess
import logging
from typing import List, Optional
from pathlib import Path

# --- 0. PRE-REQUISITE: CREATE LOG DIRECTORY ---
# The output folder must exist BEFORE the logger tries to save a file there.
LOG_DIR = Path("../output").resolve()
LOG_DIR.mkdir(parents=True, exist_ok=True)

# --- Configure Logging ---
LOG_DIR = Path("../output").resolve()

# Formatters: Detailed for the file, Clean for the terminal
file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
screen_formatter = logging.Formatter("%(message)s")

# File Handler
file_handler = logging.FileHandler(LOG_DIR / "pipeline_run.log")
file_handler.setFormatter(file_formatter)

# Screen Handler
screen_handler = logging.StreamHandler()
screen_handler.setFormatter(screen_formatter)

# Logger Initialization
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.addHandler(screen_handler)

# --- 1. CONFIGURATION ---
R_SCRIPT = "diversity_analysis.R"
MAIN_OUT_DIR = Path("../output").resolve() 

# Define analytical cohorts for pairwise comparisons
PAIRS = [
    {
        "prefix": "WA1",
        "file_1": "../output/clonotyping/WL-159-BCR/WL-159-BCR_cloned.tsv",
        "file_2": "../output/clonotyping/WL-156-BCR/WL-156-BCR_cloned.tsv",
        "name_1": "C (WA1)",
        "name_2": "CMA (WA1)"
    },
    {
        "prefix": "XBB",
        "file_1": "../output/clonotyping/WL-164-BCR/WL-164-BCR_cloned.tsv",
        "file_2": "../output/clonotyping/WL-167-BCR/WL-167-BCR_cloned.tsv",
        "name_1": "XBB.1.5 (XBB.1.5)",
        "name_2": "CMA (XBB.1.5)"
    }
]

def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Executes a terminal command, capturing its output into the Python logger."""
    cmd_str_list = [str(arg) for arg in cmd]
    
    # Open the process and intercept its standard output and errors
    process = subprocess.Popen(
        cmd_str_list, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT, 
        text=True, 
        cwd=cwd
    )
    
    # Read the output line-by-line in real-time and pass it to the logger
    if process.stdout:
        for line in process.stdout:
            # .strip() removes the extra hidden newline characters
            logging.info(line.strip())
            
    # Wait for the command to finish and check if it failed
    process.wait()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)

def main() -> None:
    logging.info("Initializing Diversity & V-Gene analysis.")
    logging.info("Performance Note: Running in single-threaded mode due to alakazam library constraints.")    
    
    # --- 2. EXECUTE PIPELINE ---
    if not Path(R_SCRIPT).exists():
        logging.error(f"R script '{R_SCRIPT}' not found in the current directory.")
        sys.exit(1)

    logging.info(f"Found {len(PAIRS)} dataset pairs to process.")

    for pair in PAIRS:
        logging.info(f"Processing Cohort: {pair['prefix'].upper()}")
        
        # Verify inputs exist prior to initiating the R environment
        if not Path(pair["file_1"]).exists() or not Path(pair["file_2"]).exists():
            logging.warning(f"Missing input files for {pair['prefix']}. Skipping pair.")
            continue
        
        cmd = [
            "Rscript", R_SCRIPT,
            pair["file_1"],
            pair["file_2"],
            pair["name_1"],
            pair["name_2"],
            str(MAIN_OUT_DIR),
            pair["prefix"]
        ]
        
        try:
            run_cmd(cmd)
        except subprocess.CalledProcessError:
            logging.error(f"Error processing {pair['prefix']}: R script failed.")
            sys.exit(1) # This tells Bash that the script actually failed


    logging.info("ALL PAIRS PROCESSED SUCCESSFULLY!")

if __name__ == "__main__":
    main()
