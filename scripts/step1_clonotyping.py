#!/usr/bin/env python3
"""
Mouse BCR 10x Genomics Processing Pipeline
Executes Immcantation suite tools (AssignGenes, MakeDb, ParseDb, singlecell-filter)
and custom clonotyping on individual single-cell datasets.

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
import pandas as pd
from typing import List, Optional
from pathlib import Path

# --- 0. PRE-REQUISITE: CREATE LOG DIRECTORY ---
# The output folder must exist BEFORE the logger tries to save a file there.
LOG_DIR = Path("../output").resolve()
LOG_DIR.mkdir(parents=True, exist_ok=True)

# --- Configure Logging ---
LOG_DIR = Path("../output").resolve()

# Create formatters: one detailed for the file, one clean for the screen
file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
screen_formatter = logging.Formatter("%(message)s")

# File Handler (Detailed)
file_handler = logging.FileHandler(LOG_DIR / "pipeline_run.log")
file_handler.setFormatter(file_formatter)

# Stream Handler (Clean Screen)
screen_handler = logging.StreamHandler()
screen_handler.setFormatter(screen_formatter)

# Initialize Logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.addHandler(screen_handler)

# --- 1. CONFIGURATION ---
REQUESTED_CORES = 4
# Auto-detect hardware capacity
SYSTEM_CORES = os.cpu_count() or 1
N_CORES = min(REQUESTED_CORES, SYSTEM_CORES)

OUTPUT_DIR = Path("../output/clonotyping").resolve()
INPUT_BASE = Path("../input").resolve() 

REF_DIR = Path("/usr/local/share/germlines/imgt/mouse/vdj")
IGDATA = "/usr/local/share/igblast"
CLONOTYPE_R = "clonotyping.R"

# BCR Specific Settings
IGBLAST_LOCI = "ig"
HEAVY_LOCI = ["IGH"]
LIGHT_LOCI = ["IGK", "IGL"]
CLONOTYPE_THRESHOLD = "0.20"  
CLONOTYPE_NORM = "len"        

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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
    # Communicate resource allocation to the user
    logging.info(f"Resource Check: Requested {REQUESTED_CORES} cores.")
    if SYSTEM_CORES < REQUESTED_CORES:
        logging.info(f"Notice: System only has {SYSTEM_CORES} cores. Auto-adjusting allocation.")
    
    logging.info(f"Action: Initializing Step 1 using {N_CORES} core(s).")
    logging.info(f"Configuration Note: To manually adjust CPU usage, modify 'REQUESTED_CORES' on line 52.")
    logging.info("")
 
    # --- 2. DATASET DISCOVERY ---
    if not INPUT_BASE.exists():
        logging.error(f"Input directory does not exist: {INPUT_BASE}")
        sys.exit(1)

    # Identify all subdirectories containing the required 10x output
    datasets = [
        d.name for d in INPUT_BASE.iterdir() 
        if d.is_dir() and (d / "filtered_contig.fasta").exists()
    ]

    if not datasets:
        logging.error(f"No valid 10x datasets found in {INPUT_BASE}")
        sys.exit(1)

    logging.info(f"Found {len(datasets)} datasets to process: {', '.join(datasets)}")

    # --- 3. MAIN PIPELINE LOOP ---
    for lib_id in datasets:
        out_dir = OUTPUT_DIR / lib_id
        
        if out_dir.exists():
            logging.info(f"Skipping {lib_id} (Folder already exists)")
            continue
            
        out_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Processing: {lib_id}")

        fasta = INPUT_BASE / lib_id / "filtered_contig.fasta"
        annots = INPUT_BASE / lib_id / "filtered_contig_annotations.csv"
        
        try:
            # Step 1a: V(D)J Annotation via IgBLAST
            run_cmd([
                "AssignGenes.py", "igblast", "-s", fasta, "--organism", "mouse", "--loci", IGBLAST_LOCI,
                "-b", IGDATA, "--format", "blast", "--nproc", str(N_CORES),
                "--outname", lib_id, "--outdir", out_dir
            ], cwd=out_dir)

            # Step 1b: Database Construction
            fmt7_file = out_dir / f"{lib_id}_igblast.fmt7"
            run_cmd([
                "MakeDb.py", "igblast", "-i", fmt7_file, "-s", fasta,
                "--10x", annots, "-r", REF_DIR, "--extended", "--failed", 
                "--outname", lib_id, "--format", "airr", "--outdir", out_dir
            ], cwd=out_dir)

            # Step 1c: Format Standardization & QC
            db_pass = out_dir / f"{lib_id}_db-pass.tsv"
            if not db_pass.exists():
                raise FileNotFoundError(f"MakeDb failed to produce {db_pass.name}")

            df = pd.read_csv(db_pass, sep='\t', dtype=str)
            
            # Extract raw 10x barcodes prior to filtering for downstream linkage
            if 'bcs_used' not in df.columns:
                df['bcs_used'] = df['sequence_id'].str.split('_contig').str[0]
                df.to_csv(db_pass, sep='\t', index=False)

            # Export distribution of locus usage for quality control
            if 'locus' in df.columns:
                locus_counts = df['locus'].value_counts().reset_index()
                locus_counts.columns = ['locus', 'count']
                locus_counts.to_csv(out_dir / f"{lib_id}_dbpass_locus_counts.tsv", sep='\t', index=False)

            # --- STEP 2: SINGLE CELL FILTERING ---
            logging.info(f"{lib_id}: Filtering Single Cells...")
            
            heavy_tsv = out_dir / f"{lib_id}_heavy.tsv"
            light_tsv = out_dir / f"{lib_id}_light.tsv"

            # Isolate heavy and light chains into separate processing streams
            run_cmd(["ParseDb.py", "select", "-d", db_pass, "-f", "locus", "-u", *HEAVY_LOCI, "-o", heavy_tsv])
            run_cmd(["ParseDb.py", "select", "-d", db_pass, "-f", "locus", "-u", *LIGHT_LOCI, "-o", light_tsv])
            
            # Retain only productive sequences (no stop codons, in-frame junctions)
            run_cmd(["ParseDb.py", "split", "-d", heavy_tsv, "-f", "productive"])
            run_cmd(["ParseDb.py", "split", "-d", light_tsv, "-f", "productive"])

            heavy_prod = out_dir / f"{lib_id}_heavy_productive-T.tsv"
            light_prod = out_dir / f"{lib_id}_light_productive-T.tsv"
            
            # Re-associate paired chains into single-cell consensus records
            run_cmd([
                "singlecell-filter", "-d", f"{heavy_prod},{light_prod}",
                "-o", out_dir, "-f", "airr"
            ], cwd=out_dir)

            # --- STEP 3: CLONOTYPING ---
            logging.info(f"{lib_id}: Running Clonotyping...")
            
            h_sc_file = out_dir / f"{lib_id}_heavy_productive-T_sc-pass.tsv"
            l_sc_file = out_dir / f"{lib_id}_light_productive-T_sc-pass.tsv"
            
            if not h_sc_file.exists() or not l_sc_file.exists():
                raise FileNotFoundError("sc-pass files missing. singlecell-filter failed.")

            # Aggregate processed heavy and light chains for the clonotyping engine
            sc_combined_path = out_dir / f"{lib_id}_sc_combined.tsv"
            pd.concat([
                pd.read_csv(h_sc_file, sep='\t', dtype=str), 
                pd.read_csv(l_sc_file, sep='\t', dtype=str)
            ], ignore_index=True).to_csv(sc_combined_path, sep='\t', index=False)

            cloned_out = out_dir / f"{lib_id}_cloned.tsv"
            
            run_cmd(["Rscript", CLONOTYPE_R, sc_combined_path, CLONOTYPE_THRESHOLD, CLONOTYPE_NORM, cloned_out, str(N_CORES)])
        except Exception as e:
            logging.error(f"Error processing {lib_id}: {e}")
            continue

    logging.info("ALL DATASETS PROCESSED SUCCESSFULLY!")

if __name__ == "__main__":
    main()
