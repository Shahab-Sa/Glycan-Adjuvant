#!/usr/bin/env python3
"""
BCR Clonotype Prism Formatter
Extracts CDR3 lengths from clonotyped BCR repertoires and formats 
the output into TSV files compatible with GraphPad Prism.

Author:          Shahab Saghaei (Research Fellow & Data Manager)
Affiliations:    Harvard Medical School | Broad Institute | Ragon Institute
Contact:         https://www.linkedin.com/in/shahabsa/

Publication:     "A glycan-based adjuvant expands the breadth and duration of 
                 protection of mRNA-based vaccines." (Nature Immunology)
Container Image: Immcantation Framework (docker.io/immcantation/suite:4.5.0)
License:         BSD 3-Clause
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path

# --- 0. LOGGING SETUP (Clean Screen, Detailed File) ---
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
OUTPUT_DIR = Path("../output/cdr3_lengths").resolve()
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Map human-readable condition labels to specific processed outputs
FILES = {
    "WA1_CMA": "../output/clonotyping/WL-156-BCR/WL-156-BCR_cloned.tsv",
    "WA1_C": "../output/clonotyping/WL-159-BCR/WL-159-BCR_cloned.tsv",
    "XBB_C": "../output/clonotyping/WL-164-BCR/WL-164-BCR_cloned.tsv",
    "XBB_CMA": "../output/clonotyping/WL-167-BCR/WL-167-BCR_cloned.tsv"
}

def process_datasets() -> None:
    # --- 2. LOAD & PROCESS DATA ---
    all_data = []
    
    for condition, path in FILES.items():
        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
        except FileNotFoundError:
            logging.warning(f"Could not find file {path}. Skipping.")
            continue
        
        required_cols = ["junction_aa", "locus", "junction", "cdr3"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"File {path} missing required columns: {', '.join(missing_cols)}")
        
        # Standardize locus identifiers into biological chain types
        df["chain"] = df["locus"].apply(lambda x: "Heavy" if x == "IGH" else ("Light" if x in ["IGK", "IGL"] else None))
        df = df.dropna(subset=["chain"])
        
        # Sanity check: Ensure junction sequences cleanly translate to CDR3 by 
        # verifying the length difference matches the two conserved boundary residues (Cys/Trp).
        junction_len = df["junction"].apply(len)
        cdr3_len_nt = df["cdr3"].apply(len)
        diff = junction_len - cdr3_len_nt
        
        if not (diff == 6).all():
            non_six_values = diff[diff != 6].unique()
            raise ValueError(f"In file {path}, junction_len - cdr3_len_nt != 6. Offending values: {non_six_values}")
        
        # Calculate biological CDR3 length in amino acids
        df["cdr3_length"] = df["junction_aa"].apply(len) - 2
        df["condition"] = condition
        
        all_data.append(df[["condition", "chain", "cdr3_length"]])
        
    if not all_data:
        logging.error("No data was successfully loaded. Exiting.")
        return

    cdr3_df = pd.concat(all_data, ignore_index=True)
    
    # --- 3. FORMAT & EXPORT FOR PRISM ---
    for chain_type in ["Heavy", "Light"]:
        df_chain = cdr3_df[cdr3_df["chain"] == chain_type]
        
        # Group lengths by condition and sort descending for visualization hierarchy
        grouped = df_chain.groupby('condition')['cdr3_length'].apply(lambda x: sorted(x, reverse=True)).to_dict()
        conditions = sorted(FILES.keys(), reverse=True)
        
        max_len = max((len(grouped[cond]) for cond in conditions if cond in grouped), default=0)
        if max_len == 0:
            continue
            
        # Pad columns with NaNs to satisfy Prism's requirement for rectangular data structures
        df_for_prism = pd.DataFrame({
            f"{cond}_{len(grouped[cond])}": grouped[cond] + [np.nan] * (max_len - len(grouped[cond]))
            for cond in conditions if cond in grouped
        })
        
        output_file = OUTPUT_DIR / f"cdr3_lengths_{chain_type.lower()}_chain.tsv"
        df_for_prism.to_csv(output_file, sep="\t", index=False)
        logging.info(f"Successfully saved Prism template: {output_file.name}")

if __name__ == "__main__":
    process_datasets()
