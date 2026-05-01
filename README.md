# BCR Repertoire Analysis Pipeline 🧬

This repository contains the complete, reproducible computational pipeline used to process and analyze murine B cell receptor (BCR) single-cell sequencing data for the publication: 

> **"A glycan-based adjuvant expands the breadth and duration of protection of mRNA-based vaccines."**
> *Nature Immunology*

**Author:** Shahab Saghaei (Research Fellow & Data Manager)  
**Affiliations:** Harvard Medical School | Broad Institute | Ragon Institute  

---

## 📌 Overview

This pipeline provides a turnkey, containerized workflow for taking raw **10x Genomics V(D)J** outputs (specifically BCR filtered contigs and annotations) and generating publication-ready metrics. 

For long-term reproducibility, the pipeline is designed to directly ingest the raw supplemental files from **NCBI GEO Accession: GSE315320**. The script automatically handles the unpacking, renaming, and organization of these flat files into a proper hierarchy independent of external database connectivity.

The automated workflow performs the following sequence:
1. **Step 0:** Automatically scans the raw NCBI GEO download folder, extracts library IDs, reconstructs the necessary directory structure, and unpacks the compressed datasets.
2. **Step 1:** Performs IgBLAST-based V(D)J gene assignment, formats data to the AIRR-C (Adaptive Immune Receptor Repertoire Community) standard, enforces heavy/light chain pairing consensus, and executes length-normalized hierarchical clonotyping incorporating light-chain splitting (`scoper`).
3. **Step 2:** Calculates biological CDR3 amino acid lengths by validating junction sequences and excluding conserved boundary residues, then formats the output matrix for GraphPad Prism integration.
4. **Step 3:** Computes bootstrapped clonal abundance, evaluates rarefied Hill diversity spectra ($q=0$ to $q=4$, including Shannon Entropy and Inverse Simpson indices) (`alakazam`), and calculates Log2 fold-change (Log2FC) for V-Gene usage across cohorts. 

---
## 🛠 Prerequisites

The core BCR pipeline is designed with a "Container-First" philosophy to ensure 100% reproducibility.

**Option 1: Containerized (Recommended)**
It eliminates local dependency conflicts by using the **Immcantation Framework** environment. You do not need to install R, Python, or IMGT databases locally for these steps.
* **Requirements:** Docker or Podman installed.
* **Target Image:** `docker.io/immcantation/suite:4.5.0`
* **Note:** This container includes the internal IMGT/V-QUEST reference directory for consistent gene assignment.

**Option 2: Local Installation (Manual)**
If you prefer to run the BCR pipeline outside of a container:
* **Requirements:** Python 3.8+ & R 4.0+
* **R Packages:** `alakazam`, `scoper`, `dplyr`, `airr`, `ggplot2`, `tidyr`.
* **Python Packages:** `pandas`, `numpy`.
* **IMGT DB:** Ensure your local IgBLAST/IMGT database is correctly configured.


---

## 🚀 Quick Start (One-Click Reproduction)

### 1. Obtain the Repository
You can obtain the code and data in one of two ways:

**Option A: Command Line** Run these commands in your terminal:  
`git clone https://github.com/Shahab-Sa/Glycan-Adjuvant.git`  
`cd Glycan-Adjuvant`  

**Option B: Manual Download** Download the ZIP file from the green **"Code"** button at the top of this page, extract it, and open your terminal **inside** the extracted folder.

### 2. Download the Raw Data
The required datasets are hosted on **NCBI GEO Accession: [GSE315320](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315320)**. 

**Option A: Direct Download (Quick Setup)**
* [Click here to directly download the supplementary RAW archive (`GSE315320_RAW.tar`)](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE315320&format=file).

**Option B: Manual Navigation**
* Navigate to the [GSE315320 accession page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315320).
* Scroll down to the "Supplementary file" section at the bottom of the page.
* Locate the `GSE315320_RAW.tar` file and click on the **(http)** link in the Download column.

Extract the contents into a new folder named `GSE315320_RAW` at the root of your cloned repository. Your structure must look like this before proceeding:
```text
Glycan-Adjuvant/
├── GSE315320_RAW/
│   ├── GSM9425132_WL-154-155-GEX-ADT_barcodes.tsv.gz
│   ├── ... (other files)
│   └── GSM9425143_WL-167-BCR_filtered_contig.fasta.gz
└── scripts/
```

### 3. Execute the BCR Pipeline
**⚠️ Terminal Path Requirement:** Before running the pipeline, your terminal must be located at the root of the `Glycan-Adjuvant` folder.

Run the following command using Docker to initiate the pipeline:

```bash
docker run --rm -it --platform linux/amd64 -v "$PWD":/data -w /data/scripts docker.io/immcantation/suite:4.5.0 bash run_pipeline.sh
```

---

## 📂 Repository Structure

```text
├── README.md
├── GSE315320_RAW/              # Extracted raw supplemental files from NCBI GEO
│   ├── GSM9425132_WL-154-155-GEX-ADT_barcodes.tsv.gz
│   ├── GSM9425134_WL-156-BCR_filtered_contig_annotations.csv.gz
│   ├── GSM9425134_WL-156-BCR_filtered_contig.fasta.gz
│   ├── ...
│   └── GSM9425143_WL-167-BCR_filtered_contig.fasta.gz
└── scripts/                    # Core pipeline architecture
    ├── run_pipeline.sh         # Master execution script (rebuilds /input automatically)
    ├── step1_clonotyping.py    # Python wrapper for Step 1
    ├── clonotyping.R           # R engine for hierarchical cloning
    ├── step2_cdr3_length.py    # Python engine for CDR3 Prism formatting
    ├── step3_diversity.py      # Python wrapper for Step 3
    └── diversity_analysis.R    # R engine for diversity & V-Gene fold-change
```
*(Note: An `input/` folder will be dynamically generated by the automated pipeline during Step 0 based on the contents of `GSE315320_RAW`)*

---

## 📊 Expected Output

Upon successful completion of the automated pipeline, a new `output/` directory will be generated in the repository root. This folder will contain:
* **`pipeline_run.log`**: Detailed execution logs.
* **`clonotyping/`**: Filtered heavy/light chain databases and finalized clonotyped TSVs.
* **`cdr3_lengths/`**: Formatted TSV templates ready for Prism visualization.
* **`bcr_summary/`**: Clone size frequency tables.
* **`diversity_analysis/`**: Rarefied Hill diversity statistics, Shannon entropy calculations, and `alakazam` visual plots.
* **`v_gene_usage/`**: Heatmaps and Fold-Change rankings for V-Gene utilization.

---

## 📚 Citations & Acknowledgements

This computational pipeline was fundamentally built upon the **Immcantation Framework**. We gratefully acknowledge and cite the authors of the suite and the specific packages utilized in our workflow:

* **Immcantation / Alakazam (Diversity, Abundance & Gene Usage):**
  > Gupta, N. T., Vander Heiden, J. A., Uduman, M., Gadala-Maria, D., Yaari, G., & Kleinstein, S. H. (2015). Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics*, 31(20), 3356–3358. https://doi.org/10.1093/bioinformatics/btv359

* **Immcantation / Scoper (Hierarchical Clonotyping):**
  > Gupta, N., Adams, K., Briggs, A., Timberlake, S., Vigneault, F., & Kleinstein, S. (2017). Hierarchical clustering can identify B cell clones with high confidence in Ig repertoire sequencing data. *The Journal of Immunology*, 2489-2499. https://doi.org/10.4049/jimmunol.1601850

---

## 📜 License
This code is released under the **BSD 3-Clause License**.
