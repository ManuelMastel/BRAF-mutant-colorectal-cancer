# Code accompanying: WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer

This repository contains the analysis code accompanying the manuscript:

**Mastel et al. _WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer_**

The study investigates how WNT-pathway activation promotes malignant transformation and immune evasion in BRAF-mutant microsatellite-stable colorectal cancer. The analyses include genetically engineered mouse models, organoid models, bulk RNA sequencing, whole-exome sequencing, single-cell RNA sequencing, CITE-seq, drug screening, and public human colorectal cancer datasets.

## Repository rationale

Rather than providing a detailed standalone analysis script for every individual figure panel, this repository provides exemplary analysis workflows for the major sequencing-based methods used in the study.

The aim is to make the code reusable and easier to follow by providing one representative analysis script per experimental method. These scripts demonstrate how the respective data types can be read in, processed, analyzed, and used for downstream visualization.

## Example analysis files

This repository includes three exemplary analysis files:

- `example_bulk_RNAseq_analysis.R`  
  Exemplary bulk RNA-seq analysis workflow for the tumor and organoid RNA-seq datasets used mainly in **Figure 2** and **Figure 3**.

- `example_CITEseq_analysis.R`  
  Exemplary CITE-seq analysis workflow for the immune microenvironment analyses shown mainly in **Figure 4** and **Figure S4**.

- `example_scRNAseq_analysis.R`  
  Exemplary single-cell RNA-seq analysis workflow for the single-cell validation analyses shown mainly in **Figure 6**.

These examples are intended to provide practical starting points for reusing the sequencing data from:

**Mastel et al. _WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer_.**

## Repository structure

```text
.
├── README.md
├── example_bulk_RNAseq_analysis.R
├── example_CITEseq_analysis.R
└── example_scRNAseq_analysis.R
```

The example scripts are not intended to reproduce every figure panel exactly. Instead, they provide method-specific workflows that can be adapted to reproduce or extend the corresponding analyses from the manuscript.

## Requirements

The analyses were performed in R. Required packages may differ between the individual example scripts.

Commonly used R packages include:

```r
tidyverse
ggplot2
dplyr
readr
DESeq2
fgsea
pheatmap
ComplexHeatmap
survival
survminer
Seurat
patchwork
maftools
pracma
```

Additional packages may be loaded within the individual scripts depending on the specific analysis workflow.

## Input data

This repository contains code for exemplary analysis workflows. Large primary datasets and processed sequencing objects are not included directly in the repository.

The scripts may require one or more of the following input files:

- raw or processed bulk RNA-seq count matrices
- sample metadata tables
- processed single-cell RNA-seq or CITE-seq count matrices
- processed Seurat objects
- antibody-derived tag count matrices for CITE-seq analyses
- whole-exome sequencing summary tables
- mutation annotation files
- public human colorectal cancer single-cell atlas data

All raw data and processed data objects are available through **ArrayExpress** and the **European Nucleotide Archive (ENA)** under the following accession numbers:

- `E-MTAB-17076` — RNA-seq of tumors
- `E-MTAB-17099` — RNA-seq of organoids
- `E-MTAB-17099` — CITE-seq
- `E-MTAB-17103` — exon sequencing
- `E-MTAB-17106` — scRNA-seq

The exemplary scripts demonstrate how data can be read in from these repositories and processed for downstream analysis.

## How to run the code

After downloading the required data from ArrayExpress or ENA, adjust the file paths in the corresponding example script.

Example:

```bash
Rscript example_bulk_RNAseq_analysis.R
```

Alternatively, open the `.R` file in RStudio and run the code interactively.

## Reproducibility notes

- The repository provides exemplary workflows for major sequencing-based methods rather than figure-by-figure reproduction scripts.
- The scripts are intended as practical starting points for reproducing and extending the analyses from the manuscript.
- File paths may need to be adjusted depending on the local location of the downloaded input data.
- Large sequencing-derived objects are not stored directly in this repository.
- Random seeds are set in scripts where stochastic analyses are performed.
- Final publication panels may differ from direct script output if figures were assembled, resized, or annotated using external layout software.

## Citation

If you use this code, please cite:

**Mastel et al. _WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer_.**

## Contact

For questions about the code or data, please contact:

**Manuel Mastel**  
German Cancer Research Center / HI-STEM  
manuel.mastel@dkfz-heidelberg.de
