# BRAF-mutant-colorectal-cancer

# Code accompanying: WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer

This repository contains the code used to generate the figures for the manuscript:

**Mastel et al. _WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer_**

The study investigates how WNT-pathway activation promotes malignant transformation and immune evasion in BRAF-mutant microsatellite-stable colorectal cancer using genetically engineered mouse models, organoid models, bulk RNA sequencing, whole-exome sequencing, single-cell RNA sequencing, CITE-seq, drug screening, and public human colorectal cancer datasets.

## Repository structure

Each figure is organized in a separate folder. Every folder contains the corresponding R analysis script and a rendered HTML file.

```text
.
├── README.md
├── Figure_1/
│   ├── Figure_1_analysis.R
│   └── Figure_1_analysis.html
├── Figure_2/
│   ├── Figure_2_analysis.R
│   └── Figure_2_analysis.html
├── Figure_3/
│   ├── Figure_3_analysis.R
│   └── Figure_3_analysis.html
├── Figure_4/
│   ├── Figure_4_analysis.R
│   └── Figure_4_analysis.html
├── Figure_5/
│   ├── Figure_5_analysis.R
│   └── Figure_5_analysis.html
├── Figure_6/
│   ├── Figure_6_analysis.R
│   └── Figure_6_analysis.html
├── Figure_S1/
│   ├── Figure_S1_analysis.R
│   └── Figure_S1_analysis.html
├── Figure_S2/
│   ├── Figure_S2_analysis.R
│   └── Figure_S2_analysis.html
├── Figure_S3/
│   ├── Figure_S3_analysis.R
│   └── Figure_S3_analysis.html
├── Figure_S4/
│   ├── Figure_S4_analysis.R
│   └── Figure_S4_analysis.html
└── Figure_S5/
    ├── Figure_S5_analysis.R
    └── Figure_S5_analysis.html


The .R files contain the executable analysis and plotting code.
The .html files are rendered reports that allow inspection of the analysis output without rerunning the code.

Figure overview
Folder	Description
Figure_1/	Genetic modeling of BRAF-mutant colorectal cancer, survival analyses, WNT activation, and mutational profiling
Figure_2/	WNT-pathway alterations, tumor initiation analyses, bulk RNA-seq, CMS classification, and pathway enrichment
Figure_3/	Isogenic organoid models, WNT activation, immune-related gene expression, transplantation assays, and drug response analyses
Figure_4/	CITE-seq and single-cell RNA-seq analysis of tumor immune microenvironment remodeling and T cell depletion experiments
Figure_5/	WNT-mediated regulation of CCL20 and functional Ccl20 overexpression experiments
Figure_6/	Human CRC survival analyses and single-cell validation of the CCL20–CCR6 axis
Figure_S1/	Additional phenotypic and molecular characterization of BRAF-mutant tumors
Figure_S2/	Genome editing validation and characterization of BRAF-mutant organoid models
Figure_S3/	Drug-screening analysis of BRAF-mutant organoids
Figure_S4/	Additional CITE-seq characterization of myeloid and T cell populations
Figure_S5/	Additional analyses of T cell depletion experiments
Requirements

The analyses were performed in R. Required packages may differ between individual figure scripts.

Commonly used R packages include:

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

Additional packages may be loaded within individual figure scripts.

Input data

This repository contains the code used for analysis and figure generation. Large primary datasets are not included directly in this repository.

Depending on the figure, the scripts may require one or more of the following input files:

processed bulk RNA-seq count matrices
processed single-cell RNA-seq or CITE-seq Seurat objects
whole-exome sequencing summary tables
mutation annotation files
drug-screening response matrices
qPCR summary tables
tumor-growth measurements
survival-analysis input tables
public human colorectal cancer mutation datasets
public human colorectal cancer single-cell atlas data

Where possible, the scripts start from processed input tables or objects used for plotting and statistical analysis.

Raw sequencing data and large processed objects should be downloaded from the data repositories listed in the manuscript once available. Data that are not publicly deposited may be requested from the corresponding author.

How to run the code

To reproduce a specific figure, open the corresponding folder and run the R script.

Example:

cd Figure_3
Rscript Figure_3_analysis.R

Alternatively, open the .R file in RStudio and run the code interactively.

The rendered .html files can be opened directly in a web browser to inspect the corresponding analysis output.

Reproducibility notes
Each figure folder contains the code required for the corresponding figure.
File paths may need to be adjusted depending on the local location of the input data.
Large sequencing-derived objects are not stored in this repository.
Random seeds are set in scripts where stochastic analyses are performed.
The final publication panels may differ slightly from the direct script output if figures were assembled, resized, or annotated using external layout software.
Citation

If you use this code, please cite:

Mastel et al. WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer.

Contact

For questions about the code or data, please contact:

Manuel Mastel
German Cancer Research Center / HI-STEM
manuel.mastel@dkfz-heidelberg.de

