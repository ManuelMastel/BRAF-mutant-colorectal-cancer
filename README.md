## Example analysis files

In addition to the figure-specific analysis scripts, this repository will include three exemplary analysis files that demonstrate how the main sequencing datasets from the manuscript can be accessed and processed.

These example files correspond to the major sequencing-based analyses presented in the paper:

- `example_bulk_RNAseq_analysis.R`  
  Exemplary bulk RNA-seq analysis workflow for the tumor and organoid RNA-seq datasets used mainly in **Figure 2** and **Figure 3**.

- `example_CITEseq_analysis.R`  
  Exemplary CITE-seq analysis workflow for the immune microenvironment analyses shown mainly in **Figure 4** and **Figure S4**.

- `example_scRNAseq_analysis.R`  
  Exemplary single-cell RNA-seq analysis workflow for the single-cell validation analyses shown mainly in **Figure 6**.

The examples are intended to provide practical starting points for reusing the sequencing data from:

**Mastel et al. _WNT-driven immune evasion promotes malignant transformation of BRAF-mutant colorectal cancer_.**

All raw data and processed data objects are available through **ArrayExpress** and the **European Nucleotide Archive (ENA)** under the following accession numbers:

- `E-MTAB-17076` — RNA-seq of tumors
- `E-MTAB-17099` — RNA-seq of organoids
- `E-MTAB-17099` — CITE-seq
- `E-MTAB-17103` — exon sequencing
- `E-MTAB-17106` — scRNA-seq

For questions about the code or data, please contact:

**Manuel Mastel**  
German Cancer Research Center / HI-STEM  
manuel.mastel@dkfz-heidelberg.de
