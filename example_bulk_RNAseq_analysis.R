#!/usr/bin/env Rscript

################################################################################
# Example bulk RNA-seq analysis
#
# Manuscript:
# Mastel et al.
# WNT-driven immune evasion promotes malignant transformation of
# BRAF-mutant colorectal cancer
#
# This script provides an exemplary workflow for the analysis of bulk RNA-seq
# data from organoid models used in the manuscript. It is intended as a reusable
# example rather than a full figure-by-figure reproduction script.
#
# The workflow includes:
#   1. Loading raw count data and sample annotation
#   2. Creating a DESeq2 object
#   3. Filtering lowly expressed genes
#   4. Running DESeq2 normalization and differential expression analysis
#   5. Exporting normalized counts
#   6. Performing variance-stabilizing transformation
#   7. Creating PCA plots
#   8. Creating exemplary differential expression result tables and volcano plots
#
# Input files:
#   - VBPNX_organoid_raw_counts_matrix.txt
#   - VBPNX_organoid_DESeq2_normalized_counts_matrix.txt
#   - VBPNX_organoid_sample_annotation.txt
#
# Data availability:
# Raw data and processed data objects are available through ArrayExpress and
# the European Nucleotide Archive (ENA).
#
# Relevant accession for organoid RNA-seq:
#   - E-MTAB-17099
#
# Contact:
# Manuel Mastel
# German Cancer Research Center / HI-STEM
# manuel.mastel@dkfz-heidelberg.de
################################################################################


################################################################################
# 1. Install packages if required
################################################################################

# This section is intentionally commented out.
# Run these commands once if the required packages are not installed.

# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("vsn")
# BiocManager::install("pheatmap")
# BiocManager::install("RColorBrewer")
# install.packages("tidyverse")
# install.packages("ggpubr")


################################################################################
# 2. Load packages
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(vsn)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(ggpubr)
})


################################################################################
# 3. Define input and output paths
################################################################################

# Set this to the folder containing the input files.
# If the files are in the same folder as this script, "." is sufficient.
data_dir <- "."

raw_counts_file <- file.path(data_dir, "VBPNX_organoid_raw_counts_matrix.txt")
normalized_counts_file <- file.path(data_dir, "VBPNX_organoid_DESeq2_normalized_counts_matrix.txt")
sample_annotation_file <- file.path(data_dir, "VBPNX_organoid_sample_annotation.txt")

output_dir <- file.path(data_dir, "bulk_RNAseq_example_output")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


################################################################################
# 4. Load input data
################################################################################

# Raw count matrix:
# Rows should correspond to genes.
# Columns should correspond to samples.
#
# The first column is expected to contain gene identifiers or gene names.
# If your file has a different structure, adjust the column handling below.

raw_counts <- read.delim(
  raw_counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_annotation <- read.delim(
  sample_annotation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Optional: load precomputed normalized counts for comparison or downstream use.
# The DESeq2-normalized counts are also generated again below from the raw counts.
precomputed_normalized_counts <- read.delim(
  normalized_counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


################################################################################
# 5. Prepare count matrix and sample annotation
################################################################################

# Use the first column as gene identifier.
# Common column names include "gene", "Gene", "gene_id", "Geneid", or "name".
# This generic approach assumes that the first column contains gene identifiers.

gene_id_column <- colnames(raw_counts)[1]

count_matrix <- raw_counts %>%
  distinct(.data[[gene_id_column]], .keep_all = TRUE) %>%
  column_to_rownames(gene_id_column)

# Convert to matrix and ensure numeric storage.
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "numeric"

# DESeq2 requires integer raw counts.
count_matrix <- round(count_matrix)

# Sample annotation must contain one row per sample.
# The sample ID column should match the column names of the count matrix.
#
# This script tries to detect a sample ID column automatically.
# If this fails, manually set sample_id_column below.

possible_sample_columns <- c("sample", "Sample", "sample_id", "SampleID", "Cells", "cell", "Cell")

sample_id_column <- possible_sample_columns[
  possible_sample_columns %in% colnames(sample_annotation)
][1]

if (is.na(sample_id_column)) {
  stop(
    "Could not detect a sample ID column in the sample annotation file. ",
    "Please rename the sample ID column to one of: ",
    paste(possible_sample_columns, collapse = ", ")
  )
}

sample_annotation <- sample_annotation %>%
  column_to_rownames(sample_id_column)

# Keep only samples present in both files and order annotation according to count matrix.
common_samples <- intersect(colnames(count_matrix), rownames(sample_annotation))

if (length(common_samples) == 0) {
  stop("No matching sample names found between count matrix and sample annotation.")
}

count_matrix <- count_matrix[, common_samples, drop = FALSE]
sample_annotation <- sample_annotation[common_samples, , drop = FALSE]

# Check that sample order matches.
stopifnot(all(colnames(count_matrix) == rownames(sample_annotation)))


################################################################################
# 6. Define experimental variables
################################################################################

# The sample annotation file should contain a column describing the biological
# condition or genotype to compare.
#
# This script tries to detect the condition column automatically.
# If needed, manually set condition_column below.

possible_condition_columns <- c("Condition", "condition", "Group", "group", "Genotype", "genotype")

condition_column <- possible_condition_columns[
  possible_condition_columns %in% colnames(sample_annotation)
][1]

if (is.na(condition_column)) {
  stop(
    "Could not detect a condition column in the sample annotation file. ",
    "Please rename the condition column to one of: ",
    paste(possible_condition_columns, collapse = ", ")
  )
}

sample_annotation[[condition_column]] <- factor(sample_annotation[[condition_column]])

# Optional: set reference condition manually.
# Adjust this if your control group has a different name.
#
# Example:
# sample_annotation[[condition_column]] <- relevel(
#   sample_annotation[[condition_column]],
#   ref = "CTRL"
# )


################################################################################
# 7. Create DESeq2 object
################################################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_annotation,
  design = as.formula(paste("~", condition_column))
)

# Remove genes with very low total counts.
dds <- dds[rowSums(counts(dds)) > 4, ]

# Run DESeq2.
dds <- DESeq(dds)


################################################################################
# 8. Export DESeq2-normalized counts
################################################################################

normalized_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

write.table(
  normalized_counts,
  file = file.path(output_dir, "DESeq2_normalized_counts_from_raw_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


################################################################################
# 9. Variance-stabilizing transformation
################################################################################

vsd <- vst(dds, blind = FALSE)

vsd_matrix <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

write.table(
  vsd_matrix,
  file = file.path(output_dir, "DESeq2_vst_transformed_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Mean-SD plot to assess variance stabilization.
pdf(file.path(output_dir, "mean_sd_plot_vst.pdf"), width = 6, height = 5)
meanSdPlot(assay(vsd), ranks = FALSE)
dev.off()


################################################################################
# 10. PCA analysis
################################################################################

pca_data <- plotPCA(
  vsd,
  intgroup = condition_column,
  ntop = 500,
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    fill = .data[[condition_column]]
  )
) +
  stat_ellipse(
    aes(color = .data[[condition_column]]),
    geom = "polygon",
    alpha = 0.18,
    level = 0.6,
    show.legend = FALSE
  ) +
  geom_point(
    shape = 21,
    size = 4,
    color = "black",
    stroke = 0.5
  ) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
  )

ggsave(
  filename = file.path(output_dir, "PCA_condition.pdf"),
  plot = pca_plot,
  width = 6,
  height = 5
)

ggsave(
  filename = file.path(output_dir, "PCA_condition.png"),
  plot = pca_plot,
  width = 6,
  height = 5,
  dpi = 300
)


################################################################################
# 11. Optional PCA by clone or additional annotation variable
################################################################################

# If the annotation file contains a clone column, create an additional PCA plot.

possible_clone_columns <- c("Clone", "clone", "Organoid", "organoid")

clone_column <- possible_clone_columns[
  possible_clone_columns %in% colnames(sample_annotation)
][1]

if (!is.na(clone_column)) {
  pca_data_clone <- plotPCA(
    vsd,
    intgroup = c(condition_column, clone_column),
    ntop = 500,
    returnData = TRUE
  )
  
  percent_var_clone <- round(100 * attr(pca_data_clone, "percentVar"))
  
  pca_clone_plot <- ggplot(
    pca_data_clone,
    aes(
      x = PC1,
      y = PC2,
      fill = .data[[condition_column]]
    )
  ) +
    geom_point(
      aes(shape = .data[[clone_column]]),
      size = 4,
      color = "black",
      stroke = 0.5
    ) +
    xlab(paste0("PC1: ", percent_var_clone[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var_clone[2], "% variance")) +
    theme_minimal(base_size = 13) +
    theme(
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
  
  ggsave(
    filename = file.path(output_dir, "PCA_condition_clone.pdf"),
    plot = pca_clone_plot,
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path(output_dir, "PCA_condition_clone.png"),
    plot = pca_clone_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
}


################################################################################
# 12. Differential expression analysis
################################################################################

# This section creates exemplary pairwise comparisons between all conditions.
# For manuscript-specific comparisons, select the biologically relevant contrasts.

condition_levels <- levels(sample_annotation[[condition_column]])

if (length(condition_levels) >= 2) {
  contrast_table <- combn(condition_levels, 2, simplify = FALSE)
  
  for (contrast_pair in contrast_table) {
    reference_condition <- contrast_pair[1]
    test_condition <- contrast_pair[2]
    
    message("Running contrast: ", test_condition, " vs ", reference_condition)
    
    res <- results(
      dds,
      contrast = c(condition_column, test_condition, reference_condition)
    )
    
    res <- lfcShrink(
      dds,
      contrast = c(condition_column, test_condition, reference_condition),
      res = res,
      type = "ashr"
    )
    
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      arrange(padj)
    
    contrast_name <- paste0(test_condition, "_vs_", reference_condition)
    
    write.table(
      res_df,
      file = file.path(output_dir, paste0("DESeq2_results_", contrast_name, ".txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    volcano_plot <- EnhancedVolcano(
      res_df,
      lab = res_df$gene_id,
      x = "log2FoldChange",
      y = "padj",
      title = contrast_name,
      subtitle = "DESeq2 differential expression",
      xlab = bquote(~Log[2] ~ "fold change"),
      ylab = bquote(~-Log[10] ~ "adjusted P value"),
      pCutoff = 0.05,
      FCcutoff = 1,
      pointSize = 2.0,
      labSize = 3.0,
      legendPosition = "right",
      legendLabSize = 10,
      legendIconSize = 3.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5
    )
    
    ggsave(
      filename = file.path(output_dir, paste0("Volcano_", contrast_name, ".pdf")),
      plot = volcano_plot,
      width = 8,
      height = 7
    )
    
    ggsave(
      filename = file.path(output_dir, paste0("Volcano_", contrast_name, ".png")),
      plot = volcano_plot,
      width = 8,
      height = 7,
      dpi = 300
    )
  }
}


################################################################################
# 13. Heatmap of most variable genes
################################################################################

# Select the 500 most variable genes based on the VST-transformed matrix.

vst_assay <- assay(vsd)

top_variable_genes <- vst_assay %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  mutate(variance = apply(vst_assay, 1, var)) %>%
  arrange(desc(variance)) %>%
  slice_head(n = min(500, nrow(.))) %>%
  pull(gene_id)

heatmap_matrix <- vst_assay[top_variable_genes, , drop = FALSE]

annotation_col <- sample_annotation %>%
  select(all_of(condition_column))

if (!is.na(clone_column)) {
  annotation_col <- sample_annotation %>%
    select(all_of(c(condition_column, clone_column)))
}

pdf(file.path(output_dir, "heatmap_top_variable_genes.pdf"), width = 8, height = 10)
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top variable genes"
)
dev.off()


################################################################################
# 14. Save R session information
################################################################################

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("Bulk RNA-seq example analysis completed successfully.")
message("Results written to: ", output_dir)