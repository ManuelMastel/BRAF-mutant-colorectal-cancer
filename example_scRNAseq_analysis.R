#!/usr/bin/env Rscript

################################################################################
# Example single-cell RNA-seq analysis
#
# Manuscript:
# Mastel et al.
# WNT-driven immune evasion promotes malignant transformation of
# BRAF-mutant colorectal cancer
#
# This script provides an exemplary workflow for the analysis of single-cell
# RNA-seq data from the T cell depletion experiment used in the manuscript.
#
# The workflow includes:
#   1. Loading gene expression counts, cell barcodes, gene features, hash counts,
#      hash features, and cell metadata
#   2. Creating a Seurat object
#   3. Adding hashtag oligo / antibody capture data
#   4. Performing RNA quality control
#   5. Using available hash classifications to retain singlets
#   6. Normalization, variable feature detection, scaling, PCA, UMAP, and clustering
#   7. Visualization of QC metrics, clusters, treatment/sample groups, and marker genes
#   8. Exemplary cell type annotation based on canonical marker genes
#   9. Saving the processed Seurat object and output plots
#
# Input files:
#   - VBPNA_TCD_gene_count_matrix_dense.txt
#   - VBPNA_TCD_barcodes.txt
#   - VBPNA_TCD_gene_features.txt
#   - VBPNA_TCD_cell_metadata_with_hash_counts.txt
#   - VBPNA_TCD_cell_to_sample_metadata.txt
#   - VBPNA_TCD_hash_count_matrix.txt
#   - VBPNA_TCD_hash_features.txt
#
# Data availability:
# Raw data and processed data objects are available through ArrayExpress and
# the European Nucleotide Archive (ENA).
#
# Relevant accession for scRNA-seq:
#   - E-MTAB-17106
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

# install.packages("Seurat")
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("ggpubr")
# install.packages("Matrix")
# install.packages("data.table")


################################################################################
# 2. Load packages
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggpubr)
  library(Matrix)
  library(data.table)
})


################################################################################
# 3. Define input and output paths
################################################################################

# Set this to the folder containing the input files.
# If the files are in the same folder as this script, "." is sufficient.
data_dir <- "."

gene_count_file <- file.path(data_dir, "VBPNA_TCD_gene_count_matrix_dense.txt")
barcode_file <- file.path(data_dir, "VBPNA_TCD_barcodes.txt")
gene_feature_file <- file.path(data_dir, "VBPNA_TCD_gene_features.txt")
cell_metadata_file <- file.path(data_dir, "VBPNA_TCD_cell_metadata_with_hash_counts.txt")
cell_to_sample_file <- file.path(data_dir, "VBPNA_TCD_cell_to_sample_metadata.txt")
hash_count_file <- file.path(data_dir, "VBPNA_TCD_hash_count_matrix.txt")
hash_feature_file <- file.path(data_dir, "VBPNA_TCD_hash_features.txt")

output_dir <- file.path(data_dir, "scRNAseq_example_output")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


################################################################################
# 4. Helper functions
################################################################################

read_count_matrix <- function(file, feature_names = NULL, barcode_names = NULL) {
  counts <- fread(file, data.table = FALSE, check.names = FALSE)
  
  # Case 1:
  # First column contains feature names and remaining columns are cell barcodes.
  if (!all(colnames(counts)[1] %in% barcode_names) && !is.numeric(counts[[1]])) {
    rownames(counts) <- counts[[1]]
    counts <- counts[, -1, drop = FALSE]
  }
  
  # Case 2:
  # Matrix has no feature-name column, so use supplied feature names.
  if (!is.null(feature_names) && nrow(counts) == length(feature_names)) {
    rownames(counts) <- feature_names
  }
  
  # If barcode names are supplied and dimensions match, use them as column names.
  if (!is.null(barcode_names) && ncol(counts) == length(barcode_names)) {
    colnames(counts) <- barcode_names
  }
  
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  counts <- Matrix(counts, sparse = TRUE)
  
  return(counts)
}


save_plot <- function(plot, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(output_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}


################################################################################
# 5. Load input tables
################################################################################

barcodes <- read.delim(
  barcode_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_features <- read.delim(
  gene_feature_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_metadata <- read.delim(
  cell_metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cell_to_sample <- read.delim(
  cell_to_sample_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

hash_features <- read.delim(
  hash_feature_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


################################################################################
# 6. Prepare barcode, gene, and hash identifiers
################################################################################

# Barcode file is expected to contain a column named "barcode".
barcode_column <- "barcode"

if (!barcode_column %in% colnames(barcodes)) {
  stop("The barcode file must contain a column named 'barcode'.")
}

cell_barcodes <- barcodes[[barcode_column]]

# Gene feature file is expected to contain feature_id and feature_name.
# Use feature_name for readability, but make names unique in case of duplicates.
if ("feature_name" %in% colnames(gene_features)) {
  gene_names <- make.unique(gene_features$feature_name)
} else if ("feature_id" %in% colnames(gene_features)) {
  gene_names <- make.unique(gene_features$feature_id)
} else {
  stop("The gene feature file must contain 'feature_name' or 'feature_id'.")
}

# Hash feature file is expected to contain hash_id and hash_name.
if ("hash_name" %in% colnames(hash_features)) {
  hash_names <- make.unique(hash_features$hash_name)
} else if ("hash_id" %in% colnames(hash_features)) {
  hash_names <- make.unique(hash_features$hash_id)
} else {
  stop("The hash feature file must contain 'hash_name' or 'hash_id'.")
}


################################################################################
# 7. Load gene expression and hash count matrices
################################################################################

gene_counts <- read_count_matrix(
  file = gene_count_file,
  feature_names = gene_names,
  barcode_names = cell_barcodes
)

hash_counts <- read_count_matrix(
  file = hash_count_file,
  feature_names = hash_names,
  barcode_names = cell_barcodes
)

# Make sure RNA and HTO matrices contain the same cells.
common_cells <- intersect(colnames(gene_counts), colnames(hash_counts))

if (length(common_cells) == 0) {
  stop("No overlapping barcodes found between RNA and hash count matrices.")
}

gene_counts <- gene_counts[, common_cells, drop = FALSE]
hash_counts <- hash_counts[, common_cells, drop = FALSE]


################################################################################
# 8. Create Seurat object
################################################################################

sc <- CreateSeuratObject(
  counts = gene_counts,
  project = "VBPNA_TCD",
  min.cells = 3,
  min.features = 200
)

# Add HTO / hashtag count data as a separate assay.
hash_counts <- hash_counts[, colnames(sc), drop = FALSE]
sc[["HTO"]] <- CreateAssayObject(counts = hash_counts)


################################################################################
# 9. Add cell metadata
################################################################################

# The metadata files are expected to contain a barcode column.
if (!"barcode" %in% colnames(cell_metadata)) {
  stop("The cell metadata file must contain a column named 'barcode'.")
}

if (!"barcode" %in% colnames(cell_to_sample)) {
  stop("The cell-to-sample metadata file must contain a column named 'barcode'.")
}

cell_metadata <- cell_metadata %>%
  distinct(barcode, .keep_all = TRUE) %>%
  column_to_rownames("barcode")

cell_to_sample <- cell_to_sample %>%
  distinct(barcode, .keep_all = TRUE) %>%
  column_to_rownames("barcode")

# Add metadata only for cells present in the Seurat object.
metadata_cells <- intersect(colnames(sc), rownames(cell_metadata))
sc <- AddMetaData(sc, metadata = cell_metadata[metadata_cells, , drop = FALSE])

sample_metadata_cells <- intersect(colnames(sc), rownames(cell_to_sample))
sc <- AddMetaData(sc, metadata = cell_to_sample[sample_metadata_cells, , drop = FALSE])


################################################################################
# 10. RNA quality control
################################################################################

DefaultAssay(sc) <- "RNA"

# Mouse mitochondrial genes usually start with "mt-".
# For human data, use "^MT-".
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")

qc_plot <- VlnPlot(
  sc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

save_plot(qc_plot, "QC_violin_before_filtering", width = 10, height = 4)

# Basic QC thresholds.
# These may need to be adjusted depending on the dataset.
min_features <- 500
max_features <- 10000
max_percent_mt <- 25

sc_filtered <- subset(
  sc,
  subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_percent_mt
)

qc_plot_filtered <- VlnPlot(
  sc_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

save_plot(qc_plot_filtered, "QC_violin_after_filtering", width = 10, height = 4)


################################################################################
# 11. Normalize HTO data and inspect hash classifications
################################################################################

DefaultAssay(sc_filtered) <- "HTO"

sc_filtered <- NormalizeData(
  sc_filtered,
  assay = "HTO",
  normalization.method = "CLR"
)

# The provided metadata already contains hash classifications.
# If hash_classification is present, use it to identify singlets.
# Otherwise, perform HTODemux.

if ("hash_classification" %in% colnames(sc_filtered@meta.data)) {
  message("Using hash classifications from metadata.")
  
  table_hash <- table(sc_filtered$hash_classification)
  write.table(
    as.data.frame(table_hash),
    file = file.path(output_dir, "hash_classification_table.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  message("No hash_classification column found. Running HTODemux.")
  
  sc_filtered <- HTODemux(
    sc_filtered,
    assay = "HTO",
    positive.quantile = 0.99
  )
  
  table_hash <- table(sc_filtered$HTO_classification.global)
  write.table(
    as.data.frame(table_hash),
    file = file.path(output_dir, "HTODemux_classification_table.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# HTO heatmap for visual inspection.
DefaultAssay(sc_filtered) <- "HTO"

hto_heatmap_file <- file.path(output_dir, "HTO_heatmap.pdf")
pdf(hto_heatmap_file, width = 8, height = 6)
print(HTOHeatmap(sc_filtered, assay = "HTO", ncells = min(5000, ncol(sc_filtered))))
dev.off()


################################################################################
# 12. Retain singlets
################################################################################

# Prefer the provided hash classification if available.
if ("hash_classification" %in% colnames(sc_filtered@meta.data)) {
  sc_singlet <- subset(sc_filtered, subset = hash_classification == "Singlet")
} else {
  Idents(sc_filtered) <- "HTO_classification.global"
  sc_singlet <- subset(sc_filtered, idents = "Singlet")
}

DefaultAssay(sc_singlet) <- "RNA"

message("Number of cells before singlet filtering: ", ncol(sc_filtered))
message("Number of singlet cells retained: ", ncol(sc_singlet))


################################################################################
# 13. Standard scRNA-seq preprocessing
################################################################################

sc_singlet <- NormalizeData(
  sc_singlet,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

sc_singlet <- FindVariableFeatures(
  sc_singlet,
  selection.method = "vst",
  nfeatures = 2000
)

top_variable_features <- head(VariableFeatures(sc_singlet), 20)

variable_feature_plot <- VariableFeaturePlot(sc_singlet)
variable_feature_plot <- LabelPoints(
  plot = variable_feature_plot,
  points = top_variable_features,
  repel = TRUE
)

save_plot(variable_feature_plot, "variable_features", width = 8, height = 5)

sc_singlet <- ScaleData(
  sc_singlet,
  features = VariableFeatures(sc_singlet)
)

sc_singlet <- RunPCA(
  sc_singlet,
  features = VariableFeatures(sc_singlet)
)

elbow_plot <- ElbowPlot(sc_singlet, ndims = 30)
save_plot(elbow_plot, "PCA_elbow_plot", width = 6, height = 5)


################################################################################
# 14. Clustering and dimensional reduction
################################################################################

n_pcs <- 10

sc_singlet <- FindNeighbors(
  sc_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)

sc_singlet <- FindClusters(
  sc_singlet,
  resolution = 0.4
)

sc_singlet <- RunUMAP(
  sc_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)

sc_singlet <- RunTSNE(
  sc_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)


################################################################################
# 15. Basic visualization
################################################################################

umap_clusters <- DimPlot(
  sc_singlet,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("Seurat clusters")

save_plot(umap_clusters, "UMAP_seurat_clusters", width = 7, height = 5)

tsne_clusters <- DimPlot(
  sc_singlet,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("Seurat clusters")

save_plot(tsne_clusters, "TSNE_seurat_clusters", width = 7, height = 5)

# Visualize sample or hash assignments if available.
metadata_columns_to_plot <- c(
  "sample",
  "hash_classification",
  "detailed_hash_classification",
  "HTO_classification",
  "HTO_classification.global"
)

metadata_columns_to_plot <- metadata_columns_to_plot[
  metadata_columns_to_plot %in% colnames(sc_singlet@meta.data)
]

for (metadata_column in metadata_columns_to_plot) {
  plot <- DimPlot(
    sc_singlet,
    reduction = "umap",
    group.by = metadata_column
  ) +
    ggtitle(metadata_column)
  
  save_plot(
    plot,
    paste0("UMAP_", metadata_column),
    width = 8,
    height = 5
  )
}


################################################################################
# 16. Marker gene visualization
################################################################################

# Marker genes shown here are examples and should be adapted to the question.
# Common markers:
#   Epcam   - epithelial / cancer cells
#   Krt8    - epithelial cells
#   Krt18   - epithelial cells
#   Ptprc   - immune cells
#   Cd3d    - T cells
#   Cd3e    - T cells
#   Cd4     - CD4 T cells
#   Cd8a    - CD8 T cells
#   Lyz2    - myeloid cells
#   Chil3   - myeloid / inflammatory myeloid cells
#   Col1a1  - fibroblasts
#   Pecam1  - endothelial cells

marker_genes <- c(
  "Epcam", "Krt8", "Krt18",
  "Ptprc", "Cd3d", "Cd3e", "Cd4", "Cd8a",
  "Lyz2", "Chil3",
  "Col1a1", "Pecam1"
)

marker_genes <- marker_genes[marker_genes %in% rownames(sc_singlet)]

if (length(marker_genes) > 0) {
  feature_plot <- FeaturePlot(
    sc_singlet,
    features = marker_genes,
    reduction = "umap",
    ncol = 3
  )
  
  save_plot(
    feature_plot,
    "UMAP_marker_genes",
    width = 12,
    height = ceiling(length(marker_genes) / 3) * 3
  )
  
  dot_plot <- DotPlot(
    sc_singlet,
    features = marker_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis()
  
  save_plot(
    dot_plot,
    "DotPlot_marker_genes",
    width = 10,
    height = 5
  )
}


################################################################################
# 17. Find cluster marker genes
################################################################################

cluster_markers <- FindAllMarkers(
  sc_singlet,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(
  cluster_markers,
  file = file.path(output_dir, "cluster_markers_FindAllMarkers.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.table(
  top10_markers,
  file = file.path(output_dir, "cluster_markers_top10_per_cluster.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


################################################################################
# 18. Exemplary cell type annotation
################################################################################

# This annotation is deliberately simple and should be reviewed manually.
# Adjust cluster-to-cell-type assignments based on marker gene expression.

sc_singlet$seurat_clusters <- as.character(sc_singlet$seurat_clusters)

# Example cluster mapping.
# Replace these assignments after inspecting marker genes.
cluster_to_celltype <- c(
  "0" = "Cancer cells",
  "1" = "Cancer cells",
  "2" = "Cancer cells",
  "3" = "Cancer cells",
  "4" = "Cancer cells",
  "5" = "Cancer cells",
  "6" = "Myeloid cells",
  "7" = "Fibroblasts",
  "8" = "T cells",
  "9" = "Endothelial cells"
)

sc_singlet$celltype <- cluster_to_celltype[sc_singlet$seurat_clusters]
sc_singlet$celltype[is.na(sc_singlet$celltype)] <- "Unassigned"

umap_celltype <- DimPlot(
  sc_singlet,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Exemplary cell type annotation")

save_plot(umap_celltype, "UMAP_celltype_annotation", width = 8, height = 5)


################################################################################
# 19. Summarize cell numbers
################################################################################

cell_number_cluster <- as.data.frame(table(sc_singlet$seurat_clusters))
colnames(cell_number_cluster) <- c("cluster", "n_cells")

write.table(
  cell_number_cluster,
  file = file.path(output_dir, "cell_numbers_by_cluster.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cell_number_celltype <- as.data.frame(table(sc_singlet$celltype))
colnames(cell_number_celltype) <- c("celltype", "n_cells")

write.table(
  cell_number_celltype,
  file = file.path(output_dir, "cell_numbers_by_celltype.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if ("sample" %in% colnames(sc_singlet@meta.data)) {
  cell_number_sample <- as.data.frame(table(sc_singlet$sample))
  colnames(cell_number_sample) <- c("sample", "n_cells")
  
  write.table(
    cell_number_sample,
    file = file.path(output_dir, "cell_numbers_by_sample.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


################################################################################
# 20. Save processed object and session information
################################################################################

saveRDS(
  sc_singlet,
  file = file.path(output_dir, "VBPNA_TCD_scRNAseq_singlet_processed_seurat_object.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("scRNA-seq example analysis completed successfully.")
message("Results written to: ", output_dir)