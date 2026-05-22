#!/usr/bin/env Rscript

################################################################################
# Example CITE-seq analysis
#
# Manuscript:
# Mastel et al.
# WNT-driven immune evasion promotes malignant transformation of
# BRAF-mutant colorectal cancer
#
# This script provides an exemplary workflow for the processing and analysis of
# CITE-seq data from subcutaneous tumor samples used in the manuscript.
#
# The workflow includes:
#   1. Loading RNA gene expression counts, ADT antibody counts, HTO hash counts,
#      cell barcodes, feature annotations, and cell metadata
#   2. Creating a Seurat object from RNA counts
#   3. Adding ADT and HTO assays
#   4. Performing RNA quality control
#   5. Using hash classifications to retain singlet cells
#   6. RNA-based normalization, variable feature detection, scaling, PCA, UMAP,
#      tSNE, and clustering
#   7. CLR normalization of ADT and HTO assays
#   8. Visualization of RNA clusters, sample identity, hash identity, and ADT
#      marker expression
#   9. Identification of RNA cluster marker genes
#  10. Exemplary cell type annotation using RNA and ADT marker information
#  11. Saving processed output files and the final Seurat object
#
# Input files:
#   - VBPN_C_A_subcut_gene_count_matrix_dense.txt
#   - VBPN_C_A_subcut_ADT_count_matrix.txt
#   - VBPN_C_A_subcut_ADT_features.txt
#   - VBPN_C_A_subcut_barcodes.txt
#   - VBPN_C_A_subcut_cell_metadata_with_hash_counts.txt
#   - VBPN_C_A_subcut_cell_to_sample_metadata.txt
#   - VBPN_C_A_subcut_gene_features.txt
#   - VBPN_C_A_subcut_hash_count_matrix.txt
#   - VBPN_C_A_subcut_hash_features.txt
#
# Data availability:
# Raw data and processed data objects are available through ArrayExpress and
# the European Nucleotide Archive (ENA).
#
# Relevant accession for CITE-seq:
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

gene_count_file <- file.path(data_dir, "VBPN_C_A_subcut_gene_count_matrix_dense.txt")
adt_count_file <- file.path(data_dir, "VBPN_C_A_subcut_ADT_count_matrix.txt")
hash_count_file <- file.path(data_dir, "VBPN_C_A_subcut_hash_count_matrix.txt")

barcode_file <- file.path(data_dir, "VBPN_C_A_subcut_barcodes.txt")
gene_feature_file <- file.path(data_dir, "VBPN_C_A_subcut_gene_features.txt")
adt_feature_file <- file.path(data_dir, "VBPN_C_A_subcut_ADT_features.txt")
hash_feature_file <- file.path(data_dir, "VBPN_C_A_subcut_hash_features.txt")

cell_metadata_file <- file.path(data_dir, "VBPN_C_A_subcut_cell_metadata_with_hash_counts.txt")
cell_to_sample_file <- file.path(data_dir, "VBPN_C_A_subcut_cell_to_sample_metadata.txt")

output_dir <- file.path(data_dir, "CITEseq_example_output")

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
# 5. Load annotation and metadata tables
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

adt_features <- read.delim(
  adt_feature_file,
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


################################################################################
# 6. Prepare barcode and feature identifiers
################################################################################

if (!"barcode" %in% colnames(barcodes)) {
  stop("The barcode file must contain a column named 'barcode'.")
}

cell_barcodes <- barcodes$barcode

if ("feature_name" %in% colnames(gene_features)) {
  gene_names <- make.unique(gene_features$feature_name)
} else if ("feature_id" %in% colnames(gene_features)) {
  gene_names <- make.unique(gene_features$feature_id)
} else {
  stop("The gene feature file must contain 'feature_name' or 'feature_id'.")
}

if ("feature_name" %in% colnames(adt_features)) {
  adt_names <- make.unique(adt_features$feature_name)
} else if ("feature_id" %in% colnames(adt_features)) {
  adt_names <- make.unique(adt_features$feature_id)
} else {
  stop("The ADT feature file must contain 'feature_name' or 'feature_id'.")
}

if ("hash_name" %in% colnames(hash_features)) {
  hash_names <- make.unique(hash_features$hash_name)
} else if ("hash_id" %in% colnames(hash_features)) {
  hash_names <- make.unique(hash_features$hash_id)
} else {
  stop("The hash feature file must contain 'hash_name' or 'hash_id'.")
}


################################################################################
# 7. Load RNA, ADT, and HTO count matrices
################################################################################

gene_counts <- read_count_matrix(
  file = gene_count_file,
  feature_names = gene_names,
  barcode_names = cell_barcodes
)

adt_counts <- read_count_matrix(
  file = adt_count_file,
  feature_names = adt_names,
  barcode_names = cell_barcodes
)

hash_counts <- read_count_matrix(
  file = hash_count_file,
  feature_names = hash_names,
  barcode_names = cell_barcodes
)

# Keep only cells present in all three modalities.
common_cells <- Reduce(
  intersect,
  list(
    colnames(gene_counts),
    colnames(adt_counts),
    colnames(hash_counts)
  )
)

if (length(common_cells) == 0) {
  stop("No overlapping barcodes found between RNA, ADT, and HTO count matrices.")
}

gene_counts <- gene_counts[, common_cells, drop = FALSE]
adt_counts <- adt_counts[, common_cells, drop = FALSE]
hash_counts <- hash_counts[, common_cells, drop = FALSE]


################################################################################
# 8. Create Seurat object and add ADT / HTO assays
################################################################################

cite <- CreateSeuratObject(
  counts = gene_counts,
  project = "VBPN_C_A_subcut_CITEseq",
  min.cells = 3,
  min.features = 200
)

# Align ADT and HTO matrices to the cells retained in the Seurat object.
adt_counts <- adt_counts[, colnames(cite), drop = FALSE]
hash_counts <- hash_counts[, colnames(cite), drop = FALSE]

cite[["ADT"]] <- CreateAssayObject(counts = adt_counts)
cite[["HTO"]] <- CreateAssayObject(counts = hash_counts)


################################################################################
# 9. Add cell metadata
################################################################################

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

metadata_cells <- intersect(colnames(cite), rownames(cell_metadata))
cite <- AddMetaData(cite, metadata = cell_metadata[metadata_cells, , drop = FALSE])

sample_metadata_cells <- intersect(colnames(cite), rownames(cell_to_sample))
cite <- AddMetaData(cite, metadata = cell_to_sample[sample_metadata_cells, , drop = FALSE])


################################################################################
# 10. RNA quality control
################################################################################

DefaultAssay(cite) <- "RNA"

# Mouse mitochondrial genes usually start with "mt-".
# For human data, use "^MT-".
cite[["percent.mt"]] <- PercentageFeatureSet(cite, pattern = "^mt-")

qc_plot_before <- VlnPlot(
  cite,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_ADT", "nCount_HTO"),
  ncol = 5,
  pt.size = 0
)

save_plot(qc_plot_before, "QC_violin_before_filtering", width = 14, height = 4)

# Basic QC thresholds.
# These values are intentionally conservative and should be adjusted after
# inspecting the QC plots for the specific dataset.
min_features <- 200
max_features <- 10000
max_percent_mt <- 25

cite_filtered <- subset(
  cite,
  subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_percent_mt
)

qc_plot_after <- VlnPlot(
  cite_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nCount_ADT", "nCount_HTO"),
  ncol = 5,
  pt.size = 0
)

save_plot(qc_plot_after, "QC_violin_after_filtering", width = 14, height = 4)


################################################################################
# 11. Normalize ADT and HTO assays
################################################################################

DefaultAssay(cite_filtered) <- "ADT"

cite_filtered <- NormalizeData(
  cite_filtered,
  assay = "ADT",
  normalization.method = "CLR",
  margin = 2
)

DefaultAssay(cite_filtered) <- "HTO"

cite_filtered <- NormalizeData(
  cite_filtered,
  assay = "HTO",
  normalization.method = "CLR",
  margin = 2
)


################################################################################
# 12. Inspect or generate hash classifications
################################################################################

# The provided metadata contains hash classifications.
# If these are present, they are used directly to select singlets.
# Otherwise, HTODemux is run on the HTO assay.

if ("hash_classification" %in% colnames(cite_filtered@meta.data)) {
  message("Using hash classifications from metadata.")
  
  hash_table <- as.data.frame(table(cite_filtered$hash_classification))
  colnames(hash_table) <- c("hash_classification", "n_cells")
  
  write.table(
    hash_table,
    file = file.path(output_dir, "hash_classification_table.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  message("No hash_classification column found. Running HTODemux.")
  
  cite_filtered <- HTODemux(
    cite_filtered,
    assay = "HTO",
    positive.quantile = 0.99
  )
  
  hash_table <- as.data.frame(table(cite_filtered$HTO_classification.global))
  colnames(hash_table) <- c("HTO_classification.global", "n_cells")
  
  write.table(
    hash_table,
    file = file.path(output_dir, "HTODemux_classification_table.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# HTO heatmap for visual inspection.
DefaultAssay(cite_filtered) <- "HTO"

pdf(file.path(output_dir, "HTO_heatmap.pdf"), width = 8, height = 6)
print(HTOHeatmap(cite_filtered, assay = "HTO", ncells = min(5000, ncol(cite_filtered))))
dev.off()


################################################################################
# 13. Retain singlets
################################################################################

if ("hash_classification" %in% colnames(cite_filtered@meta.data)) {
  cite_singlet <- subset(cite_filtered, subset = hash_classification == "Singlet")
} else {
  Idents(cite_filtered) <- "HTO_classification.global"
  cite_singlet <- subset(cite_filtered, idents = "Singlet")
}

message("Number of cells before singlet filtering: ", ncol(cite_filtered))
message("Number of singlet cells retained: ", ncol(cite_singlet))


################################################################################
# 14. RNA-based preprocessing, clustering, UMAP, and tSNE
################################################################################

DefaultAssay(cite_singlet) <- "RNA"

cite_singlet <- NormalizeData(
  cite_singlet,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cite_singlet <- FindVariableFeatures(
  cite_singlet,
  selection.method = "vst",
  nfeatures = 2000
)

top_variable_features <- head(VariableFeatures(cite_singlet), 20)

variable_feature_plot <- VariableFeaturePlot(cite_singlet)
variable_feature_plot <- LabelPoints(
  plot = variable_feature_plot,
  points = top_variable_features,
  repel = TRUE
)

save_plot(variable_feature_plot, "RNA_variable_features", width = 8, height = 5)

cite_singlet <- ScaleData(
  cite_singlet,
  features = VariableFeatures(cite_singlet)
)

cite_singlet <- RunPCA(
  cite_singlet,
  features = VariableFeatures(cite_singlet)
)

elbow_plot <- ElbowPlot(cite_singlet, ndims = 30)
save_plot(elbow_plot, "RNA_PCA_elbow_plot", width = 6, height = 5)

n_pcs <- 10

cite_singlet <- FindNeighbors(
  cite_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)

cite_singlet <- FindClusters(
  cite_singlet,
  resolution = 0.6
)

cite_singlet <- RunUMAP(
  cite_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)

cite_singlet <- RunTSNE(
  cite_singlet,
  reduction = "pca",
  dims = 1:n_pcs
)


################################################################################
# 15. Basic RNA-based visualization
################################################################################

umap_clusters <- DimPlot(
  cite_singlet,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("RNA-based Seurat clusters")

save_plot(umap_clusters, "UMAP_RNA_seurat_clusters", width = 7, height = 5)

tsne_clusters <- DimPlot(
  cite_singlet,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("RNA-based Seurat clusters")

save_plot(tsne_clusters, "TSNE_RNA_seurat_clusters", width = 7, height = 5)

metadata_columns_to_plot <- c(
  "sample",
  "hash_classification",
  "detailed_hash_classification",
  "HTO_classification",
  "HTO_classification.global"
)

metadata_columns_to_plot <- metadata_columns_to_plot[
  metadata_columns_to_plot %in% colnames(cite_singlet@meta.data)
]

for (metadata_column in metadata_columns_to_plot) {
  plot <- DimPlot(
    cite_singlet,
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
# 16. ADT marker visualization
################################################################################

DefaultAssay(cite_singlet) <- "ADT"

# Example antibody markers. Only markers present in the ADT assay are plotted.
adt_markers <- c(
  "Ms.CD45",
  "Ms.CD3",
  "Ms.CD4",
  "Ms.CD8a",
  "Ms.CD8b",
  "Ms.CD19",
  "HuMs.CD11b",
  "Ms.CD11c",
  "Ms.F4_80",
  "Ms.Ly.6C",
  "Ms.Ly.6G",
  "Ms.NK.1.1",
  "Ms.CD44",
  "Ms.CD62L",
  "Ms.CD69",
  "Ms.CD25",
  "Ms.CD279",
  "Ms.CD274",
  "Ms.CD55"
)

adt_markers <- adt_markers[adt_markers %in% rownames(cite_singlet[["ADT"]])]

if (length(adt_markers) > 0) {
  adt_feature_plot <- FeaturePlot(
    cite_singlet,
    features = adt_markers,
    reduction = "umap",
    ncol = 3
  )
  
  save_plot(
    adt_feature_plot,
    "UMAP_ADT_markers",
    width = 12,
    height = ceiling(length(adt_markers) / 3) * 3
  )
  
  adt_dot_plot <- DotPlot(
    cite_singlet,
    features = adt_markers,
    group.by = "seurat_clusters",
    assay = "ADT"
  ) +
    RotatedAxis() +
    ggtitle("ADT marker expression by RNA cluster")
  
  save_plot(
    adt_dot_plot,
    "DotPlot_ADT_markers_by_cluster",
    width = 12,
    height = 6
  )
  
  adt_vln_plot <- VlnPlot(
    cite_singlet,
    features = adt_markers[1:min(6, length(adt_markers))],
    group.by = "seurat_clusters",
    pt.size = 0,
    ncol = 3
  )
  
  save_plot(
    adt_vln_plot,
    "VlnPlot_ADT_markers_by_cluster",
    width = 12,
    height = 6
  )
}


################################################################################
# 17. RNA marker gene visualization
################################################################################

DefaultAssay(cite_singlet) <- "RNA"

rna_markers <- c(
  "Epcam", "Krt8", "Krt18",
  "Ptprc", "Cd3d", "Cd3e", "Cd4", "Cd8a",
  "Lyz2", "Chil3",
  "Col1a1", "Pecam1"
)

rna_markers <- rna_markers[rna_markers %in% rownames(cite_singlet)]

if (length(rna_markers) > 0) {
  rna_feature_plot <- FeaturePlot(
    cite_singlet,
    features = rna_markers,
    reduction = "umap",
    ncol = 3
  )
  
  save_plot(
    rna_feature_plot,
    "UMAP_RNA_marker_genes",
    width = 12,
    height = ceiling(length(rna_markers) / 3) * 3
  )
  
  rna_dot_plot <- DotPlot(
    cite_singlet,
    features = rna_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle("RNA marker gene expression by cluster")
  
  save_plot(
    rna_dot_plot,
    "DotPlot_RNA_marker_genes_by_cluster",
    width = 10,
    height = 5
  )
}


################################################################################
# 18. Find RNA cluster markers
################################################################################

DefaultAssay(cite_singlet) <- "RNA"

cluster_markers <- FindAllMarkers(
  cite_singlet,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(
  cluster_markers,
  file = file.path(output_dir, "RNA_cluster_markers_FindAllMarkers.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.table(
  top10_markers,
  file = file.path(output_dir, "RNA_cluster_markers_top10_per_cluster.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


################################################################################
# 19. Exemplary cell type annotation
################################################################################

# This annotation is deliberately simple and should be reviewed manually.
# It can be refined using both RNA marker genes and ADT protein marker expression.

cite_singlet$seurat_clusters <- as.character(cite_singlet$seurat_clusters)

# Example cluster mapping.
# Replace these assignments after inspecting the marker gene and ADT plots.
cluster_to_celltype <- c(
  "0" = "Cancer cells",
  "1" = "Cancer cells",
  "2" = "Cancer cells",
  "3" = "Cancer cells",
  "4" = "Cancer cells",
  "5" = "Cancer cells",
  "6" = "Myeloid cells",
  "7" = "T cells",
  "8" = "B cells",
  "9" = "Fibroblasts",
  "10" = "Endothelial cells"
)

cite_singlet$celltype <- cluster_to_celltype[cite_singlet$seurat_clusters]
cite_singlet$celltype[is.na(cite_singlet$celltype)] <- "Unassigned"

umap_celltype <- DimPlot(
  cite_singlet,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Exemplary cell type annotation")

save_plot(umap_celltype, "UMAP_celltype_annotation", width = 8, height = 5)


################################################################################
# 20. Summarize cell numbers
################################################################################

cell_number_cluster <- as.data.frame(table(cite_singlet$seurat_clusters))
colnames(cell_number_cluster) <- c("cluster", "n_cells")

write.table(
  cell_number_cluster,
  file = file.path(output_dir, "cell_numbers_by_cluster.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cell_number_celltype <- as.data.frame(table(cite_singlet$celltype))
colnames(cell_number_celltype) <- c("celltype", "n_cells")

write.table(
  cell_number_celltype,
  file = file.path(output_dir, "cell_numbers_by_celltype.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if ("sample" %in% colnames(cite_singlet@meta.data)) {
  cell_number_sample <- as.data.frame(table(cite_singlet$sample))
  colnames(cell_number_sample) <- c("sample", "n_cells")
  
  write.table(
    cell_number_sample,
    file = file.path(output_dir, "cell_numbers_by_sample.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if ("hash_classification" %in% colnames(cite_singlet@meta.data)) {
  cell_number_hash <- as.data.frame(table(cite_singlet$hash_classification))
  colnames(cell_number_hash) <- c("hash_classification", "n_cells")
  
  write.table(
    cell_number_hash,
    file = file.path(output_dir, "cell_numbers_by_hash_classification_singlets.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}


################################################################################
# 21. Save processed object and session information
################################################################################

saveRDS(
  cite_singlet,
  file = file.path(output_dir, "VBPN_C_A_subcut_CITEseq_singlet_processed_seurat_object.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("CITE-seq example analysis completed successfully.")
message("Results written to: ", output_dir)