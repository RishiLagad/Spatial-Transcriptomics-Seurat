# Spatial-Transcriptomics-Seurat
# Spatial Transcriptomics Analysis

This repository contains scripts and data for analyzing spatial transcriptomics data using the Seurat package in R.

## Files
- `Spatial_Transcriptomics_Analysis.R`:
- # Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load the expression matrix and metadata directly from specified paths
expr_matrix <- read.delim("C:/Users/rushi/Downloads/13322860/PDAC_Full_Merged_Matrix_2024-08-08.txt", 
                          header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim("C:/Users/rushi/Downloads/13322860/PDAC_Merged_Meta_2024-08-08.txt", 
                       header = TRUE, row.names = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "PDAC_Analysis")
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Preliminary filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Calculate the percentage of mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize mitochondrial content
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on mitochondrial content
seurat_obj <- subset(seurat_obj, subset = percent.mt < 5)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Visualization of PCA results
ElbowPlot(seurat_obj)

# Run UMAP for dimensionality reduction
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Find neighbors and identify clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Plot UMAP to visualize clusters
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
print(umap_plot)

# Heatmap of top 10 markers per cluster
top10 <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

# Save the Seurat object for future use
saveRDS(seurat_obj, file = "C:/Users/rushi/Downloads/13322860/seurat_object.rds")
git clone https://github.com/RishiLagad/Spatial-Transcriptomics-Seurat.git



## How to Use
Describe how someone can use these scripts to perform their own analysis.

## Contact
For any questions or issues, please open an issue on this repository or contact me directly.
