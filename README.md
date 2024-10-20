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
Title: Spatial Transcriptomics Analysis Using Seurat

Description: This R script is designed for comprehensive analysis of spatial transcriptomics data using the Seurat package. The script handles all steps of the analysis pipeline, from data loading to visualization of results.

Key Components of the Script:

Data Loading: The script begins by loading the expression matrix and corresponding metadata from specified file paths. It ensures the data is properly ingested into R for further processing.

Seurat Object Creation: Converts the loaded data into a Seurat object, enabling advanced single-cell analysis techniques. This includes integrating metadata that may contain critical information for subsequent analysis steps.

Quality Control and Filtering: Implements initial quality control by filtering out cells based on the number of detected genes. It then calculates the percentage of mitochondrial DNA as a proxy for cell health, removing cells with unusually high mitochondrial content to ensure data quality.

Data Normalization and Feature Selection: The script normalizes the data using a log normalization method. It identifies highly variable features, which are crucial for distinguishing different cell types or states.

Dimensionality Reduction: Performs principal component analysis (PCA) followed by Uniform Manifold Approximation and Projection (UMAP) to reduce the dimensionality of the data. This step is critical for visualizing complex data structures in a lower-dimensional space.

Clustering: Finds neighbors based on the reduced dimensions and identifies clusters using a resolution parameter that can be adjusted depending on the dataset's complexity.

Data Visualization: Generates various plots to visualize the analysis results, including violin plots for mitochondrial content, PCA elbow plots to determine the number of principal components, UMAP plots to visualize the clusters, and heatmaps to show the top 10 markers per cluster.

Data Saving: Finally, the script saves the processed Seurat object for future use, ensuring that all transformations and results can be revisited and further analyzed.

Purpose of the Script: The purpose of this script is to provide a streamlined workflow for analyzing spatial transcriptomics data, enabling researchers to quickly and effectively identify and visualize patterns of gene expression across different cell types or tissue regions.

Usage Instructions:

Update the file paths to point to your local data files.
Run the script in an R environment such as RStudio.
Adjust parameters such as the resolution for clustering based on the specific characteristics of your data.
Conclusion: This script is a key tool for researchers working with spatial transcriptomics data, providing a detailed yet adaptable approach to understanding the spatial distribution of gene expression. It leverages the powerful features of the Seurat package, ensuring robust analysis and insightful visualizations.

## Contact
For any questions or issues, please open an issue on this repository or contact me directly.
