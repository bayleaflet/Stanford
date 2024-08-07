# Script for creating Seurat object to analyze short read and long reads, generating cluster UMAP
# Save rds frequently


# Skipping creation of short read plots, as SR object is already complete.
# Skip until we get to long read analysis

# 1. Load Libraries  ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
devtools::install_github('immunogenomics/presto')

# 2. Load data ------------------------------------------------------------
## 2.1 Load object, starting with short reads -----------------------------
seu.data <- Read10X(data.dir = "/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/A01_Seurat/P8640_22261B/")
# Initialize the Seurat object with the raw (non-normalized data).
seu <- CreateSeuratObject(counts = seu.data, project = "P8640", min.cells = 3, min.features = 200)


# 3. Quality Control ------------------------------------------------------

# Filter out mitochondrial cells, calculating cell swith high mt counts
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# nCountRNA is total transcripts in cell, nFeatureRNA is total unique transcripts
seu
head(seu@meta.data, 5)

## 3.1 Visualize QC metrics with violin plot --------------------
vln_plot <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("SR_vln_plot.png", plot =vln_plot,width = 10, height=8)

## 3.2 Visualize QC metrics with scatter plot -------------------
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_plots <- plot1 | plot2
ggsave("SR_scatter_plot.png", plot =combined_plots,width = 10, height=8)



## 3.3 Required: Complete Filtration --------------------------------------
# Completes filtration, where we can change the number fields based on our decided criteria
# Small sample for P8640
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
seu



# 4. Normalize Data -------------------------------------------------------
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
seu



# 5. Identify Highly Variable Features ------------------------------------

SR_top10 <- head(VariableFeatures(seu), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = SR_top10, repel = TRUE)
combined_plots <- plot1 | plot2
ggsave("SR_variable_plots.png", plot=combined_plots, width=12, height=8)



# 6. Scale the data -------------------------------------------------------
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu


# 7. Linear Dimensional Reduction -----------------------------------------
# npcs = number of principle components, set at 20
# automatically when RunPCA runs, the reduction name is pca
seu <- RunPCA(seu, features = VariableFeatures(object = seu), npcs = 20)
print(seu[["pca"]], dims = 1:5, nfeatures = 5)

## 7.1 Plot Generation -------------------------------------------
VizDim <- VizDimLoadings(seu, dims = 1:2, reduction = "pca", nfeatures=20)
ggsave("SR_visdim_plot.png", plot =VizDim,width = 10, height=8)
DimPlot <- DimPlot(seu, reduction = "pca") + NoLegend()
ggsave("SR_dim_plot.png", plot =DimPlot,width = 10, height=8)
# Look at cells, and put in for below
seu
HeatMap1 <- DimHeatmap(seu, dims = 1, cells = 238, balanced = TRUE)
ggsave("SR_HeatMap1_plot.png", plot =HeatMap1, width = 10, height=8)
HeatMap2 <- DimHeatmap(seu, dims = 1:12, cells = 238, balanced = TRUE)
ggsave("SR_HeatMap2_plot.png", plot =HeatMap2, width = 10, height=8)
elbow <- ElbowPlot(seu)
ggsave("SR_elbow_plot.png", plot =elbow, width = 10, height=8)






# 8. Cluster Cells --------------------------------------------------------
# Review Elbow plot to adjust dims (adjusted to 20, re-review)
seu <- FindNeighbors(seu, dims = 1:20, reduction = "pca")
seu <- FindClusters(seu, resolution = 1.0)
head(Idents(seu), 5)
seu



# 9. Non - Linear Dimensional Reduction -----------------------------------
# review dims
# Automatic nape is umap
seu <- RunUMAP(seu, dims = 1:20, reduction='pca')
seu

# See the plot without labels
umap_plot<-DimPlot(seu, reduction = "umap")
ggsave("SR_umap_reduction.png", plot =umap_plot, width = 10, height=8)

# SR_seu.rds no longer a usuable copy due to errors, and more being overwritten
# saveRDS(seu, file = "SR_seu.rds")
# Final filetype after LR analysis will be SR_LR_seu.rds
# Short read analysis only is named: SR_seu_backup.rds

saveRDS(seu, file = "SR_seu_backup.rds")
# Keep saved copy separate

seu <- LoadSeuratRds("SR_seu_backup.rds")
seu
getwd()



# 10. Find All Markers: Downstream analysis for dge ----------------------
seu.markers <- FindAllMarkers(seu, only.pos = TRUE)
seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
FeaturePlot(seu, features = c(top10$gene))
seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seu, features = top10$gene) + NoLegend()
# Mostly focusing on dot plot and not dge


## 10.1 Generate Dot Plot -------------------------------------------------


#Anuja: please save differential expression
# Generate Plot with genes of interest
goi <-   c("CD68", "CD14","FCGR3A","CLEC9A","CD1C", "HLA-DQA1", "HLA-DQB1","KIT", "TPSB2",
           "CD3D", "IL7R", "CD8A", "NKG7", "GNLY", "CD79A", "MS4A1", "JCHAIN", "IGHG1",
           "EPCAM", "TFF3", "MUC2",
           "PECAM1", "PLVAP","VWF","COL1A1","LUM", "DCN")
theme_dot_2 <- ggplot2::theme(
  axis.text = element_text(size = 14, face = "bold"),
  axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1),
  axis.title = element_blank(),
  plot.title = element_text(hjust = 0.5),
  legend.text = element_text(size = 14, face = "bold"),
  axis.line = element_line(color = "black"),
  legend.position = "top",
  legend.box = "horizontal")

pdf("dotplot_clusters.pdf", width=20, height=10)
DotPlot(seu, features = goi)+
  geom_point(
    aes(size = pct.exp),
    shape = 21,colour = "black",
    stroke = 0.5) +
  scale_colour_gradient2(low ="#5385BC", mid = "#FFFFFF", high = "#A50026") + coord_flip() +
  guides(size = guide_legend(override.aes = list(
    shape = 21,
    colour = "black",
    fill = "white"
  ))) + theme_dot_2
dev.off()


## 10.2 Annotate Clusters -------------------------------------------------
head(seu@meta.data)
DimPlot(seu, reduction="umap", group.by = "RNA_snn_res.1")

table(seu@meta.data$RNA_snn_res.1, seu@meta.data$RNA_snn_res.1)

clusters <- levels(as.factor(seu@meta.data$RNA_snn_res.1))
clusters
cell_types <- c("CD4_T", "macrophages", "tumor", "CD8T_NK", "stromal")
seu@meta.data$celltypes <- plyr::mapvalues(x = seu@meta.data$RNA_snn_res.1,
                                                  from = clusters,
                                                  to = cell_types)
head(seu@meta.data)
table(seu@meta.data$RNA_snn_res.1, seu@meta.data$celltypes)

Idents(seu) <- seu@meta.data$celltypes

#Anuja:
#c0:CD4_T
#c1:macrophages
#c2:tumor
#c3:CD8T_NK
#c4:stromal

## 10.3 Save new plot with labels, save object ----------------------------
umap_plot<-DimPlot(seu, reduction = "umap", label=TRUE)
ggsave("SR_umap_reduction.png", plot =umap_plot, width = 10, height=8)
saveRDS(seu, file = "SR_seu_backup.rds")


# -------------------------------------------------------------------------
#==========================================================================#
# Longread Assay Initiation -----------------------------------------------
#==========================================================================#
# -------------------------------------------------------------------------


seu <- LoadSeuratRds("SR_seu_backup.rds")
setwd("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/")
dest_dir='/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/master.R_out/'



# 11.1 Extract barcodes from Seu obj and LR
short_barcodes <- colnames(seu)

# Load the sparse matrix (assuming it's saved as a Matrix Market file)
longread_matrix <- readMM("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/P8640_22261B.matrix.mtx")

# Load the genes and cell barcodes
genes <- readLines("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/P8640_22261B.features.tsv")
cells <- readLines("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/P8640_22261B.barcodes.tsv")
# Add a '-1' to every cell
cells <- paste0(cells, "-1")

# Set the row and column names for the matrix
rownames(longread_matrix) <- genes
colnames(longread_matrix) <- cells
# Convert to a data frame
longread_matrix <- as.data.frame(as.matrix(longread_matrix))
longread_barcodes <- colnames(longread_matrix)



seu


# 12. Filters LR matrix to only include barcodes present in seu obj -------
matched_barcodes <- intersect(longread_barcodes, short_barcodes)
longread_matrix <- longread_matrix[, matched_barcodes]

# 13. Create LR assay and add to Seu obj ----------------------------------
longread_assay <- CreateAssayObject(counts=longread_matrix)
seu[["longread"]] <-longread_assay

seu
head(seu@meta.data)
DefaultAssay(seu) <- "longread"
head(seu@meta.data)

# 14. Quality Control -----------------------------------------------------
# Different quality control. measures required for longread analysis, no mt for example

## 14.1 Visualize current QC metrics with violin plot ---------------------
LR_vln_plot <- VlnPlot(seu, features = c("nFeature_longread", "nCount_longread"), ncol = 2)
# adjusted y axis for comparison of other pipelines
p1 <- LR_vln_plot[[1]] + scale_y_continuous(limits = c(0, 8000))
p2 <- LR_vln_plot[[2]] + scale_y_continuous(limits = c(0, 160000))
LR_formatted_plots <- p1 | p2
ggsave(file.path(dest_dir, "LR_vln_plot.png"), plot =LR_vln_plot, width = 10, height=8)

## 14.2 Visualize QC metrics with scatter plot ----------------------------
LR_plot1 <- FeatureScatter(seu, feature1 = "nCount_longread", feature2 = "nFeature_longread")
ggsave(file.path(dest_dir,"LR_scatter_plot.png"), plot =LR_plot1, width = 10, height=8)


seu

## 14.3 Perform Filtration ------------------------------------------------
# Completes filtration, where we can change the number fields based on our decided criteria
# Be very cautious with usage of filtration with this step for long reads. First view the scatter plot to see if this is really needed.
# seu <- subset(seu, subset = nFeature_longread > 200 & nFeature_longread < 5000)
seu

# 14.4 Normalize the new assay --------------------------------------------
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 5000)
seu

# 15. Identify Highly Variable Features -----------------------------------

LR_top10 <- head(VariableFeatures(seu), 10)

# Plot variable features
LR_plot1 <- VariableFeaturePlot(seu)
LR_plot2 <- LabelPoints(plot = LR_plot1, points = LR_top10, repel = TRUE)
LR_combined_plots <- LR_plot1 | LR_plot2
ggsave(file.path(dest_dir,"LR_variable_plots.png"), plot=LR_combined_plots, width=12, height=8)


# 16. Scale the data ------------------------------------------------------
all.LR.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.LR.genes)
seu


# 17. Linear Dimensional Reduction ----------------------------------------
# npcs = number of principle components, set at 20
seu <- RunPCA(seu, features = VariableFeatures(object = seu), npcs = 20, reduction.name = "lr_pca", reduction.key = "lr_pca")
seu
print(seu[["lr_pca"]], dims = 1:5, nfeatures = 5)

## 17.1 Plot Generation ------------------------------------------
LR_VizDim <- VizDimLoadings(seu, dims = 1:2, reduction = "lr_pca", nfeatures=20)
ggsave(file.path(dest_dir,"LR_visdim_plot.png"), plot =LR_VizDim,width = 10, height=8)
LR_DimPlot <- DimPlot(seu, reduction = "lr_pca") + NoLegend()
ggsave(file.path(dest_dir,"LR_dim_plot.png"), plot =LR_DimPlot,width = 10, height=8)
# Look at cells, and put in for below
seu
png(file.path(dest_dir,"LR_HeatMap1_plot.png"), width = 14, height = 8, units = "in", res = 200)
DimHeatmap(seu, dims = 1, cells = 238, balanced = TRUE, reduction = "lr_pca")
dev.off()
png(file.path(dest_dir,"LR_HeatMap2_plot.png"), width = 14, height = 8, units = "in", res = 200)
DimHeatmap(seu, dims = 1:12, cells = 238, balanced = TRUE, reduction = "lr_pca")
dev.off()
LR_elbow <- ElbowPlot(seu, ndims=20, reduction="lr_pca")
ggsave(file.path(dest_dir,"LR_elbow_plot.png"), plot =LR_elbow, width = 10, height=8)

# 18. Cluster Cells -------------------------------------------------------
# Review Elbow plot to adjust dims (adjusted to 20, re-review)
# Somewhere in 18 or 19 the active.indent item is overwritten in these from SR to LR
seu <- FindNeighbors(seu, dims = 1:20, reduction = 'lr_pca')
seu <- FindClusters(seu, resolution = 1.0)
head(Idents(seu), 5)
seu


# 19. Non - Linear Dimensional Reduction ----------------------------------
# Always review dimensions, reduction, and reduction name
seu <- RunUMAP(seu, dims = 1:20, reduction.name = "lr_umap", reduction = "lr_pca")
seu
LR_umap_plot<-DimPlot(seu, reduction = "lr_umap")
ggsave(file.path(dest_dir,"LR_umap_reduction.png"), plot =LR_umap_plot, width = 10, height=8)

saveRDS(seu, file = "FLR_SR_LR_seu.rds")
# Keep saved copy separate
seu


# 20. Verify Data ---------------------------------------------------------

head(seu@meta.data)
# Write out comparison table of SR to LR with celltypes
compare_table <- table(seu@meta.data$celltypes, seu@meta.data$longread_snn_res.1)
compare_table
write.csv(compare_table, file.path(dest_dir,'compare.csv'))
# Write out comparison table of SR to LR without celltypes
compare2_table <- table(seu@meta.data$RNA_snn_res.1, seu@meta.data$longread_snn_res.1)
compare2_table
write.csv(compare2_table, file.path(dest_dir,'compare_no_ann.csv'))


# 21. Create merged UMAP --------------------------------------------------

sr_plot <- DimPlot(seu, group.by = "celltypes",
                   label= TRUE,
                   reduction = "umap" ) +
  ggtitle("RNA_snn_res.1")

lr_plot <- DimPlot(seu,
                   group.by = "celltypes",
                   label = TRUE,
                   reduction = "lr_umap") +
  ggtitle("longread_snn_res.1")

# Set custom axis limits
lr_plot_custom <- lr_plot +
  xlim(-15, 5) +
  ylim(-10, 10)
new_merged_plot <- sr_plot | lr_plot

ggsave(file.path(dest_dir,"FLR_LR_UMAP.png"), plot = new_merged_plot,width = 14, height=8)



# 22. Find all markers ----------------------------------------------------
seu
markers <- FindAllMarkers(seu)
write.csv(markers, "LR_differential_expression_markers.csv")


# 23. Generate a feature plot ---------------------------------------------
feature_plot <- DotPlot(seu, assay = "RNA", group.by = "RNA_snn_res.1", features = c("CD3D", "TFF3", "CD68", "ACTA2", "DCN"))
ggsave(file.path(dest_dir,"Feature_dot_plot.png"), plot = feature_plot, width = 16, height=8)

# 24. Save seu object -----------------------------------------------------

saveRDS(seu, file = "FLR_SR_LR_seu.rds")
``



