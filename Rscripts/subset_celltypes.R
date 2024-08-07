# Flames Analysis R Script

# 1. Load Libraries  ------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(devtools)
devtools::install_github('immunogenomics/presto')
install.packages("readr")
library(readr)



dest_dir = '/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/subset.R_out'
last_dir = '/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/master.R_out'

# 2. Load object ----------------------------------------------------------
setwd("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/")
seu <- LoadSeuratRds("SR_LR_seu.rds")

# 3. Create subsets of seu object -----------------------------------------

macrophage_subset <- subset(seu, subset = celltypes == "macrophages")
# Verify Subset Created properly
table(seu@meta.data$celltypes)
table(macrophage_subset@meta.data$celltypes)
# Find markers if needed for downstream analysis
cluster1.macrophage <- FindMarkers(macrophage_subset, ident.1 = "1")
head(cluster1.macrophage, n = 20)


# Create T-cell object
Tcell_subset <- subset(seu, subset = celltypes %in% c("CD4_T", "CD8T_NK"))
table(seu@meta.data$celltypes)
table(Tcell_subset@meta.data$celltypes)

# Create tumor-epithelial object
tum_epi_subset <- subset(seu, subset = celltypes %in% c("tumor", "stromal"))
# Verify subset created properly
table(seu@meta.data$celltypes)
table(tum_epi_subset@meta.data$celltypes)


# 3. Include goi --------------------- -------------------------------

# What genes of interest are you looking at? Include them here
iso_of_i <- c("SPP1", "APOC1", "NUPR1", "APOE","FN1","CAPG","CD9","TREM2","LGALS3","FTL","IL1B",
              "GNLY","KRT8","TSPAN31")



# 5. Create plots ---------------------------------------------------------
# Use isoform of interest to generate plots with each subset
# May change features as desired

# Adjust assay to short read
DefaultAssay(macrophage_subset) <- "RNA"
DefaultAssay(Tcell_subset) <- "RNA"
DefaultAssay(tum_epi_subset) <- "RNA"

mac_ioi_feature <- FeaturePlot(macrophage_subset, features = iso_of_i)
mac_ioi_feature
tcell_ioi_feature <- FeaturePlot(Tcell_subset, features = iso_of_i)
tcell_ioi_feature
tum_epi_ioi_feature <- FeaturePlot(tum_epi_subset, features = iso_of_i)
tum_epi_ioi_feature


mac_ioi_heatmap <- DoHeatmap(macrophage_subset, features = iso_of_i)
mac_ioi_heatmap
tcell_ioi_heatmap <- DoHeatmap(Tcell_subset, features = iso_of_i)
tcell_ioi_heatmap
tum_epi_ioi_heatmap <- DoHeatmap(tum_epi_subset, features = iso_of_i)
tum_epi_ioi_heatmap


# Longread plotting -------------------------------------------------------

# Because longread name is longer, we must adjust matching features accordingly
# Extract all feature names from the longread assay
all_features <- rownames(GetAssayData(seu, assay = "longread", slot = "data"))

# Function to match features that start with any of the iso_of_i entries
get_matching_features <- function(features, iso_of_i) {
  matching_features <- unique(unlist(sapply(iso_of_i, function(x) {
    grep(paste0("^", x), features, value = TRUE)
  })))
  return(matching_features)
}

# Get the matching features
matching_features <- get_matching_features(all_features, iso_of_i)

# Since matching_features is returning a bunch, let's narrow it down to 4 for now. Can adjust accordingly
# Limit the number of features to display to 4
limited_features <- matching_features[1:min(4, length(matching_features))]

DefaultAssay(macrophage_subset) <- "longread"
DefaultAssay(Tcell_subset) <- "longread"
DefaultAssay(tum_epi_subset) <- "longread"

LR_mac_ioi_feature <- FeaturePlot(macrophage_subset, features = limited_features)
LR_mac_ioi_feature
LR_tcell_ioi_feature <- FeaturePlot(Tcell_subset, features = limited_features)
LR_tcell_ioi_feature
LR_tum_epi_ioi_feature <- FeaturePlot(tum_epi_subset, features = limited_features)
LR_tum_epi_ioi_feature

LR_mac_ioi_heatmap <- DoHeatmap(macrophage_subset, features= matching_features)
LR_mac_ioi_heatmap
LR_tcell_ioi_heatmap <- DoHeatmap(Tcell_subset, features=matching_features)
LR_tcell_ioi_heatmap
LR_tum_epi_ioi_heatmap <- DoHeatmap(tum_epi_subset, features = matching_features)
LR_tum_epi_ioi_heatmap


# 6. Do plots for goi: ex SPP1 --------------------------------------------

# Extract all feature names from the longread assay
DefaultAssay(macrophage_subset) <- 'longread'
spp1_isoforms <- rownames(macrophage_subset)[startsWith(rownames(macrophage_subset), 'SPP1')]
spp1_isoforms<- spp1_isoforms[order(spp1_isoforms)]
apoc_isoforms <- rownames(macrophage_subset)[startsWith(rownames(macrophage_subset), 'APOC1')]
apoc_isoforms <- apoc_isoforms[order(apoc_isoforms)]

# Generate FeaturePlot and Heatmap for SPP1 and APOC genes
spp1_feature_plot <- FeaturePlot(macrophage_subset, features = spp1_isoforms, reduction ='umap',
                                ncol = 2)

spp1_heatmap <- DoHeatmap(macrophage_subset, features = spp1_isoforms, assay='longread')

apoc_feature_plot <- FeaturePlot(macrophage_subset, features = apoc_isoforms, reduction = 'umap', ncol = 2)
apoc_heatmap <- DoHeatmap(macrophage_subset, features = apoc_isoforms, assay='longread')

# Display the plots
print(spp1_feature_plot)
print(spp1_heatmap)
ggsave(file.path(dest_dir, "LR_spp1_feature_plot.png"), plot =spp1_feature_plot, width = 10, height=15)
ggsave(file.path(dest_dir, "LR_spp1_heatmap.png"), plot =spp1_heatmap, width = 10, height=8)

print(apoc_feature_plot)
print(apoc_heatmap)
ggsave(file.path(dest_dir, "LR_apoc_feature_plot.png"), plot =apoc_feature_plot, width = 10, height=15)
ggsave(file.path(dest_dir, "LR_apoc_heatmap.png"), plot =apoc_heatmap, width = 10, height=8)


# 6. Save objects ---------------------------------------------------------

saveRDS(macrophage_subset, "macrophage.rds")
saveRDS(Tcell_subset, "Tcell.rds")
saveRDS(tum_epi_subset, "tumor_epi.rds")



# 3. Other analysis --------------------- -------------------------------


# Uses iso_of_i to filter through the diff expression markers to find a match
marker_data <- read_csv(file.path(last_dir, "LR_differential_expression_markers.csv"))

# Initialize an empty data frame to store matched rows
all_matched_isos <- data.frame()

# Iterate through each pattern in iso_of_i and find matching rows
for (pattern in iso_of_i) {
  matched_isos <- subset(marker_data, grepl(pattern, gene, ignore.case = TRUE))
  all_matched_isos <- rbind(all_matched_isos, matched_isos)
}

setwd("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/subset.R_out/")
write_csv(all_matched_isos, "isoforms_of_interest.csv")
print(all_matched_isos)

isoforms_of_interest <- read_csv("isoforms_of_interest.csv")

