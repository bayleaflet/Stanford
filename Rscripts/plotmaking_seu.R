
# Standardization for comparison
# Load all objects
flames_seu <- LoadSeuratRds("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/SR_LR_seu.rds")
custom_seu <- LoadSeuratRds('/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/02_custom_analysis/CA_SR_LR_seu.rds')
flair_seu <- LoadSeuratRds('/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/FLR_SR_LR_seu.rds')


setwd('/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/comparisons/')
goi <- c("SPP1", "APOC1", "NUPR1", "APOE", "FN1", "CAPG",
         "CD9", "TREM2", "LGALS3", "FTL", "IL1B", "TSPAN31", "GNLY")
# Construct regular expressions to match exact gene prefixes
goi_patterns <- paste0("^", goi, "-")

# Get the expression matrix for the pipeline ------------------------------

flames_expr_matrix <- GetAssayData(flames_seu, slot = "data")
custom_expr_matrix <- GetAssayData(custom_seu, slot = "data")
flair_expr_matrix <- GetAssayData(flair_seu, slot = "data")


# Initialize a data frame to store the results ----------------------------


flames_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)
custom_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)
flair_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)


# Calculate metrics for each filtered isoform, flames ---------------------

for (isoform in rownames(flames_expr_matrix)) {
  # Check if the isoform matches any exact gene prefix from goi
  if (any(sapply(goi_patterns, function(pattern) grepl(pattern, isoform)))) {
    flames_isoform_data <- flames_expr_matrix[isoform, ]
    num_cells <- sum(flames_isoform_data > 0)
    total_transcripts <- sum(flames_isoform_data)
    avg_transcripts_per_cell <- if (num_cells > 0) total_transcripts / num_cells else 0

    # Append to isoform_metrics
    flames_isoform_metrics <- rbind(flames_isoform_metrics, data.frame(
      Isoform = isoform,
      NumCellsExpressingIsoform = num_cells,
      TotalTranscripts = total_transcripts,
      AvgTranscriptsPerCell = avg_transcripts_per_cell,
      stringsAsFactors = FALSE
    ))
  }
}
# Calculate metrics for each filtered isoform, custom ---------------------
# Custom
for (isoform in rownames(custom_expr_matrix)) {
  # Check if the isoform matches any exact gene prefix from goi
  if (any(sapply(goi_patterns, function(pattern) grepl(pattern, isoform)))) {
    custom_isoform_data <- custom_expr_matrix[isoform, ]
    num_cells <- sum(custom_isoform_data > 0)
    total_transcripts <- sum(custom_isoform_data)
    avg_transcripts_per_cell <- if (num_cells > 0) total_transcripts / num_cells else 0

    # Append to isoform_metrics
    custom_isoform_metrics <- rbind(custom_isoform_metrics, data.frame(
      Isoform = isoform,
      NumCellsExpressingIsoform = num_cells,
      TotalTranscripts = total_transcripts,
      AvgTranscriptsPerCell = avg_transcripts_per_cell,
      stringsAsFactors = FALSE
    ))
  }
}

# Calculate metrics for each filtered isoform, flair ---------------------
for (isoform in rownames(flair_expr_matrix)) {
  # Check if the isoform matches any exact gene prefix from goi
  if (any(sapply(goi_patterns, function(pattern) grepl(pattern, isoform)))) {
    flair_isoform_data <- flair_expr_matrix[isoform, ]
    num_cells <- sum(flair_isoform_data > 0)
    total_transcripts <- sum(flair_isoform_data)
    avg_transcripts_per_cell <- if (num_cells > 0) total_transcripts / num_cells else 0

    # Append to isoform_metrics
    flair_isoform_metrics <- rbind(flair_isoform_metrics, data.frame(
      Isoform = isoform,
      NumCellsExpressingIsoform = num_cells,
      TotalTranscripts = total_transcripts,
      AvgTranscriptsPerCell = avg_transcripts_per_cell,
      stringsAsFactors = FALSE
    ))
  }
}


# Sort the data frame by Isoform name in ascending order ------------------
flames_isoform_metrics <- flames_isoform_metrics[order(flames_isoform_metrics$Isoform), ]
custom_isoform_metrics <- custom_isoform_metrics[order(custom_isoform_metrics$Isoform), ]
flair_isoform_metrics <- flair_isoform_metrics[order(flair_isoform_metrics$Isoform), ]


# Write the sorted data frame to a CSV file -------------------------------
write.csv(flames_isoform_metrics, file = "flames_isoform_metrics_goi_sorted_by_isoform.csv", row.names = FALSE)
write.csv(custom_isoform_metrics, file = "custom_isoform_metrics_goi_sorted_by_isoform.csv", row.names = FALSE)
write.csv(flair_isoform_metrics, file = "flair_isoform_metrics_goi_sorted_by_isoform.csv", row.names = FALSE)


# Plotmaking --------------------------------------------------------------


# Create graph
# Assuming isoform_metrics has columns: Isoform, NumCellsExpressingIsoform, TotalTranscripts, AvgTranscriptsPerCell

# Extract gene names (assuming isoform names are like "GENE-ENST...")
flames_isoform_metrics$Gene <- sapply(strsplit(flames_isoform_metrics$Isoform, "-"), `[`, 1)
custom_isoform_metrics$Gene <- sapply(strsplit(custom_isoform_metrics$Isoform, "-"), `[`, 1)
flair_isoform_metrics$Gene <- sapply(strsplit(flair_isoform_metrics$Isoform, "-"), `[`, 1)

# Categorize AvgTranscriptsPerCell into ranges
flames_isoform_metrics$AvgTranscriptsPerCellRange <- cut(flames_isoform_metrics$AvgTranscriptsPerCell,
                                                  breaks = c(-Inf, 1, 2, 3, Inf),
                                                  labels = c("0-1", "1-2", "2-3", "3+"))
custom_isoform_metrics$AvgTranscriptsPerCellRange <- cut(custom_isoform_metrics$AvgTranscriptsPerCell,
                                                         breaks = c(-Inf, 1, 2, 3, Inf),
                                                         labels = c("0-1", "1-2", "2-3", "3+"))
flair_isoform_metrics$AvgTranscriptsPerCellRange <- cut(flair_isoform_metrics$AvgTranscriptsPerCell,
                                                         breaks = c(-Inf, 1, 2, 3, Inf),
                                                         labels = c("0-1", "1-2", "2-3", "3+"))

# Dot Plot: Avg Transcripts per Cell
fl1 <- ggplot(flames_isoform_metrics, aes(x = Gene, y = Isoform)) +
  geom_point(aes(size = as.numeric(AvgTranscriptsPerCellRange), color = NumCellsExpressingIsoform)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Flames Average Transcripts per Cell per Isoform", x = "Gene", y = "Isoform", size = "Avg Transcripts/Cell Range", color = "Num. Cells Expressing Isoform") +
  scale_size_continuous(labels = c("0-1", "1-2", "2-3", "3+")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("flames_avg_transcripts_per_cell.png", width = 10, height = 8)

# Dot Plot: Avg Transcripts per Cell
custom1<-ggplot(custom_isoform_metrics, aes(x = Gene, y = Isoform)) +
  geom_point(aes(size = as.numeric(AvgTranscriptsPerCellRange), color = NumCellsExpressingIsoform)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Custom Average Transcripts per Cell per Isoform", x = "Gene", y = "Isoform", size = "Avg Transcripts/Cell Range", color = "Num. Cells Expressing Isoform") +
  scale_size_continuous(labels = c("0-1", "1-2", "2-3", "3+")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("custom_avg_transcripts_per_cell.png", width = 10, height = 8)

# Dot Plot: Avg Transcripts per Cell
flr1<- ggplot(flair_isoform_metrics, aes(x = Gene, y = Isoform)) +
  geom_point(aes(size = as.numeric(AvgTranscriptsPerCellRange), color = NumCellsExpressingIsoform)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Flair Average Transcripts per Cell per Isoform", x = "Gene", y = "Isoform", size = "Avg Transcripts/Cell Range", color = "Num. Cells Expressing Isoform") +
  scale_size_continuous(labels = c("0-1", "1-2", "2-3", "3+")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("flair_avg_transcripts_per_cell.png", width = 10, height = 8)

avg_trscr_plot <- fl1 | custom1 | flr1
ggsave("combined_avg_trsc.pdf", width = 30, height = 15)



# 3 Dimensional plotting for comparing pipelines --------------------------

library(dplyr)
library(plotly)

# Merge the datasets
merged_data <- flames_isoform_metrics %>%
  inner_join(custom_isoform_metrics, by = c("Gene", "Isoform"), suffix = c("_flames", "_custom")) %>%
  inner_join(flair_isoform_metrics, by = c("Gene", "Isoform")) %>%
  rename(AvgTranscriptsPerCell_flair = AvgTranscriptsPerCell, NumCellsExpressingIsoform_flair = NumCellsExpressingIsoform)

# Create a 3D scatter plot
p <- plot_ly(merged_data,
             x = ~AvgTranscriptsPerCell_flames,
             y = ~AvgTranscriptsPerCell_custom,
             z = ~AvgTranscriptsPerCell_flair,
             text = ~paste("Gene:", Gene, "<br>Isoform:", Isoform),
             color = ~NumCellsExpressingIsoform_flames,
             colors = colorRamp(c("blue", "red")),
             size = ~NumCellsExpressingIsoform_flames,
             marker = list(sizemode = 'diameter'),
             type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Plot: Avg Transcripts per Cell",
         scene = list(xaxis = list(title = 'FLAMES'),
                      yaxis = list(title = 'Custom'),
                      zaxis = list(title = 'FLAIR')))

# Save the plot as an HTML file
htmlwidgets::saveWidget(p, "3d_avg_transcripts_per_cell.html")


## Unfiltered Isoform Metrics

## Unfiltered Isoform Metrics ----------------------------------------------


# Initialize a data frame to store the results for all isoforms


# Flames ------------------------------------------------------------------
flames_unfiltered_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)

# Calculate metrics for each isoform without filtering for goi
for (isoform in rownames(flames_expr_matrix)) {
  flamesu_isoform_data <- flames_expr_matrix[isoform, ]
  flamesu_num_cells <- sum(flamesu_isoform_data > 0)
  flamesu_total_transcripts <- sum(flamesu_isoform_data)
  flamesu_avg_transcripts_per_cell <- if (flamesu_num_cells > 0) flamesu_total_transcripts / flamesu_num_cells else 0

  # Append to unfiltered_isoform_metrics
  flames_unfiltered_isoform_metrics <- rbind(flames_unfiltered_isoform_metrics, data.frame(
    Isoform = isoform,
    NumCellsExpressingIsoform = flamesu_num_cells,
    TotalTranscripts = flamesu_total_transcripts,
    AvgTranscriptsPerCell = flamesu_avg_transcripts_per_cell,
    stringsAsFactors = FALSE
  ))
}

# Sort the data frame by Isoform name in ascending order
flames_unfiltered_isoform_metrics <- flames_unfiltered_isoform_metrics[order(flames_unfiltered_isoform_metrics$Isoform), ]



# Custom ------------------------------------------------------------------
# Initialize a data frame to store the results for all isoforms
custom_unfiltered_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)

# Calculate metrics for each isoform without filtering for goi
for (isoform in rownames(custom_expr_matrix)) {
  customu_isoform_data <- custom_expr_matrix[isoform, ]
  customu_num_cells <- sum(customu_isoform_data > 0)
  customu_total_transcripts <- sum(customu_isoform_data)
  customu_avg_transcripts_per_cell <- if (customu_num_cells > 0) customu_total_transcripts / customu_num_cells else 0

  # Append to unfiltered_isoform_metrics
  custom_unfiltered_isoform_metrics <- rbind(custom_unfiltered_isoform_metrics, data.frame(
    Isoform = isoform,
    NumCellsExpressingIsoform = customu_num_cells,
    TotalTranscripts = customu_total_transcripts,
    AvgTranscriptsPerCell = customu_avg_transcripts_per_cell,
    stringsAsFactors = FALSE
  ))
}
# Sort the data frame by Isoform name in ascending order
custom_unfiltered_isoform_metrics <- custom_unfiltered_isoform_metrics[order(custom_unfiltered_isoform_metrics$Isoform), ]


# Flair -------------------------------------------------------------------
# Initialize a data frame to store the results for all isoforms
flair_unfiltered_isoform_metrics <- data.frame(
  Isoform = character(),
  NumCellsExpressingIsoform = integer(),
  TotalTranscripts = integer(),
  AvgTranscriptsPerCell = numeric(),
  stringsAsFactors = FALSE
)

# Calculate metrics for each isoform without filtering for goi
for (isoform in rownames(flair_expr_matrix)) {
  flairu_isoform_data <- flair_expr_matrix[isoform, ]
  flairu_num_cells <- sum(flairu_isoform_data > 0)
  flairu_total_transcripts <- sum(flairu_isoform_data)
  flairu_avg_transcripts_per_cell <- if (flairu_num_cells > 0) flairu_total_transcripts / flairu_num_cells else 0

  # Append to unfiltered_isoform_metrics
  flair_unfiltered_isoform_metrics <- rbind(flair_unfiltered_isoform_metrics, data.frame(
    Isoform = isoform,
    NumCellsExpressingIsoform = flairu_num_cells,
    TotalTranscripts = flairu_total_transcripts,
    AvgTranscriptsPerCell = flairu_avg_transcripts_per_cell,
    stringsAsFactors = FALSE
  ))
}

# Sort the data frame by Isoform name in ascending order
flair_unfiltered_isoform_metrics <- flair_unfiltered_isoform_metrics[order(flair_unfiltered_isoform_metrics$Isoform), ]



# Continue with plotmaking ------------------------------------------------



# Define a function to classify isoforms as known or novel
classify_isoform <- function(isoform) {
  if (grepl("ENST", isoform)) {
    return("Known")
  } else {
    return("Novel")
  }
}

# Apply classification to each isoform and create a new column
flames_unfiltered_isoform_metrics$IsoformType <- sapply(flames_unfiltered_isoform_metrics$Isoform, classify_isoform)
custom_unfiltered_isoform_metrics$IsoformType <- sapply(custom_unfiltered_isoform_metrics$Isoform, classify_isoform)
flair_unfiltered_isoform_metrics$IsoformType <- sapply(flair_unfiltered_isoform_metrics$Isoform, classify_isoform)

# Count the total number of each type
flames_isoform_counts <- table(flames_unfiltered_isoform_metrics$IsoformType)
custom_isoform_counts <- table(custom_unfiltered_isoform_metrics$IsoformType)
flair_isoform_counts <- table(flair_unfiltered_isoform_metrics$IsoformType)

# Convert the counts to a data frame for plotting
flames_isoform_counts_df <- as.data.frame(flames_isoform_counts)
colnames(flames_isoform_counts_df) <- c("IsoformType", "Count")
custom_isoform_counts_df <- as.data.frame(custom_isoform_counts)
colnames(custom_isoform_counts_df) <- c("IsoformType", "Count")
flair_isoform_counts_df <- as.data.frame(flair_isoform_counts)
colnames(flair_isoform_counts_df) <- c("IsoformType", "Count")



# # Plotting Count results --- --------------------------------------

flm_count<- ggplot(flames_isoform_counts_df, aes(x = IsoformType, y = Count, fill = IsoformType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5,
            fontface = "bold") +
  labs(title = "Total Isoforms Detected with Flames", x = "Isoform Type", y = "Count") +
  coord_cartesian(ylim = c(0, 40000)) + # Set y-axis limit
  theme_minimal()

ggsave("flames_isoform_counts_plot.png", plot = flm_count, width = 10, height = 6)

custom_count <- ggplot(custom_isoform_counts_df, aes(x = IsoformType, y = Count, fill = IsoformType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5,
            fontface = "bold") +
  labs(title = "Total Isoforms Detected with Custom", x = "Isoform Type", y = "Count") +
  coord_cartesian(ylim = c(0, 40000)) + # Set y-axis limit
  theme_minimal()

ggsave("custom_isoform_counts_plot.png", plot = custom_count, width = 10, height = 6)

flr_count<-ggplot(flair_isoform_counts_df, aes(x = IsoformType, y = Count, fill = IsoformType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 5,
            fontface = "bold") +
  labs(title = "Total Isoforms Detected with Flair", x = "Isoform Type", y = "Count") +
  coord_cartesian(ylim = c(0, 40000)) + # Set y-axis limit
  theme_minimal()

ggsave("flair_isoform_counts_plot.png", plot = flr_count, width = 10, height = 6)

# ALl three
combined_counts <- flm_count | custom_count | flr_count
ggsave("Count_comparison.pdf", plot=combined_counts, width=16, height = 10)

# Plotting Quality control -----------------------
flames_vln_plot <- VlnPlot(flames_seu, features = c("nFeature_longread", "nCount_longread"), ncol = 2, group.by = "celltypes")
# adjusted y axis for comparison of other pipelines
flames_p1 <- flames_vln_plot[[1]] + scale_y_continuous(limits = c(0, 8000))
flames_p2 <- flames_vln_plot[[2]] + scale_y_continuous(limits = c(0, 160000))
flames_formatted_plots <- flames_p1 | flames_p2


custom_vln_plot <- VlnPlot(custom_seu, features = c("nFeature_longread", "nCount_longread"), ncol = 2, group.by = "celltypes")
# adjusted y axis for comparison of other pipelines
custom_p1 <- custom_vln_plot[[1]] + scale_y_continuous(limits = c(0, 8000))
custom_p2 <- custom_vln_plot[[2]] + scale_y_continuous(limits = c(0, 160000))
custom_formatted_plots <- custom_p1 | custom_p2


flair_vln_plot <- VlnPlot(flair_seu, features = c("nFeature_longread", "nCount_longread"), ncol = 2, group.by = "celltypes")
# adjusted y axis for comparison of other pipelines
flair_p1 <- flair_vln_plot[[1]] + scale_y_continuous(limits = c(0, 8000))
flair_p2 <- flair_vln_plot[[2]] + scale_y_continuous(limits = c(0, 160000))
flair_formatted_plots <- flair_p1 | flair_p2


merged_qc <- flames_formatted_plots | custom_formatted_plots | flair_formatted_plots
ggsave("QC_Comparison.pdf", plot = merged_qc, width = 24, height = 8)


# Umap Plots --------------------------------------------------------------

# Flames
flm_lr_plot <- DimPlot(flames_seu,
                   group.by = "celltypes",
                   label = TRUE,
                   reduction = "lr_umap") +

  ggtitle("Flames longread_snn_res.1")

# Set custom axis limits
flm_lr_plot <- flm_lr_plot +
  xlim(-15, 7) +
  ylim(-10, 10)

# Custom
cu_lr_plot <- DimPlot(custom_seu,
                   group.by = "celltypes",
                   label = TRUE,
                   reduction = "lr_umap") +
  ggtitle("Custom longread_snn_res.1")

# Set custom axis limits
cu_lr_plot <- cu_lr_plot +
  xlim(-15, 7) +
  ylim(-10, 10)

# Flair

flr_lr_plot <- DimPlot(flair_seu,
                   group.by = "celltypes",
                   label = TRUE,
                   reduction = "lr_umap") +
  ggtitle("Flair longread_snn_res.1")

# Set custom axis limits
flr_lr_plot <- flr_lr_plot +
  xlim(-15, 7) +
  ylim(-10, 10)

merged_plot <- flm_lr_plot | cu_lr_plot | flr_lr_plot

ggsave("Umap_Comparison.pdf", width = 20, height = 8)



# Compare Subsetted objects -----------------------------------------------
flames_seu_sub <- LoadSeuratRds("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/subset.R_out/macrophage.rds")
custom_seu_sub <- LoadSeuratRds('/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/02_custom_analysis/subset.R_out/macrophage.rds')
flair_seu_sub <- LoadSeuratRds('/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/03_flair_analysis/macrophage.rds')


iso_of_i <- c("SPP1", "APOC1", "NUPR1", "APOE","FN1","CAPG","CD9","TREM2","LGALS3","FTL","IL1B",
              "GNLY","KRT8","TSPAN31")

flames_all_features <- rownames(GetAssayData(flames_seu_sub, assay = "longread", layer  = "data"))
custom_all_features <- rownames(GetAssayData(flames_seu_sub, assay = "longread", layer  = "data"))
flair_all_features <- rownames(GetAssayData(flames_seu_sub, assay = "longread", layer  = "data"))

# Function to match features that start with any of the iso_of_i entries
get_matching_features <- function(features, iso_of_i) {
  matching_features <- unique(unlist(sapply(iso_of_i, function(x) {
    grep(paste0("^", x), features, value = TRUE)
  })))
  return(matching_features)
}

# Get the matching features
flames_matching_features <- get_matching_features(flames_all_features, iso_of_i)
custom_matching_features <- get_matching_features(custom_all_features, iso_of_i)
flair_matching_features <- get_matching_features(flair_all_features, iso_of_i)

# Since matching_features is returning a bunch, let's narrow it down to 4 for now. Can adjust accordingly
# Limit the number of features to display to 4
flames_limited_features <- flames_matching_features[1:min(4, length(flames_matching_features))]
custom_limited_features <- custom_matching_features[1:min(4, length(custom_matching_features))]
flair_limited_features <- flair_matching_features[1:min(4, length(flair_matching_features))]

DefaultAssay(flames_seu_sub) <- "longread"
DefaultAssay(custom_seu_sub) <- "longread"
DefaultAssay(flair_seu_sub) <- "longread"



# This section isn't returning the properly scaled data. Needs revisiting.
flames_mac_ioi_feature <- FeaturePlot(flames_seu_sub, features = flames_limited_features)
custom_mac_ioi_feature <- FeaturePlot(custom_seu_sub, features = custom_limited_features)
flair_mac_ioi_feature <- FeaturePlot(flair_seu_sub, features = flair_limited_features)



merged <- flames_mac_ioi_feature | custom_mac_ioi_feature | flair_mac_ioi_feature
ggsave("spp1_Comparison.pdf", plot=merged, width = 20, height = 10)

