# List of gene name prefixes to filter
gene_prefixes <- c("SPP1", "APOC1", "NUPR1", "APOE", "FN1", "CAPG",
                   "CD9", "TREM2", "LGALS3", "FTL", "IL1B", "KRT8", "TSPAN31", "GNLY")

# Function to read GTF file
read_gtf <- function(file) {
  gtf <- read.table(file, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="")
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  return(gtf)
}

# Function to filter GTF file based on gene name prefixes
filter_gtf <- function(gtf, prefixes) {
  # Extract gene_name from attribute column
  gtf$gene_name <- sub(".*gene_name \"([^\"]+)\".*", "\\1", gtf$attribute)

  # Filter rows where gene_name starts with any of the prefixes
  filtered_gtf <- gtf[sapply(gtf$gene_name, function(x) any(sapply(prefixes, function(prefix) startsWith(x, prefix)))), ]

  return(filtered_gtf)
}

# Function to write filtered GTF to a new file
write_gtf <- function(gtf, file) {
  write.table(gtf[,1:9], file=file, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# Main script
input_file <- "/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/genes/genes.gtf"
output_file <- "/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/iso_of_i.gtf"

gtf <- read_gtf(input_file)
filtered_gtf <- filter_gtf(gtf, gene_prefixes)
write_gtf(filtered_gtf, output_file)

cat("Filtering completed. Filtered GTF file saved as", output_file, "\n")
