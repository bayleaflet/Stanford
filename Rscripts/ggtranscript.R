# Run from scrna-4a environment
# conda activate scrna-4a


library(ggtranscript)
library(magrittr)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(stringr)
library(GenomicRanges)

setwd("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis")

# For faster processing of the script, utilize filter_gtf.R to filter gtf for genes of interest
# Import custom annotated file
gtf <- import("//mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/subset.R_out/iso_of_i.gtf")
gtf_df <- as.data.frame(gtf)
class(gtf)

tib_gtf <- gtf %>% dplyr::as_tibble()
class(tib_gtf)


# Generate plot for SPP1 --------------------------------------------------


# Filter gtf for gene of interest, here "SPP1"
goi <- "SPP1"
spp1_annotation_from_gtf <- tib_gtf %>%
  filter(!is.na(gene_name), gene_name == goi) %>%
  select(seqnames, start, end, strand, type, gene_name, transcript_name, transcript_id, transcript_biotype)

# Extract exons
spp1_exons <- spp1_annotation_from_gtf %>% filter(type == "exon")

# Read isoforms_of_interest.csv created previously in subset_celltypes.R
isoforms_of_interest <- read.csv("isoforms_of_interest.csv")

# Extract transcript_id values from the gene column
isoform_transcript_ids <- isoforms_of_interest$gene %>%
  str_extract("ENST\\d+")

# Filter spp1_exons based on transcript_id
filtered_spp1_exons <- spp1_exons %>%
  filter(transcript_id %in% isoform_transcript_ids)

# Plot
spp1_plot <- filtered_spp1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype),
    size = 1  # Adjust size if needed
  ) +
  geom_intron(
    data = to_intron(filtered_spp1_exons, "transcript_name"),
    aes(strand = strand)
  ) +
  theme_minimal() +
  labs(title = "Exons of SPP1",
       x = "Genomic Position",
       y = "Transcript")
setwd("/mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/C01_Combined/01_flames_analysis/ggranscript_out/")
ggsave("spp1_exons.png", plot = spp1_plot, height = 10, width = 8)


# Repeat Process for APOC1 ------------------------------------------------

# Now switching to generic names. When using this script, only adjust goi and final plot name.
# Keep in mind, if you do not change the goi or plot name, the plot will be overwritten!


goi <- "APOC1"
default_annotation_from_gtf <- tib_gtf %>%
  filter(!is.na(gene_name), gene_name == goi) %>%
  select(seqnames, start, end, strand, type, gene_name, transcript_name, transcript_id, transcript_biotype)

# Extract exons
default_exons <- default_annotation_from_gtf %>% filter(type == "exon")

# Only need to load isoforms file one time. If it hasn't been loaded yet in previous steps, do this step
# Read isoforms_of_interest.csv created previously in subset_celltypes.R
# isoforms_of_interest <- read.csv("isoforms_of_interest.csv")

# Extract transcript_id values from the gene column
isoform_transcript_ids <- isoforms_of_interest$gene %>%
  str_extract("ENST\\d+")

# Filter default_exons based on transcript_id
filtered_default_exons <- default_exons %>%
  filter(transcript_id %in% isoform_transcript_ids)

# Plot
default_plot <- filtered_default_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = transcript_biotype),
    size = 1  # Adjust size if needed
  ) +
  geom_intron(
    data = to_intron(filtered_default_exons, "transcript_name"),
    aes(strand = strand)
  ) +
  theme_minimal() +
  labs(title = "Exons of APOC1",
       x = "Genomic Position",
       y = "Transcript")

ggsave("APOC1_exons.png", plot = default_plot, height = 10, width = 8)


