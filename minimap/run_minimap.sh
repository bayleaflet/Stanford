#!/bin/bash


# Author: Baylee Christensen - Intern, Ji Research Group
# Creation Date: 7/24/24
# Description: Maps longreads to reference genome
# Environment: conda activate long_reads

MINIMAP_EXE='/mnt/ix1/Resources/tools/minimap2/v2.24/minimap2'
SAMTOOLS_EXE='/mnt/ix1/Resources/tools/samtools/v1.15.1/bin/samtools'
GENOME_REF='/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/fasta/genome.fa'
#If doing isoform analysis, provide bed12 annotations for gene transcripts for better alignment around splice junctions
TRANSCRIPT_BED='/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/genes/gene_transcripts.bed'


align_and_process() {
    local sample_name="$1"
    local fq_path="$2"
    echo "Start date and time: "
    date +"%Y-%m-%d %H:%M:%S"
    # Extract filename without extension
    local fq_name=$(basename "$fq_path" .fastq.gz)

    # Minimap2 alignment -> sam | remove supplementary alignments and convert to bam
    $MINIMAP_EXE -ax splice -k14 -t 40 --secondary=no --junc-bed $TRANSCRIPT_BED $GENOME_REF "$fq_path" \
          | $SAMTOOLS_EXE view -Sbh -F 2048 - > "${sample_name}.primary.bam"

    # Sort/index bam
    $SAMTOOLS_EXE sort "${sample_name}.primary.bam" > "${sample_name}.primary.st.bam"
    $SAMTOOLS_EXE index "${sample_name}.primary.st.bam"

    # Cleanup
    rm "${sample_name}.primary.bam"
    echo "End date and time: "
    date +"%Y-%m-%d %H:%M:%S"
}


