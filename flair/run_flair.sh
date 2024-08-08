#!/bin/bash

# Author: Baylee Christensen, Intern, Ji Research Group
# Date: 8/8/2024

# Default paths
GENOME=/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38

# Function to display usage
usage() {
    echo "Usage: $0 -f FQFILE"
    exit 1
}

# Parse command-line arguments
while getopts ":f:" opt; do
    case ${opt} in
        f )
            FQFILE=$OPTARG
            ;;
        \? )
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        : )
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check if FQFILE was provided
if [ -z "$FQFILE" ]; then
    echo "FQFILE is required."
    usage
fi

# Script start time
date +"%Y-%m-%d %H:%M:%S"

# FLAIR align
flair align --threads 32 --reads $FQFILE --genome ${GENOME}/fasta/genome.fa \
            --nvrna --junction_bed ${GENOME}/genes/gene_transcripts.bed
echo "aligner complete"
date +"%Y-%m-%d %H:%M:%S"

# FLAIR correct
flair correct --threads 12 -q flair.aligned.bed --gtf ${GENOME}/genes/genes.gtf --genome ${GENOME}/fasta/genome.fa \
              --nvrna
echo "flair correct complete"
date +"%Y-%m-%d %H:%M:%S"

# FLAIR collapse
flair collapse --threads 16 -q flair_all_corrected.bed --reads $FQFILE  \
               --genome ${GENOME}/fasta/genome.fa --gtf ${GENOME}/genes/genes.gtf --annotation_reliant generate \
               --generate_map --trust_ends --isoformtss --no_gtf_end_adjustment --keep_intermediate
echo "flair collapse complete"
date +"%Y-%m-%d %H:%M:%S"

# FLAIR quantify
flair quantify --isoforms flair.collapse.isoforms.fa --reads_manifest reads_mani_ex.tsv --threads 12
echo "flair quantify complete"
date +"%Y-%m-%d %H:%M:%S
