#!/bin/bash
# Run in /venvs2/anaconda3/envs/groupenvs/flames venv
# Run after running flexiplex

# Default paths
GENOME_FA=/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/fasta/genome.fa
GENES_GTF=/mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/genes/genes.gtf
MINIMAP2_DIR=/mnt/ix1/Resources/tools/minimap2/v2.24
FLAMES_PY=/mnt/ix1/Resources/tools/FLAMES/v1.11.0/python
CWD=$(pwd)

# Function to display usage
usage() {
    echo "Usage: $0 -f FASTQPATH"
    exit 1
}

# Parse command-line arguments
while getopts ":f:" opt; do
    case ${opt} in
        f )
            FASTQPATH=$OPTARG
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

# Check if FASTQPATH was provided
if [ -z "$FASTQPATH" ]; then
    echo "FASTQPATH is required."
    usage
fi

# Create output directory
OUTDIR=${CWD}/out
mkdir -p $OUTDIR

# Run the FLAMES pipeline
python3 ${FLAMES_PY}/sc_long_pipeline.py -i $FASTQPATH -a $GENES_GTF \
    --config_file ${CWD}/config2_file.json --outdir ${CWD} --genomefa $GENOME_FA --minimap2_dir $MINIMAP2_DIR
