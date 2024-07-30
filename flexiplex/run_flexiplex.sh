#!/bin/bash

# Author: Baylee Christensen - Intern, Ji Research Group
# Creation Date: 7/24/24
# Description: Run flexiplex on your sample using a whitelist of existing barcodes. Unzips the fastq, applied flexiplex using whitelist, merges everything back up.

# Set global variables
flexi_path="/mnt/ix2/Experimental_tools/flexiplex/v1.01/flexiplex"

set_globals() {
    sample=$1
    fastqpath=$2
    barcodes=$3
}

# Function to process the fastq file
# Pass in directory for fastq file

process_fastq() {
    local sample=$1
    local fastqpath=$2
    local barcodes=$3

    # Check if the fastqpath file exists
    if [ ! -f "$fastqpath" ]; then
        echo "Error: File $fastqpath not found!"
        return 1
    fi

    # Run flexiplex on fastq file
    gunzip -c "$fastqpath" | \
    "$flexi_path" -d 10x5v2 -p 8 -k "$barcodes" | \
    gzip > "${sample}_merged_flex.fastq.gz"

    # Check if the processing was successful
    if [ $? -ne 0 ]; then
        echo "Error: Processing failed!"
        return 1
    fi

    echo "Processing completed successfully."
}

# Main function to execute the script
main() {
    if [ "$#" -ne 3 ]; then
        echo "Usage: $0 <sample> <fastqpath> <barcodes>"
        exit 1
    fi

    set_globals "$@"
    process_fastq "$sample" "$fastqpath" "$barcodes"
}

# Call main
main "$@"

