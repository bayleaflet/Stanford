#!/bin/bash

# Author: Baylee Christensen, Intern, Ji Research Group
# Date: 7/04/2024
# Incomplete script for running scNanoGPS. Ideally this would run from start to finish but there is issues with the curation step. Never got to the reporter step.

sample=''
scnano=''
working=''
SAMTOOLSPATH=''
# Cloned github repo to ${working}

##################################
#            NANOQC              #
##################################

# Command for running NanoQC to compute histogram of raw read lengths of all reads.

python3 ${scnano}/other_utils/read_length_profiler.py -i ${working}/${sample}/${sample}.fastq.gz -d ${working}/${sample}/
echo "Quality check completed"

# OUTPUT: read_length.png,  Time elapsed 53:51

# Optional: polyAtail is either in first or last 100nt range of each nanopore read. This checks qualities of both, and draws the per-base quality score boxplot.
# Output is named first_tail.fastq.gz and last_tail.fastq.gz, check data with fastqc

# python3 ${working}/scNanoGPS/other_utils/prepare_read_qc.py -i ${working}/${sample}.fastq.gz -d ${working}/${sample}/
#
#
##################################
#           SCANNER              #
##################################

# Command that scans for both TruSeq Read1 and polyA tail of the reads.
# Extracts CBs and UMIs that are neighbored by these sequence blocks and outputs two files: fastq holding insert sequences and barcode_list.tsv storing reads info like names,cbs, umis, etc.

python3 ${scnano}/scanner.py -i ${working}/${sample}/${sample}.fastq.gz -d ${working}/${sample}/ -t 48 --lUMI=10
echo "Scanner completed"
# Ran on 48 nodes on gari, Getting warning message: /venvs2/anaconda3/envs/scnanogps/lib/python3.8/site-packages/Bio/pairwise2.py:278: BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module. warnings.warn(

# OUTPUT:${working}/${sample}/scanner.log.txt, where it details a time elapsed 20 59 :56.20, and more, other files generated:  barcode_list.tsv.

##################################
#          ASSIGNER              #
##################################

# Designed for CB collapsing and estimation of optimal CB number without 'guidance'.
# Can modify the forced_no number if we know the number length of the cell barcode

python3 ${scnano}/assigner.py -t 48 -i ${working}/${sample}/barcode_list.tsv.gz -d ${working}/${sample}
echo "Assigner completed"
# OUTPUT: ${working}/${sample}/notforce/assigner.log.txt. Time Elapsed: 0 : 8 : 37.4, and generated files have CB prefix

##################################
#           CURATOR              #
##################################

# Step used for demultiplexing, filtering, ref genome mapping and re-mapping, and UMI collapsing.

python3 ${scnano}/curator.py -t 48 -d ${working}/${sample} --ref_genome /mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/fasta/genome.fa
echo "Curator completed"

# OUTPUT: curator.log.txt

##################################
#           REPORTER             #
##################################

# Reporter scripts that generates multi-omics profiles from same single cells with long read sequencing data.

# Single cell gene expression profile
# python3 ${working}/scNanoGPS/reporter_expression.py -t 24 --gtf ${working}/scNanoGPS/genes.gtf --featurecounts /mnt/ix1/Resources/tools \
# -d ${working}/${sample} --min_gene_no 250

# Single cell Isoform Profile
# python3 ${working}/scNanoGPS/reporter_isoform.py -t 16 --liqa_ref ${working}/scNanoGPS/genes.liqa.refgene -d ${working}/${sample}/ --gtf ${working}/scNanoGPS/genes.gtf
#
#
# Single cell SNV profile
# python3 ${working}/scNanoGPS/reporter_SNV.py -t 24 --ref_genome /mnt/ix1/Resources/10X_Resources/cellranger_rna/refdata-cellranger-3.0.0/GRCh38/fasta/genome.fa \
# --annovar ${working}/scNanoGPS/annovar/ --annovar_db ${working}/scNanoGPS/annovar/humandb/ -d ${working}/${sample}/

# Generate final summary table
