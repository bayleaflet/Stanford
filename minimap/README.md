# Minimap Processing

This script processes a fastq file from the previous step, flexiplex, maps to the reference genome, then indexes it.

## Prerequisites
Make sure you have the following tools installed:
- `gzip`
Make sure you are in the right environment! I used:
conda activate long_reads

## Usage

To run the script, you need to provide the following arguments in order:
1. `sample`: The sample name.
2. `fastqpath: The path to the fastq file from your previous flexiplex run

### Example Command

```bash

nohup nice bash run_minimap.sh
                P8640_22261B_mp
                /mnt/ix1/Projects_lite/20240603_BC_P8640/scRNA/B01_LongRead/01.0_run_flexiplex/P8640_22261B_merged_flex.fastq.gz
                >out.minimap.tmp
                2>err.minimap.tmp
                &



### Output
Your output will look something like this:
P8640_22261B_mp.primary.st.bam
P8640_22261B_mp.primary.st.bam.bai
