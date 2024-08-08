Author: Baylee Christensen, Intern, Ji Research Group
Date: 7/17/24
FLAMES PIPELINE


Note: This particular pipeline is ran in two separate steps due to the differing environments being used.


Initiate run_flames.sh:

1. Activate environment:
```conda activate /venvs2/anaconda3/envs/groupenvs/flames```

2. Review the configuration file in first_run, and create your own (if needed).


3. Initiate script
# Keep in mind, you should be using your flexiplex fastq if you are performing scRNA analysis.
``` nohup nice bash run_flames.sh -f /your/path/to/fastq -c /your/path/to/config_file.json 2>stderr.txt >stdout.txt & ```


Outputs of master.sh:
1. transcript_count.csv.gz,
2. Intermediate files not used
-  transcript_count.bad_coverage.csv.gz
-   isoform_FSM_annotation.csv
-   align2genome.bam
-   transcript_assembly.fa
-   realign2transcript.bam

Initiate matrix_seurat_prep.py:

Next, you will need to assemble the output into a usable matrix that can be loaded into Seurat.
1. Activate the correct environment:
``` conda activate python-3 ```


2. Run matrix_seurat_prep.py
```nohup nice python matrix_seurat_prep.py -t /path/to/transcript_count.csv.gz/ -f /path/to/features.tsv.gz 2>err.msp.tmp >out.msp.tmp &```
The transcript count file should come from running run_flames.sh.
Optional: -o for output, otherwise output file will be named matrix.csv
Also, for example, My features path is here: /mnt/ix1/Seq_Runs/20200326_NV1_0834/22214B_gene_expr/outs/filtered_feature_bc_matrix/features.tsv.gz



Lastly, if you are interested in viewing your matrix to ensure the output is what your expected, please use view_matrix.R within R studio and adjust paths to view the matrix.

Further loading of this matrix (into Seurat) can be found in ../C01_Combined/01_flames_analysis
