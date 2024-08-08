Author: Baylee Christensen, Intern, Ji Research Group
Date: 8/8/2024

FLAIR PIPELINE

1. Activate the correct environment
```conda activate /venvs2/anaconda3/envs/groupenvs/flair```
2. To complete the flair quantify step, you will need a reads manifest. My name is called reads_mani_ex.tsv for an example.

3. Pass in arguments and run script
```nohup nice bash run_flair.sh -f /your/path/to/fq >out.flr.tmp 2>err.flr.tmp ``

To complete processing of output into a matrix, you should consult.. for the scripts she used to create the Seurat ready matrix. The script I have in here isn't quite correct.
Another slight issue I had when running this script is I had to manually create test.quant.temp for it to complete fully..
