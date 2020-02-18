# Repository of reproducible code for DOT1L project

This repository is a collection of jupyter notebooks and custom python scripts used in the analysis presented in the paper: "DOT1L Methyltransferase Activity Preserves SOX2-Enhancer Accessibility And Prevents Activation of Repressed Genes In Murine Stem Cells".

## Setting up the environment for reproducing the results

1. Recreate the conda environments in your local machine. You find the requirements in the conda_environment folder;
2. Add the folder AUXILIARY_SCRIPTS to both your PATH variable and to your PYTHONPATH variable;
3. Install snakePipes v.1.1.1 (to do that, follow the detailed explanation available at https://snakepipes.readthedocs.io/en/latest/content/setting_up.html )_

## Workflow 

1. Download the fastq files from GEO: GSE135318
2. Run the appropriate snakepipes pipelines to recreate the standard output folder structure.
3. After changing the paths according to your local file system, run the notebooks.


