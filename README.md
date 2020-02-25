# genome_transcriptome_to_gff_file_pipeline

This repository serves the replicability of the manuscript: “Nutrient-driven genome evolution revealed by comparative genomics of chrysomonad flagellates”, which explains the meaning of the working steps.

**Requirements:**
  - Snakemake (see https://snakemake.readthedocs.io/en/stable/)
  - Conda (optional, recommended for building an executable environment, https://docs.conda.io/en/latest/)
    In the terminal type:    
    ```
    conda env create -f conda_env_list.yml
    ```

**Getting Started:**
The directory structure is designed according to the following structure:
```
project folder   
│└───subfolder1: code
│      here should be the scripts
└───subfolder2: raw
│      here should be your input data e.g. sequencing reads
│      optional: subfolder with pacbio/nanopore reads
│      optional: subfolder with transcriptomic data
└───subfolder3: results
       folder will automatically created
```
At the beginning of the main script, the samples and folders can be set. An assignment of different input sources (e.g. illumina and pacbio reads) will be set automatically as long as the file path contain the sample name. 

Before running the script, activate the working environment:
```
    source activate gff_creation_env
```
**Run the script:**
```
    snakemake -s map_RNASeq2genome_2gff_general.py
```
