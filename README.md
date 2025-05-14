# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Dependencies and Setup
1. Create a Conda environment for running snakemake
```
conda create -c conda-forge -c bioconda -c r -n snakemake snakemake r-base r-tidyverse r-remotes
```

## Steps

1. Before sequencing, make an Illumina sample sheet (see `example_files/SampleSheet.csv`). Use this sample sheet to conduct the sequencing run on your Illumina sequencer.
2. Create a directory structure that contains `SampleSheet.csv` and a directory named `raw` that contains all of your `fastq.gz` files. See example below:

```
.
├── SampleSheet.csv
└── raw
    ├── CH-15342_S1_L001_R1_001.fastq.gz
    ├── CH-15342_S1_L001_R2_001.fastq.gz
    ├── CH-15343_S2_L001_R1_001.fastq.gz
    ├── CH-15343_S2_L001_R2_001.fastq.gz
    ├── CH-15344_S3_L001_R1_001.fastq.gz
    └── CH-15344_S3_L001_R2_001.fastq.gz
```

3. Run the `preprocess.R` script to create the `samples.csv` and `units.csv` required by the snakemake pipeline
```
preprocess.R -f SampleSheet.csv
```

4. 
