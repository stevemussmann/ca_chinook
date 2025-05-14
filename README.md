# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Dependencies and First-time Setup
1. Create a Conda environment for running snakemake
```
conda create -c conda-forge -c bioconda -c r -n snakemake snakemake r-base r-tidyverse r-remotes
```

2. Clone the `eriqande/california-chinook-microhaps` to your computer. 
```
git clone https://github.com/eriqande/california-chinook-microhaps.git
```

3. Within the `california-chinook-microhaps` folder that you cloned to your computer, make a directory named `data`
```
cd /path/to/california-chinook-microhaps ## Replace `/path/to` in this command with the actual path to the california-chinook-microhaps folder
mkdir data
```


## Steps

1. Before sequencing, make an Illumina sample sheet (see `example_files/SampleSheet.csv`). Use this sample sheet to conduct the sequencing run on your Illumina sequencer.
2. Create a folder within the `california-chinook-microhaps/data` directory that is named according to your sequencing run number (e.g., run001). Within this folder, place the `SampleSheet.csv` file and a folder named `raw` that contains all of your `fastq.gz` files. See example below:

```
.
└── run001
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
