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


## Running the Pipeline

1. Before sequencing, make an Illumina sample sheet (see `example_files/SampleSheet.csv`). Use this sample sheet to conduct the sequencing run on your Illumina sequencer. Optionally, you can also demultiplex using bcl2fastq if you did not use the sample sheet to automatically demultiplex files on the sequencer. [Go to bcl2fastq instructions.](#bcl2fastq)
2. Create a folder within the `california-chinook-microhaps/data` directory that is named according to your sequencing run number (e.g., run001). Within this folder, place the `SampleSheet.csv` file and a folder named `raw` that contains all of your `fastq.gz` files. See example below:

```
.
└── run001/
    ├── SampleSheet.csv
    └── raw/
        ├── CH-15342_S1_L001_R1_001.fastq.gz
        ├── CH-15342_S1_L001_R2_001.fastq.gz
        ├── CH-15343_S2_L001_R1_001.fastq.gz
        ├── CH-15343_S2_L001_R2_001.fastq.gz
        ├── CH-15344_S3_L001_R1_001.fastq.gz
        ├── CH-15344_S3_L001_R2_001.fastq.gz
        └── etc.
```

3. Activate your conda environment
```
conda activate snakemake
```

3. Run the `preprocess.R` script to create the two files required by the snakemake pipeline: `samples.csv` and `units.csv`. Watch carefully for errors, and correct any problems until you successfully create the two files.
```
preprocess.R -f SampleSheet.csv
```

4. This pipeline reads from / writes to many files simultaneously. Check the number of files that can be open simultaneously on your computer with the `ulimit` command.
```
ulimit -n
```
If this reports a low number (e.g., 256), then set it to something sufficiently high (e.g., 4096) with the following command:
```
ulimit -n 4096
```

5. From the `california-chinook-microhaps` directory, run the following command to execute the pipeline. You can adjust the number of cores (currently set at 4) up or down depending upon how many physical processor cores your computer has. This can take a while, especially when running the pipeline for the first time. 
```
cd /path/to/california-chinook-microhaps

snakemake --config run_dir=data/run001 --configfile config/Chinook/config.yaml --use-conda --cores 4
```

6. When finished, the pipeline will output microhaplotypes to the `california-chinook-microhaps/data/run001/Chinook/microhaplot` folder.


## Demultiplexing with bcl2fastq (optional) <a name="bcl2fastq"></a>

