# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Dependencies and First-time Setup
1. Create a Conda environment for running snakemake
```
conda create -c conda-forge -c bioconda -c r -n snakemake snakemake r-base r-tidyverse r-remotes r-devtools r-optparse vcftools
```

2. Clone the `eriqande/mega-simple-microhap-snakeflow` to your computer. 
```
git clone https://github.com/eriqande/mega-simple-microhap-snakeflow.git
```

3. Within the `mega-simple-microhap-snakeflow` folder that you cloned to your computer, make a directory named `data`
```
cd /path/to/mega-simple-microhap-snakeflow ## Replace `/path/to` in this command with the actual path to the mega-simple-microhap-snakeflow folder

mkdir data
```

4. Install Rstudio on your computer.


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

## Processing the Snakemake Pipeline Output
1. When finished, the pipeline will output microhaplotypes to the `mega-simple-microhap-snakeflow/data/run001/Chinook/microhaplot` folder. In this folder, run the `fixMicrohaplot.sh` script from the scripts folder in this repository. This will rename some functions within the `ui.R` and `server.R` scripts in this folder so that they are compatible with the most recent versions of the `ggiraph` R package.
```
cd /path/to/mega-simple-microhap-snakeflow/data/run001/Chinook/microhaplot
fixMicrohaplot.sh
```

2. Open Rstudio and then open the `server.R` script in Rstudio. Click the "Run App" button to launch the shiny program. Wait until the program finishes loading, then select `FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds` from the "Select Data Set" dropdown box. Then click the "Table" button (bottom, right of center of window) and finally click the "Download" button (in the top right portion of the window). Save the file with an informative name (e.g., `run001_rosa_microhap_snplicon.csv`). Leave Rstudio and the shiny server open.

3. Repeat the previous step to export `FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds` with an informative name (e.g., `run001_lfar_wrap_vgll3six6.csv`)

4. Get the .csv file from the `mega-simple-microhap-snakeflow` output that contains read counts for the sex-linked marker. This will be the `ordered-read-counts-table.csv` file in the `mega-simple-microhap-snakeflow/data/run001/Chinook/idxstats/target_fastas/ROSA/rosawr` folder. This .csv file should contain a column named `sdy_I183` which has read counts for the sex-linked marker. The number of reads sequenced per individual will be converted into female/male calls.

5. Get the .vcf file that contains the ROSA genotypes from the `mega-simple-microhap-snakeflow` output. This will be in the `mega-simple-microhap-snakeflow/data/run001/Chinook/vcfs/ROSA/target_fasta/rosawr` folder and will be named `variants-bcftools.vcf`.

6. Filter the ROSA genotypes for sequencing depth and quality in vcftools. This is accomplished by running the `vcftools.sh` script on the ROSA .vcf genotypes file, specifying your input file (e.g., `variants-bcftools.vcf`) and your run number (e.g., `run001`) on the command line. For example:
```
./vcftools.sh variants-bcftools.vcf run001
```
This example will output a file named `CH_run001_greb1_q20dp5.vcf` which will be used as input for the R script.

7. 

## Demultiplexing with bcl2fastq (optional) <a name="bcl2fastq"></a>
If you demultiplex with bcl2fastq, do not use the `--no-lane-splitting` option. This is because the [mega-simple-microhap-snakeflow](https://github.com/eriqande/mega-simple-microhap-snakeflow) pipeline expects a fastq file name format of `LibraryName_S1_L001_R1_001.fastq.gz`. Using the `--no-lane-splitting` option removes the 'L001' portion of the name and will cause the pipeline to fail. You may also need to reduce the number of allowed `--barcode-mismatches` if using 6-bp length barcodes.

```
bcl2fastq --barcode-mismatches 0
```
