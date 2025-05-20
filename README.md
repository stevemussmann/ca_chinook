# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Table of Contents
1. [Dependencies and First-time Setup](#installation)
    * [Basic account configuration, conda installation, etc.](#condainstall)
    * [Setup for mega-simple-microhap-snakeflow](#mega)
    * [Setup files from this repository](#myscripts)
    * [Setting up Rstudio](#rstudio)
2. [Running the Pipeline](#pipeline)
3. [Processing the Snakemake Pipeline Output](#processing)
4. [Pipeline Outputs](#output)
5. [Optional: Demultiplex with bcl2fastq](#bcl2fastq)

<hr>

## Dependencies and First-time Setup <a name="installation"></a>

### basic account configuration, conda installation, etc. <a name="condainstall"></a>
Before proceeding, make sure you have Windows Subsystem for Linux (WSL) installed on your computer. The installation code below has been tested in Ubuntu 22.04. The installation process described below (i.e., Conda installation, creation of local/src and local/bin folders, modification of .bashrc, etc.) only needs to be done once for setup of this pipeline unless you move to a new computer. Miniconda needs to be installed if it is not already configured on your computer. If you already have Miniconda installed, then skip to step 5 in this section.

1. Launch Windows Subsystem for Linux (WSL) and download the miniconda installer:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2. Install miniconda accepting all defaults. Answer 'yes' when asked if you want to initialize conda. 

```
bash https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

3. Exit and relaunch WSL before proceeding. When your terminal window reopens, run the following command:
```
conda config --set auto_activate_base false
```

4. Once again, exit and relaunch WSL before proceeding.

5. Create location in which you will install software for genotyping. Then change directories into the `~/local/src` folder that you just made.
```
mkdir -p ~/local/src ~/local/bin
```

6. Edit your bash profile to look for executables in the ~/local/bin folder.
```
echo 'export PATH=$PATH:$HOME/local/bin' >> ~/.bashrc
```

7. Also add one more thing to your .bashrc to make some later steps a little easier.
```
echo "alias wget='wget --no-check-certificate'" >> ~/.bashrc
```

8. Exit and relaunch WSL one last time before proceeding.


### setup mega-simple-microhap-snakeflow <a name="mega"></a>
1. Create a Conda environment for running snakemake
```
conda create -c conda-forge -c bioconda -c r -n snakemake snakemake r-base r-tidyverse r-remotes r-devtools r-optparse vcftools zlib liblzma-devel
```

3. Clone the `eriqande/mega-simple-microhap-snakeflow` to your computer. 
```
cd ~/local/src
git clone https://github.com/eriqande/mega-simple-microhap-snakeflow.git
```

4. Patch the `mega-simple-microhap-snakeflow` pipeline so that it doesn't attempt to install the `microhaplotextract` library every time the pipeline runs (or ever). We will instead install this package manually.
```
cd ~/local/src/mega-simple-microhap-snakeflow/workflow/script
rm create_microhaplot_folder.R
wget https://raw.githubusercontent.com/stevemussmann/ca_chinook/refs/heads/main/patches/create_microhaplot_folder.R
```

5. Download the `microhaplotextract` R package (this is the 'just-for-extracting' branch of the microhaplot github repository), activate the snakemake Conda environment that you made in Step 1, and install the `microhaplotextract` from the source zip file.
```
cd ~/local/src/mega-simple-microhap-snakeflow/
wget https://github.com/eriqande/microhaplot/archive/just-for-extracting.zip
conda activate snakemake
R --slave -e "devtools::install_local('just-for-extracting.zip', upgrade='never')"
```

6. We're also going to install the 'emc-edits' branch of the `microhaplotopia` R package in a similar way.
```
wget https://github.com/eriqande/microhaplotopia/archive/emc-edits.zip
R --slave -e "devtools::install_local('emc-edits.zip', upgrade='never')"
```

7. We also want to copy some files to the `mega-simple-microhap-snakeflow` directory so that it doesn't attempt to download and/or build them upon the first run of the pipeline. Copy the entire `resources` folder from the AFTC 'rando' server to `~/local/src/mega-simple-microhap-snakeflow`

8. Within the `mega-simple-microhap-snakeflow` folder that you cloned to your computer, make a directory named `data`. This is the folder where you will 
```
mkdir ~/local/src/mega-simple-microhap-snakeflow/data
```

### setup files from this repository <a name="myscripts"></a>


### setting up Rstudio <a name="rstudio"></a>
1. Install R and Rstudio on your computer from apps-to-go.

2. Install microhaplot in Rstudio

<hr>

## Running the Pipeline <a name="pipeline"></a>

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

<hr>

## Processing the Snakemake Pipeline Output <a name="processing"></a>
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

<hr>

## Pipeline Outputs <a name="output"></a>

<hr>

## Demultiplexing with bcl2fastq (optional) <a name="bcl2fastq"></a>
If you demultiplex with bcl2fastq, do not use the `--no-lane-splitting` option. This is because the [mega-simple-microhap-snakeflow](https://github.com/eriqande/mega-simple-microhap-snakeflow) pipeline expects a fastq file name format of `LibraryName_S1_L001_R1_001.fastq.gz`. Using the `--no-lane-splitting` option removes the 'L001' portion of the name and will cause the pipeline to fail. You may also need to reduce the number of allowed `--barcode-mismatches` if using 6-bp length barcodes.

```
bcl2fastq --barcode-mismatches 0
```
