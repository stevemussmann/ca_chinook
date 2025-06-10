# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Table of Contents
1. [Dependencies and First-time Setup](#installation)
    * [Basic account configuration, conda installation, etc.](#condainstall)
    * [Setup for mega-simple-microhap-snakeflow](#mega)
    * [Setup files from this repository](#myscripts)
    * [Setting up Rstudio](#rstudio)
    * [Pipeline updates](#update)
2. [Running the Snakemake Pipeline](#pipeline)
3. [Processing the Snakemake Pipeline Output](#processing)
4. [Processed Outputs](#output)
5. [Optional: Demultiplex with bcl2fastq](#bcl2fastq)
6. [Running rubias](#rubias)

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
bash Miniconda3-latest-Linux-x86_64.sh
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

2. To make life easier while cloning github repositories, run the following:
```
git config --global http.sslverify false
```

4. Clone the `eriqande/mega-simple-microhap-snakeflow` to your computer. 
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

7. We also want to copy some files to the `mega-simple-microhap-snakeflow` directory so that it doesn't attempt to download and/or build them upon the first run of the pipeline. Copy the entire `resources` folder from the AFTC 'rando' server (found within the `GTseq_processing/CA_Chinook_microhaplotype_files` directory) to `~/local/src/mega-simple-microhap-snakeflow`

8. Within the `mega-simple-microhap-snakeflow` folder that you cloned to your computer, make a directory named `data`. This is the folder where you will put all of your data when running the pipeline. 
```
mkdir ~/local/src/mega-simple-microhap-snakeflow/data
```

### setup files from this repository <a name="myscripts"></a>
1. Clone this repository to the ~/local/scripts folder that you made earlier.
```
cd ~/local/src
git clone https://github.com/stevemussmann/ca_chinook.git
```

2. Make all scripts executable
```
cd ~/local/src/ca_chinook/scripts
chmod u+x *.pl
chmod u+x *.R
chmod u+x *.sh
```

3. Link scripts in your ~/local/bin folder so they can be executed from anywhere on your computer.
```
cd ~/local/bin
for file in ~/local/src/ca_chinook/scripts/*.pl; do ln -s $file; done;
for file in ~/local/src/ca_chinook/scripts/*.R; do ln -s $file; done;
for file in ~/local/src/ca_chinook/scripts/*.sh; do ln -s $file; done;
```


### setting up Rstudio <a name="rstudio"></a>
1. Install R, Rstudio, and Rtools on your computer from apps-to-go. Packages should be named something like `IFW-R 4.4.2`, `IFW-RStudio-2024.09.1`, and `IFW-Rtools44`.

2. Open Rstudio. In Rstudio, install the `devtools` library if you do not already have it.
```
install.packages("devtools", dependencies=TRUE)
```

3. Download microhaplot for Rstudio. Go to [https://github.com/ngthomas/microhaplot](https://github.com/ngthomas/microhaplot). Click on the `<> Code` button, then `Download Zip`. If given an option, download this file to your `Downloads` folder in Windows.
  
4. Now install microhaplot in Rstudio. Replace `username` in the command below with your Windows username. 
```
devtools::install_local('C:\Users\username\Downloads\microhaplot-master.zip', upgrade='never')
```

<hr>

### Pipeline updates <a name="update"></a>
1. To update the ca_chinook scripts to the latest versions, navigate to your `~/local/src/ca_chinook` directory and run `git pull`. This should pull in the latest versions of any scripts used for processing.

```
cd ~/local/src/ca_chinook
git pull
```

<hr>

## Running the Snakemake Pipeline <a name="pipeline"></a>

1. Before sequencing, make an Illumina sample sheet (see `example_files/SampleSheet.csv`). Use this sample sheet to conduct the sequencing run on your Illumina sequencer. Optionally, you can also demultiplex using bcl2fastq if you did not use the sample sheet to automatically demultiplex files on the sequencer. [Go to bcl2fastq instructions.](#bcl2fastq)
2. Create a folder within the `mega-simple-microhap-snakeflow/data` directory that is named according to your sequencing run number (e.g., run001). Within this folder, place the `SampleSheet.csv` file and a folder named `raw` that contains all of your `fastq.gz` files. See example below:

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

5. From the `mega-simple-microhap-snakeflow` directory, run the following command to execute the pipeline. You can adjust the number of cores (currently set at 4) up or down depending upon how many physical processor cores your computer has. This command will take a while, perhaps an hour or more, especially when running the pipeline for the first time. 
```
cd ~/local/src/mega-simple-microhap-snakeflow/
snakemake --config run_dir=data/run001 --configfile config/Chinook/config.yaml --use-conda --cores 4
```

<hr>

## Processing the Snakemake Pipeline Output <a name="processing"></a>
This section explains how to obtain and modify files so that they can be read by some R scripts to produce the final genotypes file. I wrote a script that should handle retrieving the snakemake pipeline output from WSL, copying it to your OneDrive, and automatically handle some of the file filtering and file conversions. For example, if you are working on the first run (run001) of p134, then you can execute this script with the following command:
```
caChinookCopyAndProcess.sh p134 run001
```
If the script works, you can skip to [step 8](#step8). Otherwise, start at step 1. Most procedures in this section can alternatively be accomplished by dragging and dropping files in the Windows interface, but VCF filtering and another file modification require the command line. 

1. It now may be easier to copy some files to a new location to make them easier to access. For example, you could make a folder in your OneDrive to hold these files. I suggest opening your OneDrive folder in WSL and running the following command. Substitute the name of your run for the example `run001` in the first line of the command block below.
```
RUN="run001"
mkdir -p p134/$RUN/snakemake_output p134/$RUN/processing/microhaplot
```
This will create the following directory structure. The `snakemake_output` folder is where you will copy the unmodified outputs of the snakemake pipeline. The `processing` folder is where you will run a series of scripts on select output files to filter and convert the microhaplotype data to a usable form.
```
p134
└── run001
    ├── processing
    │   └── microhaplot
    └── snakemake_output
```

2. Copy three folders from your `mega-simple-microhap-snakeflow/data/$RUN/Chinook` folder to the new location. These are `idxstats`, `microhaplot`, and `vcfs`. 
```
WINDOWSUSER=`powershell.exe '$env:UserName' | sed 's/\r//g'`
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/idxstats/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/.
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/microhaplot/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/.
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/vcfs/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/.
```

3. Go to the directory that contains the .vcf file with ROSA genotypes. This will be in the `p134/$RUN/snakemake_output/vcfs/ROSA/target_fasta/rosawr` folder and will be named `variants-bcftools.vcf`.
```
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/vcfs/ROSA/target_fasta/rosawr
```

4. Filter the ROSA genotypes for sequencing depth and quality in vcftools. This is accomplished by running the `vcftools.sh` script on the ROSA .vcf genotypes file, specifying your input file (e.g., `variants-bcftools.vcf`) and your run number (e.g., `run001`) on the command line. This example will output a file named `CH_run001_greb1_q20dp5.recode.vcf` which will be used as input for the R script. For example:
```
vcftools.sh variants-bcftools.vcf run001
```

5. Copy the `CH_run001_greb1_q20dp5.recode.vcf` file to `p134/$RUN/processing/microhaplot`
```
cp CH_run001_greb1_q20dp5.recode.vcf /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/processing/.
```

6. Get the .csv file that contains read counts for the sex-linked marker. This will be the `ordered-read-counts-table.csv` file in the `p134/$RUN/Chinook/idxstats/target_fastas/ROSA/rosawr` folder. This .csv file should contain a column named `sdy_I183` which has read counts for the sex-linked marker. The number of reads sequenced per individual will be converted into female/male calls. Copy this file to the `p134/$RUN/processing` folder.
```
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/idxstats/target_fastas/ROSA/rosawr
cp ordered-read-counts-table.csv /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/processing/.
```

7. Go into the `p134/$RUN/snakemake_output/microhaplot` folder. Run the `fixMicrohaplot.sh` script from the scripts folder in this repository. This will rename some functions within the `ui.R` and `server.R` scripts in this folder so that they are compatible with the most recent versions of the `ggiraph` R package which is a dependency of the `microhaplot` R package.
```
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/microhaplot
fixMicrohaplot.sh
```

8. <a name="step8"></a>Open Rstudio and then open the modified `server.R` script in Rstudio. This will be in the `p134/run001/snakemake_output/microhaplot` folder in your OneDrive that was created by the `caChinookCopyAndProcess.sh` script. Click the "Run App" button to launch the shiny program. Wait until the program finishes loading, then select `FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds` from the "Select Data Set" dropdown box. Then click the "Table" button (bottom, right of center of window) and finally click the "Download" button (in the top right portion of the window). Save the file with an informative name (e.g., `run001_lfar_wrap_vgll3six6.csv`) in the `p134/$RUN/processing/microhaplot` folder that you made in step 1 of this section. Leave Rstudio and the shiny server open.
   
9. Repeat the previous step to export `FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds` with an informative name (e.g., `run001_rosa_microhap_snplicon.csv`) in the `p134/$RUN/processing/microhaplot` folder that you made in step 1 of this section.

10. Checkpoint: Your `processing` folder should be arranged as you see below. Verify this is accurate before proceeding.
```
processing
├── microhaplot
│   ├── run001_lfar_wrap_vgll3six6.csv
│   └── run001_rosa_microhap_snplicon.csv
├── CH_run001_greb1_q20dp5.recode.vcf
└── ordered-read-counts-table.csv
```

11. Open a new Linux shell from inside your `p134/$RUN/processing` directory. An easy way to do this is to navigate to this directory in Windows Explorer, hold the `shift` key while right-clicking somewhere in the empty space in this window, then choose "Open Linux shell here" from the right-click menu. Activate your snakemake conda environment and run the `microhaplotopia.R` script from this directory. This will convert the files you generated and/or copied in the past several steps into a usable format.
```
conda activate snakemake
microhaplotopia.R -r CH_run001_greb1_q20dp5.recode.vcf
```
Most settings can be left as default. The only thing you should have to specify is the name of the VCF file you created in step 8. However, below I provide a comprehensive list of settings that can be modified (if desired) from the command line. 
```
Usage: ./microhaplotopia.R [options]


Options:
        -a ALLELEBALANCE, --alleleBalance=ALLELEBALANCE
                minimum allele balance threshold (default = 0.35)

        -d HAPDEPTH, --hapDepth=HAPDEPTH
                minimum haplotype sequencing depth (default = 4)

        -D TOTDEPTH, --totDepth=TOTDEPTH
                minimum locus sequencing depth (default = 8)

        -g GREB1ROSAOUT, --greb1rosaOut=GREB1ROSAOUT
                output file for greb1rosa haplotype string (default = greb1rosa_all_hapstr.txt)

        -G GREBINFO, --grebInfo=GREBINFO
                info needed to generate greb1rosa haplotype string (default = ~/local/src/ca_chinook/example_files/greb1_roha_alleles_reordered_wr.txt)

        -f FINALOUT, --finalOut=FINALOUT
                name for final output .csv file (default = haps_2col_final.csv)

        -m MICROHAPLOT, --microhaplot=MICROHAPLOT
                microhaplot directory name (default = microhaplot)

        -o OUTPUT, --output=OUTPUT
                output directory name (default = output)

        -r ROSA, --rosa=ROSA
                ROSA VCF file (No default; required)

        -R REPORTS, --reports=REPORTS
                reports directory name (default = reports)

        -s SDY, --sdy=SDY
                sdy read count file (default = ordered-read-counts-table.csv)

        -x MAXFEMALE, --maxFemale=MAXFEMALE
                maximum sdy reads allowed to call female (default = 2)

        -X MAXUNKNOWN, --maxUnknown=MAXUNKNOWN
                maximum number of sdy reads for unknown (default = 5)

        -h, --help
                Show this help message and exit
```

Outputs from the `microhaplotopia.R` script are discussed in the next section.

12. Convert the locus names in your final genotype .csv file. Change directories into the `processing/output` directory and run the `caChinookRenameLoci.pl` script. This will output a new file with the loci named according to the AmpliconName field in [this table](https://github.com/eriqande/california-chinook-microhaps/blob/main/inputs/Calif-Chinook-Amplicon-Panel-Information.csv). The new file will have `lociRenamed` inserted into its file name before the .csv extension.
```
caChinookRenameLoci.pl -f haps_2col_final.csv
```
The output will be `haps_2col_final.lociRenamed.csv`. 

<hr>

## Processed Outputs <a name="output"></a>
The pipeline writes several files to the `output` and `reports` directories. After running `microhaplotopia.R` you should have the following in your `processing` folder:
```
processing
├── CH_run001_greb1_q20dp5.recode.vcf
├── microhaplot
│   ├── run001_lfar_wrap_vgll3six6.csv
│   └── run001_rosa_microhap_snplicon.csv
├── ordered-read-counts-table.csv
├── output
│   ├── greb1rosa_all_hapstr.txt
│   ├── haps_2col_final.csv
│   └── sdy_calls.csv
└── reports
    ├── extra_alleles.csv
    ├── extra_alleles_individuals.csv
    ├── extra_alleles_locus.csv
    ├── extra_alleles_source.csv
    ├── individual_depth.csv
    ├── locus_depth.csv
    ├── missingData.png
    ├── nxa.csv
    └── nxaCount.csv
```
### Output Directory
This directory contains the main outputs, including the final genotypes file.
1. `greb1rosa_all_hapstr.txt` is a tab-delimited file that contains haplotype strings for the greb1rosa loci. See example below:
```
Indiv   hapstr
4330-001        LNNLLNLLLLLL
4330-002        ENNEENEEEEEH
4330-003        LNNLLNLLLLLL
4330-004        LNNLLNLLLLLL
4330-005        ?NNLLNLLL??L
```

2. `sdy_calls.csv` is a comma-delimited file that contains the genetic sex calls (sdy_sex) for all individuals, and the number of times the sdy marker was identified among reads generated for each individual (sdy_I183).
```
Indiv,sdy_I183,sdy_sex
4330-001,1379,Male
4330-002,6,Male
4330-003,1139,Male
4330-004,2,Female
4330-005,1113,Male
```

3. `haps_2col_final.csv` is a comma-delimited file that contains all genotype data for all individuals. The column percMicroHap reports the percentage of loci that amplified successfully for each individual.
```
Indiv,sdy_sex,hapstr,percMicroHap,NC_037099.1:62937268-62937373_1,NC_037099.1:62937268-62937373_2,NC_037104.1:55923357-55923657_1,NC_037104.1:55923357-55923657_2,...
4330-001,Male,LNNLLNLLLLLL,92.74611398963731,T,T,A,A,...
4330-002,Male,ENNEENEEEEEH,95.85492227979275,T,T,A,A,...
4330-003,Male,LNNLLNLLLLLL,98.96373056994818,C,T,A,A,...
4330-004,Female,LNNLLNLLLLLL,100,T,T,A,A,...
4330-005,Male,?NNLLNLLL??L,87.56476683937824,T,T,A,A,...
```

### Reports Directory
1. `extra_alleles.csv`
2. `extra_alleles_individuals.csv`
3. `extra_alleles_locus.csv`
4. `extra_alleles_source.csv`
5. `individual_depth.csv`
6. `locus_depth.csv`
7. `missingData.png`
8. `nxa.csv`
9. `nxaCount.csv`


<hr>

## Demultiplexing with bcl2fastq (optional) <a name="bcl2fastq"></a>
If you demultiplex with bcl2fastq, do not use the `--no-lane-splitting` option. This is because the [mega-simple-microhap-snakeflow](https://github.com/eriqande/mega-simple-microhap-snakeflow) pipeline expects a fastq file name format of `LibraryName_S1_L001_R1_001.fastq.gz`. Using the `--no-lane-splitting` option removes the 'L001' portion of the name and will cause the pipeline to fail. You may also need to reduce the number of allowed `--barcode-mismatches` if using 6-bp length barcodes.

```
bcl2fastq --barcode-mismatches 0
```

<hr>

## Running rubias <a name="rubias"></a>
### First time setup
With your `snakemake` conda environment active, run the following:
```
R --slave -e "install.packages('rubias', dependencies=TRUE, repos='http://cran.rstudio.com')"
```

### Running the script
Place your baseline and mixture files together in the same directory. Run rubias on these files with the following example command. Replace .csv file names with your actual file names:
```
rubias.R -m mixtureFile.csv -b baselineFile.csv
```

### Outputs
The `output_final` directory will be created in the folder from which you executed the `rubias.R` script. It will contain the following outputs:
1. `all_top_repgroup_sumPofZ.csv`
2. `all_top3pops.csv`
3. `all_toppop.csv`
4. `final_duplicates.csv`
