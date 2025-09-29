# California Chinook Salmon Genotyping
Scripts and documentation for California Chinook microhaplotypes

## Table of Contents

### Genotyping
1. [Dependencies and First-time Setup](#installation)
    * [Basic account configuration, conda installation, etc](#condainstall)
    * [OPTIONAL: Installing bcl2fastq through conda](#bcl2fastqinstall)
    * [Setup for mega-simple-microhap-snakeflow](#mega)
    * [Setup files from this repository](#myscripts)
    * [Pipeline updates](#update)
3. [Running the Snakemake Pipeline](#pipeline)
4. [Processing the Snakemake Pipeline Output](#processing)
5. [Processed Outputs](#output)
6. [Optional: Demultiplex with bcl2fastq](#bcl2fastq)

### Genetic Stock ID
1. [Installing rubias](#installrubias)
2. [Running rubias](#rubias)
3. [Rubias outputs](#rubiasout)

### CKMRsim
1. [Install CKMRsim](#installCKMRsim)

### Colony
1. [Install Colony](#installColony)

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
   
### OPTIONAL: Installing bcl2fastq through conda <a name="bcl2fastqinstall"></a>

1. Create a new conda environment and install bcl2fastq into it.
```
conda create -n bcl2fastq -c conda-forge -c bioconda -c dranew bcl2fastq=2.19.0
```

2. You should now be able to launch the conda environment with the following command:
```
conda activate bcl2fastq
```

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

<hr>

### Pipeline updates <a name="update"></a>
1. To update the ca_chinook scripts to the latest versions, navigate to your `~/local/src/ca_chinook` directory, run `git pull`, verify that all scripts are executable, and make links in the `~/local/bin` directory. This should pull in the latest versions of any scripts used for processing.

```
cd ~/local/src/ca_chinook
git pull
cd ~/local/src/ca_chinook/scripts

chmod u+x *.pl
chmod u+x *.R
chmod u+x *.sh

cd ~/local/bin
for file in ~/local/src/ca_chinook/scripts/*.pl; do ln -s $file; done;
for file in ~/local/src/ca_chinook/scripts/*.R; do ln -s $file; done;
for file in ~/local/src/ca_chinook/scripts/*.sh; do ln -s $file; done;
```
You may receive a failure message if a symbolic link already exists for a script. These error messages can be ingored.
```
## Example error message:
# ln: failed to create symbolic link './caChinookPipeline.sh': File exists
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

3. From here, you should be able to run the entire pipeline using the `caChinookPipeline.sh` script. This requires three inputs given after the script name, as shown below.
```
Example Usage: caChinookPipeline.sh <projectNumber> <runNumber> <cores>
```
In this example, we will use p134 as the projectNumber, run001 as the runNumber, and 8 for the number of cores. This last number (cores) will be provided to the Snakemake pipeline. You can adjust the number of cores up or down depending upon how many physical processor cores your computer has. If running a Windows machine, I recommend checking the 'System Information' window for your computer for the number of processor cores. I recommend that you don't exceed this number. I have gone as high as 16 cores (the number of physical processor cores in my laptop). You can run the entire pipeline with the command shown below:

```
caChinookPipeline.sh p134 run001 8
```
If the pipeline runs successfully you can skip ahead to [Processed Outputs](#output). If this script doesn't work, or if you just want to run through the pipeline manually, continue with step 4. 

4. Activate your conda environment
```
conda activate snakemake
```

5. Run the `preprocess.R` script to create the two files required by the snakemake pipeline: `samples.csv` and `units.csv`. Watch carefully for errors, and correct any problems until you successfully create the two files.
```
preprocess.R -f SampleSheet.csv
```

6. This pipeline reads from / writes to many files simultaneously. Check the number of files that can be open simultaneously on your computer with the `ulimit` command.
```
ulimit -n
```
If this reports a low number (e.g., 256), then set it to something sufficiently high (e.g., 4096) with the following command:
```
ulimit -n 4096
```

7. From the `mega-simple-microhap-snakeflow` directory, run the following command to execute the pipeline. You can adjust the number of cores (currently set at 8) up or down depending upon how many physical processor cores your computer has. If running a Windows machine, I recommend checking the 'System Information' window for your computer for the number of processor cores. I recommend that you don't exceed this number. I have gone as high as 16 cores (the number of physical processor cores in my laptop). The snakemake command will take a while, perhaps an hour or more, especially when running the pipeline for the first time.
```
cd ~/local/src/mega-simple-microhap-snakeflow/
snakemake --config run_dir=data/run001 --configfile config/Chinook/config.yaml --use-conda --cores 8
```

<hr>

## Processing the Snakemake Pipeline Output <a name="processing"></a>
This section explains how to obtain and modify files so that they can be read by some R scripts to produce the final genotypes file. I wrote a script that should handle retrieving the snakemake pipeline output from WSL, copying it to your OneDrive, and automatically handle some of the file filtering and file conversions. For example, if you are working on the first run (run001) of p134, then you can execute this script with the following command:
```
caChinookCopyAndProcess.sh p134 run001
```
If the script works, you can skip to [step 9](#step9). Otherwise, start at step 1. Most procedures in this section can alternatively be accomplished by dragging and dropping files in the Windows interface, but VCF filtering and another file modification require the command line. 

1. It now may be easier to copy some files to a new location to make them easier to access. For example, you could make a folder in your OneDrive to hold these files. I suggest opening your OneDrive folder in WSL and running the following command. Substitute the name of your run for the example `run001` in the first line of the command block below.
```
RUN="run001"
mkdir -p p134/$RUN/snakemake_output p134/$RUN/processing/
```
This will create the following directory structure. The `snakemake_output` folder is where you will copy the unmodified outputs of the snakemake pipeline. The `processing` folder is where you will run a series of scripts on select output files to filter and convert the microhaplotype data to a usable form.
```
p134
└── run001
    ├── processing
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

5. Copy the `CH_run001_greb1_q20dp5.recode.vcf` file to `p134/$RUN/processing`
```
cp CH_run001_greb1_q20dp5.recode.vcf /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/processing/.
```

6. Get the .csv file that contains read counts for the sex-linked marker. This will be the `ordered-read-counts-table.csv` file in the `p134/$RUN/Chinook/idxstats/target_fastas/ROSA/rosawr` folder. This .csv file should contain a column named `sdy_I183` which has read counts for the sex-linked marker. The number of reads sequenced per individual will be converted into female/male calls. Copy this file to the `p134/$RUN/processing` folder.
```
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/snakemake_output/idxstats/target_fastas/ROSA/rosawr
cp ordered-read-counts-table.csv /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/processing/.
```

7. Copy the .rds file containing microhaplotype data for the LFAR, WRAP, vgll3, and six6 loci into the `p134/$RUN/processing` folder.
```
cd ~/local/src/mega-simple-microhap-snakeflow/data/run001trim/Chinook/microhaplot
cp FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/p134/$RUN/processing/.
```

8. Now do the same for the .rds file containing the rest of the microhaplotype loci.
```
cp FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/.
``` 

9. <a name="step9"></a>Checkpoint: Your `processing` folder should be arranged as you see below. Verify this is accurate before proceeding.
```
processing
├── CH_run001_greb1_q20dp5.recode.vcf
├── FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds
├── FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds
└── ordered-read-counts-table.csv
```

10. Open a new Linux shell from inside your `p134/$RUN/processing` directory. An easy way to do this is to navigate to this directory in Windows Explorer, hold the `shift` key while right-clicking somewhere in the empty space in this window, then choose "Open Linux shell here" from the right-click menu. Activate your snakemake conda environment and run the `microhaplotopia.R` script from this directory. This will convert the files you generated and/or copied in the past several steps into a usable format. You will also need to provide the `microhaplotopia.R` script with your run name (e.g., `run001`). 
```
conda activate snakemake
microhaplotopia.R -r CH_run001_greb1_q20dp5.recode.vcf -R run001
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
                minimum locus sequencing depth (default = 10)

        -f FINALOUT, --finalOut=FINALOUT
                name for final output .csv file (default = haps_2col_final.csv)

        -g GREB1ROSAOUT, --greb1rosaOut=GREB1ROSAOUT
                output file for greb1rosa haplotype string (default = greb1rosa_all_hapstr.txt)

        -G GREBINFO, --grebInfo=GREBINFO
                info needed to generate greb1rosa haplotype string (default = ~/local/src/ca_chinook/example_files/greb1_roha_alleles_reordered_wr.txt)

        -L LFAR, --lfar=LFAR
                rds file of loci mapped to full genome (default = FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds)

        -m MICROHAPLOT, --microhaplot=MICROHAPLOT
                microhaplot directory name (default = microhaplot)

        -o OUTPUT, --output=OUTPUT
                output directory name (default = output)

        -O REPORTS, --reports=REPORTS
                reports directory name (default = reports)

        -r ROSA, --rosa=ROSA
                ROSA VCF file (No default; required)

        -R RUN, --run=RUN
                Run number (No default; required)

        -s SDY, --sdy=SDY
                sdy read count file (default = ordered-read-counts-table.csv)

        -S SNPLICON, --snplicon=SNPLICON
                rds file of loci mapped to target fastas (default = FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds)

        -x MAXFEMALE, --maxFemale=MAXFEMALE
                maximum sdy reads allowed to call female (default = 2)

        -X MAXUNKNOWN, --maxUnknown=MAXUNKNOWN
                maximum number of sdy reads for unknown (default = 5)

        -h, --help
                Show this help message and exit
```

Outputs from the `microhaplotopia.R` script are discussed in the next section.

11. [OPTIONAL] Convert the locus names in your final genotype .csv file. Change directories into the `processing/output` directory and run the `caChinookRenameLoci.pl` script. This will output a new file with the loci named according to the AmpliconName field in [this table](https://github.com/eriqande/california-chinook-microhaps/blob/main/inputs/Calif-Chinook-Amplicon-Panel-Information.csv). The new file will have `lociRenamed` inserted into its file name before the .csv extension.
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
├── FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds
├── FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds
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

Load your conda environment:
```
conda activate bcl2fastq
```

Run bcl2fastq:
```
bcl2fastq --barcode-mismatches 0
```

<hr>

## Genetic Stock ID
### First time setup <a name="installrubias"></a>
Activate your `snakemake` conda environment (if not already active).
```
conda activate snakemake
```

Run the following to install rubias:
```
R --slave -e "install.packages('rubias', dependencies=TRUE, repos='http://cran.rstudio.com')"
```

<hr>

### Running rubias <a name="rubias"></a>
Place your baseline and mixture files together in the same directory. Run rubias on these files with the following example command. Replace .csv file names with your actual file names:
```
rubias.R -m mixtureFile.csv -b baselineFile.csv
```

<hr>

### Rubias Outputs <a name="rubiasout"></a>
The `output_final` directory will be created in the folder from which you executed the `rubias.R` script. It will contain the following outputs:
1. `all_top_repgroup_sumPofZ.csv`
2. `all_top3pops.csv`
3. `all_toppop.csv`
4. `final_duplicates.csv`

<hr>

## CKMRsim
### First time setup <a name="installCKMRsim"></a>
Activate your `snakemake` conda environment (if not already active).
```
conda activate snakemake
```

Download and install CKMRsim package
```
wget -O CKMRsim.zip https://github.com/eriqande/CKMRsim/archive/master.zip
R --slave -e "devtools::install_local('CKMRsim.zip', upgrade='never')"
R --slave -e "install.packages('gRbase', dependencies=TRUE, repos='http://cran.rstudio.com')"
```

## Colony
### First time setup <a name="installColony"></a>
Go to your home directory
```
cd ~/local/src
```

Download the most recent Colony package for Linux. If the command below doesn't work, select most recent colony package for linux from https://www.zsl.org/about-zsl/resources/software/colony/thank-you
```
wget https://cms.zsl.org/sites/default/files/2025-07/colony2_Lnx_15_07_2025.zip
```

Unzip the file
```
unzip colony2_Lnx_15_07_2025.zip
```

Link the `colony2p.ifort.impi2018.out` file in your `~/local/bin folder` as `colony`
```
ln -s ~/local/src/colony2_Lnx_15_07_2025/colony2p.ifort.impi2018.out colony
```

Setup a conda environment for colony that includes the Intel MPI runtime library
```
conda create -n colony
conda activate colony
conda install -c https://software.repos.intel.com/python/conda/ -c conda-forge impi_rt
```

