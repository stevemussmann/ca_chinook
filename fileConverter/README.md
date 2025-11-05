# fileConverter

Program for filtering microhaplotype data and converting to Colony (and other) file formats. 

This program is a work in progress. This README.md file will be updated regarding program functionality as features are added. Currently the program serves two main functions:
1. Filtering microhaplotype data output by the snakemake pipeline to remove loci and individuals with excessive missing data.
2. Converting the data to Colony input format. It will also output a filtered .csv file if requested.

## Python Version Compatibility
This program has only been tested in Python v3.12.10. However, it should be compatible with Python v3.8+.

## Dependencies
- pandas

## Installation
Activate the `snakemake` conda environment that you made for processing microhaplotype data and install the `pandas` dependency.

```
conda activate snakemake
conda install pandas
```

`microhapConvert.py` was downloaded when you installed this repository on your computer. To make it easily accessible from the command line, make the script executable and install it in your path.

```
cd ~/local/src/ca_chinook/fileConverter
chmod u+x microhapConvert.py
cd ~/local/bin/
ln -s $HOME/local/src/ca_chinook/fileConverter/microhapConvert.py
```

You can test if the converter is properly installed by attempting to print the help menu:
```
microhapConvert.py --help
```

## Input Requirements
Briefly, the haps_2col_final.csv file produced at the end of the `caChinookPipeline.sh` script is acceptable as input, with some small modifications (see below). 

In other words, the input should be in csv format with two columns for each locus. Headings for the two alleles should end in _1 and _2 (e.g., Ots_100884-287_1 and Ots_100884-287_2). The first column should be titled `indiv` and contain the individual IDs for all samples. You do not need to remove the special columns from the file that were inserted by the genotyping pipeline (`sdy_sex`, `hapstr`, `rosa_pheno`, and `percMicroHap`) because these will be automatically identified and removed by the file converter. 

Modifications: You do, however, need to add some information for certain programs. 
* The second column should be titled exactly `Population ID` and this column should contain population information for all individuals.
* The third column in the file should be titled exactly `colony2`. In this column, identify all potential offspring as `offspring` and all candidate parents by their sex (`male` or `female`).
* Additional columns are needed for SNPPIT and Sequoia conversions (see SNPPIT and Sequoia sections at the end).

## Program Options
Required Inputs:
* **`-f` / `--infile`:** Specify the input file in .csv format.
* **`-r` / `--runname`:** Provide a unique name for the Colony run. Output files from Colony will receive this name.

Required for SNPPIT conversion only:
* **`-Z` / `--snppitmap`:** Specify a tab-delimited map in which the first column lists each population, the second column lists its status as POP or OFFSPRING, and the third column lists the potential parental POP(s) for each OFFSPRING.

Filtering options:
* **`-i` / `--pmissind`:** Enter the maximum allowable proportion of missing data for an individual (default = 0.3).
* **`-l` / `--pmissloc`:** Enter the maximum allowable proportion of missing data for a locus (default = 0.3).
* **`-R` / `--removeloci`:** Specify a list of loci to remove (input = plain text file, one locus per line).
* **`-m` / `--mono`:** Remove monomorphic loci from final output (default = True).

Arguments that apply to colony-format outputs only:
* **`-d` / `--droperr`:** Enter the assumed allelic dropout rate (default = 0.0005).
* **`-g` / `--genoerr`:** Enter the assumed genotyping error rate (default = 0.0005).
* **`-I` / `--inbreed`:** Accepts only 0 and 1 as input. 0 = inbreeding absent; 1 = inbreeding present (default = 0).
* **`-L` / `--runlength`:** Accepts only 1, 2, 3, or 4 as input. 1/2/3/4 = Short/Medium/Long/VeryLong run (default = 2).
* **`-M` / `--pmale`:** Enter the assumed probability of father being among candidate parents (default = 0.5). Value is ignored if no candidate fathers provided in the dataset.
* **`-F` / `--pfemale`:** Enter the assumed probability of mother being among candidate parents (default = 0.5). Value is ignored if no candidate mothers provided in the dataset.

File format conversion output options (at least one required):
* **`-c` / `--csv`:** Write filtered csv format file.
* **`-C` / `--colony`:** Write colony format file.
* **`-q` / `--sequoia`:** Prints a [sequoia](https://jiscah.github.io/) formatted genotype file.
* **`-z` / `--snppit`:** Prints a file in [snppit](https://github.com/eriqande/snppit) format (-Z option is also required for snppit conversion as specified above).

## Running the program
First activate your snakemake pipeline.
```
conda activate snakemake
```

As an example, the following command could be run to produce both a filtered .csv file and a file in Colony format:
```
microhapConvert.py -f haps_2col_final.csv -r example -c -C
```

## Outputs
Converted and filtered files are written to a `convertedFiles/` output directory that is created in the directory from which you executed the conversion program. The table below provides the output file names from each option. If only an extension is provided, the the input file basename is retained and the new file extension is applied to the output. This results in the .csv format option producing a file with the same name as the input file, but because it is written to the `convertedFiles/` output directory the original file will be retained. A file named `locusDictionary.json` will also be written. This contains the necessary information to translate the Colony haplotype numbers back into nucleotide format. 

<div align="center">
  
| Format       | Output file name and/or extension  |  Program Option  |
| :----------- | :--------------------------------: | :--------------: |
| Colony       | Colony2.Dat                        | `-C`             |
| CSV          | .csv                               | `-c`             |
| Sequoia      | .sequoia; sequoia.lh.txt           | `-q`             |
| SNPPIT       | .snppit                            | `-z`             |
  
</div>

### Sequoia
The Sequoia conversion relies upon some of the optional SNPPIT columns that are also used for the SNPPIT file conversion (see below). Use the POPCOLUMN_SEX column to specify sex data for all individuals. Only case insensitive versions of `f`, `female`, `m`, and `male` will be recognized. All other values and blank cells will be converted to unknown sex data value in sequoia (3). 

The OFFSPRINGCOLUMN_BORN_YEAR is used to specify the birth year for all individuals. You can enter birth year data in this column even for the 'parental' populations. This will not cause any problems for the SNPPIT file conversion as listed below.

The code for creating the life history data file has not yet been robustly tested, so there could be bugs.

Files can be read into sequoia with the following commands:
```
library("sequoia")

# genotypes file
geno <- as.matrix(read.csv("filename.sequoia", sep="\t", header=FALSE, row.names=1))

# life history file
lh <- read.csv("sequoia.LH.txt", sep="\t", header=TRUE)
```

### SNPPIT
The SNPPIT conversion has a few special requirements that are not needed for other file formats. Firstly, a special tab-delimited snppit map file is required as supplemental input. Essentially, each line of this file is intended to contain all of the information of lines starting with the POP and OFFSPRING keywords, as seen on [pages 22-23 of the SNPPIT program documentation](https://github.com/eriqande/snppit/blob/master/doc/snppit_doc.pdf). However, note that the columns of the snppitmap are in a different order than they appear in the final snppit-formatted file (i.e., 'popname'\tab'POP' rather than 'POP'\tab'popname').

Secondly, the user can utilize SNPPIT's 'optional' columns as seen on [pages 24-26 of the SNPPIT program documentation](https://github.com/eriqande/snppit/blob/master/doc/snppit_doc.pdf) by including the relevant data in their input .xlsx file. To do this, add columns to your input .xlsx file with headings that exactly match the optional columns used by SNPPIT. For example, if you want to use the POPCOLUMN_SEX option in SNPPIT, then include a column named exactly POPCOLUMN_SEX in the 'Final Genotypes' worksheet of your input .xlsx file, and fill this column with the appropriate values for each individual, as necessary.

Generally, the values you input in these optional columns should exactly match the values as they would appear in the final SNPPIT file. Additional details are as follows:
* Values in the POPCOLUMN_SEX column are somewhat flexible. Case-insensitive versions of 'f' and 'female' or 'm' and 'male' will be converted to 'F' and 'M' respectively. Blank cells and any other values entered in this column will be converted to missing data ('?').
* Values in columns containing year data (POPCOLUMN_REPRO_YEARS, OFFSPRINGCOLUMN_BORN_YEAR, OFFSPRINGCOLUMN_SAMPLE_YEAR) must be valid four-digit integers.
* Currently there are no data validation measures implemented for the other two optional columns (POPCOLUMN_SPAWN_GROUP and OFFSPRINGCOLUMN_AGE_AT_SAMPLING) so please make sure anything you enter in these columns is exactly as you want it to appear in the final file.
* Unnecessary data can be left as blank cells. An example of this is that you would not need to enter sex data for OFFSPRING groups in the POPCOLUMN_SEX column.
