# colonyConverter

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

`colonyConverter.py` was downloaded when you installed this repository on your computer. To make it easily accessible from the command line, make the script executable and install it in your path.

```
cd ~/local/src/ca_chinook/colonyConverter
chmod u+x colonyPrep.py
cd ~/local/bin/
ln -s $HOME/local/src/ca_chinook/colonyConverter/colonyPrep.py
```

You can test if the converter is properly installed by attempting to print the help menu:
```
colonyPrep.py --help
```

## Input Requirements
Briefly, the haps_2col_final.csv file produced at the end of the `caChinookPipeline.sh` script is acceptable as input, with some small modifications (see below). 

In other words, the input should be in csv format with two columns for each locus. Headings for the two alleles should end in _1 and _2 (e.g., Ots_100884-287_1 and Ots_100884-287_2). The first column should be titled `indiv` and contain the individual IDs for all samples. You do not need to remove the special columns from the file that were inserted by the genotyping pipeline (`sdy_sex`, `hapstr`, `rosa_pheno`, and `percMicroHap`) because these will be automatically identified and removed by the file converter. 

Modifications: You do, however, need to add some information for Colony2. The second column in the file should be titled exactly `colony2`. In this column, identify all potential offspring as `offspring` and all candidate parents by their sex (`male` or `female`). 

## Program Options
Required Inputs:
* **`-f` / `--infile`:** Specify the input file in .csv format.
* **`-r` / `--runname`:** Provide a unique name for the Colony run. Output files from Colony will receive this name.

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

## Running the program
First activate your snakemake pipeline.
```
conda activate snakemake
```

As an example, the following command could be run to produce both a filtered .csv file and a file in Colony format:
```
colonyPrep.py -f haps_2col_final.csv -r example -c -C
```

## Outputs
Converted and filtered files are written to a `convertedFiles/` output directory that is created in the directory from which you executed the conversion program. The table below provides the output file names from each option. If only an extension is provided, the the input file basename is retained and the new file extension is applied to the output. This results in the .csv format option producing a file with the same name as the input file, but because it is written to the `convertedFiles/` output directory the original file will be retained. A file named `locusDictionary.json` will also be written. This contains the necessary information to translate the Colony haplotype numbers back into nucleotide format. 

<div align="center">
  
| Format       | Output file name and/or extension  |  Program Option  |
| :----------- | :--------------------------------: | :--------------: |
| Colony       | Colony2.Dat                        | `-C`             |
| CSV          | .csv                               | `-c`             |
  
</div>
