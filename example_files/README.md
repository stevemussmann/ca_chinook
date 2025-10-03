# Example Files

This directory contains example files used for processing of data in this pipeline. 

## SampleSheet.csv
Description:
* This file provides an example for demultiplexing data through the Illumina MiSeq or through Illumina's `bcl2fastq` program.
* It is also provided to a preprocessing script to generate the `samples.csv` and `units.csv` files needed for the snakemake pipeline.

Recommendations and Troubleshooting for a MiSeq run:
* Because this is a comma-delimited file, the MiSeq will refuse to properly process the comma-delimited data in the 'comment' field. Delete the data from this field (i.e., `"ROSA, FullPanel"`) before uploading a SampleSheet to the MiSeq. Add these data back into the comment field before running the snakemake pipeline.
* Make sure the 'Sample_Well' field uses 3-character alphanumeric designations for all wells (e.g., A01 instead of A1). Otherwise the MiSeq will refuse to read the SampleSheet properly.
* If desired, you can output all Chinook microhaplotype panel options by placing the following text (including quotes) in the 'comment' field: `"LFAR, WRAP, VGLL3SIX6, ROSA, TRANSITION, SNPS84, FullPanel"`. However, only `"ROSA, FullPanel"` is necessary to obtain the full microhaplotype panel output. 

Recommendations and Troubleshooting for a NextSeq run:
* This SampleSheet format will not be compatible with the NextSeq 1000/2000 models. Demultiplex these runs using Illumina's `bcl2fastq` program.
* If processing a NextSeq run, only use the `"ROSA, FullPanel"` options for the snakemake pipeline. Adding additional Chinook microhaplotype panel options will dramatically increase the amount of hard drive space required for processing by the snakemake pipeline. 

## greb1_roha_alleles_reordered_wr.txt
Description:
* This file was provided by Anthony Clemento and is used by the pipeline to ensure all loci in the ROSA haplotype string are correctly coded and output in the correct order.
