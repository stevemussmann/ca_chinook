#!/bin/bash

# exit and print usage if command line input not provided.
if [ $# -lt 2 ]
then
	echo "Usage: caChinookCopyAndProcess.sh <projectNumber> <runNumber>\n"
	exit
fi

PRJ=$1
RUN=$2
WINDOWSUSER=`powershell.exe '$env:UserName' | sed 's/\r//g'` #get user's windows username

# make directories that data will be copied into
mkdir -p /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/microhaplot

# copy data into onedrive snakemake_output directory
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/idxstats/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/.
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/microhaplot/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/.
cp -r ~/local/src/mega-simple-microhap-snakeflow/data/$RUN/Chinook/vcfs/ /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/.

# change directories to location of ROSA vcf file
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/vcfs/ROSA/target_fasta/rosawr

# filter vcf file
caChinookVcftools.sh variants-bcftools.vcf $RUN

# copy filtered vcf file to new location in onedrive
cp CH_${RUN}_greb1_q20dp5.recode.vcf /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/.

# copy file containing sdy read counts to new location in onedrive
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/idxstats/target_fastas/ROSA/rosawr
cp ordered-read-counts-table.csv /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/.

# fix R shiny files so they will run properly with recent version of ggiraph
cd /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output/microhaplot
fixMicrohaplot.sh

exit
