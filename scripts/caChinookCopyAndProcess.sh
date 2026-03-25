#!/bin/bash

# exit and print usage if command line input not provided.
if [ $# -lt 2 ]
then
	echo -e "\nUsage: caChinookCopyAndProcess.sh <projectNumber> <runNumber>\n"
	exit
fi

PRJ=$1
RUN=$2
WINDOWSUSER=`powershell.exe '$env:UserName' | sed 's/\r//g'` #get user's windows username

# make directories that data will be copied into
mkdir -p /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/snakemake_output /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing

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

# copy .rds files into proper location
cd $HOME/local/src/mega-simple-microhap-snakeflow/data/run001trim/Chinook/microhaplot
cp FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/. # copy loci mapped to full genome
cp FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds /mnt/c/Users/$WINDOWSUSER/OneDrive\ -\ DOI/$PRJ/$RUN/processing/. # copy loci mapped to target fastas

exit
