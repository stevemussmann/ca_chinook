#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: vcftools.sh <input.vcf> <runNumber>"
	echo ""
	echo "Input help: your input file will probably be named \"variants-bcftools.vcf\""
	echo "Output help: your output will be named using the run number."
	echo ""
	echo "For example, if you input run number \"run001\", then your output file will be named CH_run001_greb1_q20dp5.vcf"
	echo ""
	exit
fi

VCF=$1
QUAL=20
DP=5
OUT="CH_${2}_greb1_q${QUAL}dp${DP}"

vcftools --vcf $VCF --minQ $QUAL --minDP $DP --remove-indels --recode --recode-INFO-all --out $OUT

exit
