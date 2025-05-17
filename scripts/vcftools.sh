#!/bin/bash

VCF=$1
QUAL=20
DP=5

vcftools --vcf $VCF --minQ $QUAL --minDP $DP --remove-indels --recode --recode-INFO-all --out CH_test_greb1_q20dp5

exit
