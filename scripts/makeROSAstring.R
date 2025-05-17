#!/usr/bin/env Rscript

## Nearly all of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

# set working directory in Rstudio and capture path in WD
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
WD <- getwd()

library("tidyverse")
library("vcfR")

# read files
vcfFile <- "CH_test_greb1_q20dp5.recode.vcf" #your input vcf file
locs <- read.table(file="greb1_roha_alleles_reordered_wr.txt", header=T, sep="\t", stringsAsFactors = F) %>% as_tibble()

v <- read.vcfR(vcfFile)
vt <- vcfR2tidy(v)

# Note: E (early), L (late), W (winter), N (non-winter), H (heterozygous) and ? (missing)
# Extract the individual genotype calls at each SNP
gtypes <- vt$gt %>% select(ChromKey, POS, Indiv, gt_GT_alleles, gt_DP, gt_AD) %>%
  left_join(vt$fix %>% select(ChromKey, CHROM, POS)) %>%
  unite("CHROMPOS", c("CHROM", "POS"), sep="_", remove = T) %>%
  bind_rows(.id = "vcf_file") %>%
  filter(CHROMPOS %in% locs$CHROMPOS) %>%
  group_by(Indiv) %>%
  arrange(match(CHROMPOS, locs$CHROMPOS), .by_group = TRUE) %>%
  ungroup()

ptypes <- gtypes %>%
  select(CHROMPOS, Indiv, gt_GT_alleles) %>%
  complete(CHROMPOS, Indiv) %>%
  replace_na(list(gt_GT_alleles = "REF")) %>%
  left_join(locs) %>%
  group_by(Indiv) %>%
  distinct(Indiv, CHROMPOS, .keep_all = TRUE) %>%
  arrange(haporder, .by_group = TRUE) %>%
  summarise(hapstr = paste(pheno, collapse = "")) %>%
  mutate(hapstr = str_replace(hapstr, "NA", "?"))
#View(ptypes)

ptype_sum <- ptypes %>% group_by(hapstr) %>%
  summarise(total = n())
#View(ptype_sum)

write.table(ptypes, file="greb1rosa_all_16MAY2025.txt", quote=F, sep="\t", row.names=F)

#quit()