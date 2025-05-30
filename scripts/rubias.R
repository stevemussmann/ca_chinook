#!/usr/bin/env Rscript

## Most of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

library("tidyverse", quietly=TRUE)
library("rubias")
library("optparse", quietly=TRUE)

## command line option parsing
option_list = list(
  make_option(
    c("-b", "--baseline"), 
    type="character", 
    default="SWFSC-chinook-reference-baseline.select.csv", 
    help="name for baseline genotypes .csv file (default = SWFSC-chinook-reference-baseline.select.csv)", 
    metavar="mixture"
  ),
  make_option(
    c("-m", "--mixture"), 
    type="character", 
    default="haps_2col_final.lociRenamed.csv", 
    help="name for mixture genotypes .csv file (default = haps_2col_final.lociRenamed.csv)", 
    metavar="mixture"
  ),
  make_option(
    c("-o", "--output"), 
    type="character", 
    default="output_final", 
    help="output directory name (default = output_final)", 
    metavar="output"
  )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$mixture)){
#  print_help(opt_parser)
#  stop("mixture genotypes file must be specified (option -m).", call.=FALSE)
#}

# set working directory in Rstudio
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
WD <- getwd()

################################################################################
## SET FILE AND PATH NAMES HERE
################################################################################

out <- opt$output
outDir <- file.path(WD, out)
dir.create(outDir, showWarnings=FALSE)

################################################################################

baseline <- read_csv(opt$baseline)
#baseline <- read_csv("subRoSA_baseline_final_rubias183_collapse_MDF_MDS.csv")

# And we need our genotypes file from microhaplotopia
#genos <- read_csv("Results_Interlab_Plates_1_2_3_GeM_DWR_20250521.csv")
genos <- read_csv(opt$mixture)

# We will  need some required fields for the rubias file
req <- tibble(sample_type = "mixture", repunit = NA, collection = "final")
mixfile <- genos %>% cbind(req, .)
mixfile_char <- mixfile %>% mutate_at(vars(-(sample_type:indiv)), as.character)
#View(mixfile_char)

# Found a bug in close_matching_samples that doesn't like all fish being mixture
mixfile_char[1,1] = "reference"
dups <- close_matching_samples(mixfile_char, 5, min_frac_non_miss = 0.8, min_frac_matching = 0.9)
mixfile_char[1,1] = "mixture"

dups

# Huh weird, looks like those consecutive samples match
# Going to leave it in there but will write it out so we have a record
dupsOut <- file.path(outDir, "final_duplicates.csv")
write_csv(dups, dupsOut)

# We need to make the mixfile match the loci in the baseline
common <- intersect(colnames(baseline) , colnames(mixfile_char))
baseline_ready <- baseline %>% select(all_of(common))
mixfile_char_ready <- mixfile_char %>% select(all_of(common))

# calculate the genetic mixture
mix_est <- infer_mixture(reference = baseline_ready, 
                         mixture = mixfile_char_ready,
                         #method = "BR",
                         #reps = 50000,
                         #burn_in = 10000,
                         gen_start_col = 5)

# Most likely single populatioon
topassign <- mix_est$indiv_posteriors %>% 
  group_by(indiv) %>% 
  arrange(desc(PofZ)) %>% 
  slice(1)


topOut <- file.path(outDir, "all_toppop.csv")
write_csv(topassign, topOut)

# Top 3 likely populations
top3assign <- mix_est$indiv_posteriors %>% 
  group_by(indiv) %>% 
  arrange(desc(PofZ)) %>% 
  slice(1:3)

top3Out <- file.path(outDir, "all_top3pops.csv")
write_csv(top3assign, top3Out)

# Sum over pop assignments to reporting_group for individuals posteriors
rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ)) %>% 
  group_by(indiv) %>% 
  arrange(desc(rep_pofz)) %>% 
  slice(1)

repOut <- file.path(outDir, "all_top_repgroup_sumPofZ.csv")
write_csv(rep_indiv_ests, repOut)

quit()