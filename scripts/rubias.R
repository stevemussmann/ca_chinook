#!/usr/bin/env Rscript

## Most of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

suppressPackageStartupMessages(library("tidyverse", quietly=TRUE))
library("rubias")
library("optparse", quietly=TRUE)

## command line option parsing
option_list = list(
  make_option(
    c("-b", "--baseline"), 
    type="character", 
    default="SWFSC-chinook-reference-baseline.select.csv", 
    help="name for baseline genotypes .csv file (default = SWFSC-chinook-reference-baseline.select.csv)", 
    metavar="baseline"
  ),
  make_option(
    c("-m", "--mixture"), 
    type="character", 
    default="haps_2col_final.csv", 
    help="name for mixture genotypes .csv file (default = haps_2col_final.csv)", 
    metavar="mixture"
  ),
  make_option(
    c("-o", "--output"), 
    type="character", 
    default="output_final", 
    help="output directory name (default = output_final)", 
    metavar="output"
  ),
  make_option(
  	c("-p", "--missing"),
	type="numeric",
	default=50.0,
	help="minimum genotype percent threshold to retain individual (default = 50.0%)",
	metavar="percent"
  )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#if (is.null(opt$mixture)){
#  print_help(opt_parser)
#  stop("mixture genotypes file must be specified (option -m).", call.=FALSE)
#}

WD <- getwd()

################################################################################
## SET FILE AND PATH NAMES HERE
################################################################################

out <- opt$output
outDir <- file.path(WD, out)
dir.create(outDir, showWarnings=FALSE)

################################################################################

# read baseline data
baseline <- read_csv(opt$baseline, col_types = cols(.default = "c"))

# And we need our genotypes file from microhaplotopia
genos <- read_csv(opt$mixture, col_types = cols(.default = "c")) # read in all columns as character data as default 

# We will  need some required fields for the rubias file
req <- tibble(sample_type = "mixture", repunit = NA, collection = "final")
mixfile <- genos %>% cbind(req, .)
mixfile_char <- mixfile %>% mutate_at(vars(-(sample_type:indiv)), as.character)

# filter to remove individuals missing high proportion of loci
if( "percMicroHap" %in% names(mixfile_char) ){
	mixfile_char <- mixfile_char %>% mutate(percMicroHap = as.numeric(percMicroHap)) # make sure percMicroHap value is numeric
	mixfile_char <- mixfile_char %>% filter(percMicroHap >= opt$missing)
}

# extract ROSA info to use later for final assignments
rosa_columns <- c("indiv", "hapstr", "hapstr_dist", "canonical_rosa_pheno")
rosa_data <- mixfile_char %>% select(all_of(rosa_columns)) 

# extract LFAR markers to use later for final assignments
lfar_columns <- c("indiv", "NC_037130.1:864908-865208", "NC_037130.1:1062935-1063235")
lfar_data <- mixfile_char %>% select(contains(lfar_columns)) %>% mutate(lfar_na_count = rowSums(is.na(.)))

# drop columns that may or may not exist
columns_to_drop <- c("sdy_model_sex", "hapstr", "hapstr_dist", "canonical_rosa_pheno", "percMicroHap", "perc_Xtra")
mixfile_char <- mixfile_char %>% select(-any_of(columns_to_drop))

# Found a bug in close_matching_samples that doesn't like all fish being mixture
mixfile_char[1,1] = "reference"
dups <- close_matching_samples(mixfile_char, 5, min_frac_non_miss = 0.8, min_frac_matching = 0.9)
mixfile_char[1,1] = "mixture"

# write out records of duplicates
dupsOut <- file.path(outDir, "final_duplicates.csv")
write_csv(dups, dupsOut)

# We need to make the mixfile match the loci in the baseline
common <- intersect(colnames(baseline) , colnames(mixfile_char))
baseline_ready <- baseline %>% select(all_of(common))
mixfile_char_ready <- mixfile_char %>% select(all_of(common))

## test if there are any loci remaining that are common to both inputs and exit if <=4
# This is intended to catch cases in which no locus names are shared by the two inputs
if( length( common ) <= 4){
	cat("\nThese are the only columns shared by the mixture and baseline files:\n")
	print(common)
	stop("Check that your mixture and baseline files use the same locus names.\n\n")
}

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

# combine topassign and rep_indiv_ests with rosa_data
temp <- left_join(rep_indiv_ests, topassign, by = "indiv") %>% 
  select(-c(mixture_collection.x, mixture_collection.y, repunit.y, log_likelihood, z_score, missing_loci)) %>%
  rename(repunit = repunit.x) %>% rename(collection_PofZ = PofZ)
temp2 <- left_join(temp, rosa_data, by="indiv")
classification_data <- left_join(temp2, lfar_data, by="indiv")

finalClass <- classification_data %>% mutate(final_reporting_group = case_when(
  ## GSI = Winter
  repunit == "winter" ~ "Winter",
  
  ## GSI = Late fall
  repunit == "sacfall" & collection == "ColemanLF" & (canonical_rosa_pheno == "Spring" | canonical_rosa_pheno == "Winter" | canonical_rosa_pheno == "Sp-Win") ~ "FRLspring", #ROSA = early/early
  repunit == "sacfall" & collection == "ColemanLF" & (canonical_rosa_pheno == "Sp-Fall" | canonical_rosa_pheno == "Sp-Win") ~ "Lfall", #ROSA = early/late
  repunit == "sacfall" & collection == "ColemanLF" & canonical_rosa_pheno == "Fall" ~ "Lfall", #ROSA = late/late
  # sacfall & ColemanLF with missing ROSA will be handled by .default = "Unknown"
  
  ## GSI = Fall
  repunit == "sacfall" & collection != "ColemanLF" & (canonical_rosa_pheno == "Spring" | canonical_rosa_pheno == "Winter" | canonical_rosa_pheno == "Sp-Win") & lfar_na_count < 4 ~ "FRLspring", # early/early
  repunit == "sacfall" & collection != "ColemanLF" & (canonical_rosa_pheno == "Sp-Fall" | canonical_rosa_pheno == "Sp-Win") & lfar_na_count < 4 ~ "Fall", #ROSA = early/late
  repunit == "sacfall" & collection != "ColemanLF" & canonical_rosa_pheno == "Fall" & lfar_na_count < 4 ~ "Fall", #ROSA = late/late
  # sacfall missing ROSA will be handled by .default = "Unknown"
  repunit == "sacfall" & collection != "ColemanLF" & (canonical_rosa_pheno == "Spring" | canonical_rosa_pheno == "Winter" | canonical_rosa_pheno == "Sp-Win") & lfar_na_count == 4 ~ "FRLspring", #early/early
  repunit == "sacfall" & collection != "ColemanLF" & canonical_rosa_pheno == "Fall" & lfar_na_count == 4 ~ "Fall or Late fall", #ROSA = late/late
  # sacfall missing ROSA and missing LFAR will be handled by .default = "Unknown"
  
  ## GSI = MDspring
  repunit == "sacspring" & collection == "MillDeerSp" & (canonical_rosa_pheno == "Spring" | canonical_rosa_pheno == "Winter" | canonical_rosa_pheno == "Sp-Win") ~ "MDspring", # ROSA = early/early
  repunit == "sacspring" & collection == "MillDeerSp" & (canonical_rosa_pheno == "Sp-Fall" | canonical_rosa_pheno == "Sp-Win") ~ "MDspring", #ROSA = early/late
  repunit == "sacspring" & collection == "MillDeerSp" & canonical_rosa_pheno == "Fall" ~ "Fall", #ROSA = late/late
  # sacspring & MDspring with missing ROSA will be handled by .default = "Unknown"
  
  ## GSI = Bspring
  repunit == "sacspring" & collection == "ButteSp" & (canonical_rosa_pheno == "Spring" | canonical_rosa_pheno == "Winter" | canonical_rosa_pheno == "Sp-Win") ~ "Bspring", # early/early
  repunit == "sacspring" & collection == "ButteSp" & (canonical_rosa_pheno == "Sp-Fall" | canonical_rosa_pheno == "Sp-Win") ~ "Bspring", #ROSA = early/late
  repunit == "sacspring" & collection == "ButteSp" & canonical_rosa_pheno == "Fall" ~ "Fall", #ROSA = late/late
  # sacspring & ButteSp with missing ROSA will be handled by .default = "Unknown"
  
  # all other cases
  .default = "Unassigned"
))

# write final classifications
finalOut <- file.path(outDir, "final_classifications.csv")
write_csv(finalClass, finalOut)

cat("\n\nFinal output written to", finalOut, "\n\n")

# add back in number of present and missing loci
rep_indiv_ests <- left_join(rep_indiv_ests, topassign %>% select(indiv, n_non_miss_loci, n_miss_loci, missing_loci), by = "indiv")

repOut <- file.path(outDir, "all_top_repgroup_sumPofZ.csv")
write_csv(rep_indiv_ests, repOut)

quit()
