#!/usr/bin/env Rscript

## Most of the code in this script was borrowed extensively from code provided 
## by Anthony Clemento at NOAA and Bryan Nguyen at California DWR

## Packages that are necessary. 
library("CKMRsim")
suppressPackageStartupMessages(library("tidyverse", quietly=TRUE))
library("optparse", quietly=TRUE)

## command line option parsing
option_list = list(
  make_option(
    c("-a", "--auto"),
    action="store_true",
    default=FALSE,
    help="Turn on automatic log-likelihood threshold calculation"
  ),
  make_option(
    c("-d", "--outdir"), 
    type="character", 
    default="ckmrsim_output", 
    help="output directory name (default = ckmrsim_output)", 
    metavar="OUTPUT_DIRECTORY"
  ),
  make_option(
    c("-l", "--logl"),
    type="numeric",
    default=10,
    help="minimum log-likelihood threshold for assigning pairwise relationships (default = 10)",
    metavar="LOGL_THRESHOLD"
  ),
  make_option(
    c("-L", "--loci"), 
    type="numeric", 
    default=120, 
    help="minimum number of shared loci for assigning pairwise relationships (default = 120)", 
    metavar="SHARED_LOCI"
  ),
  make_option(
    c("-n", "--name"), 
    type="character", 
    #default="ckmrsim_output", 
    help="project name (required)", 
    metavar="PROJECT_NAME"
  ),
    make_option(
    c("-o", "--offspring"),
    type="character", 
    #default="haps_2col_final.csv", 
    help="name for unknown genotypes .csv file (required)", 
    metavar="OFFSPRING_GENOTYPES_CSV"
  ),
  make_option(
    c("-p", "--parents"),
    type="character", 
    #default="SWFSC-chinook-reference-baseline.select.csv", 
    help="name for parental genotypes .csv file (required)", 
    metavar="PARENTS_GENOTYPES_CSV"
  )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# test for use of required options (name, offspring, parents)
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Project name must be specified (option -n).", call.=FALSE)
}

if (is.null(opt$offspring)){
  print_help(opt_parser)
  stop("Offspring .csv file must be specified (option -o).", call.=FALSE)
}

if (is.null(opt$parents)){
  print_help(opt_parser)
  stop("Parents .csv file must be specified (option -p).", call.=FALSE)
}

# function to test if file exists
fileExists <- function(file) {
  if (!file.exists(file)) {
    cat("\n")
    stop(paste("\n", file, "does not exist.\n\n"))
  }
}

# test if input files exist
fileExists(opt$offspring)
fileExists(opt$parents)

# make output directory
OUTDIR <- opt$outdir
if( !dir.exists(OUTDIR)){
  dir.create(OUTDIR) 
}

# read genotypes files
parents <- read_tsv(opt$parents, col_types = cols(.default = col_character())) %>% rename(ID = indiv)
offspring <- read_tsv(opt$offspring, col_types = cols(.default = col_character())) %>% rename(ID = indiv)

# get lists of parents and offspring
parent_ids <- parents$ID
offspring_ids <- offspring$ID

# set loglikelihood threshold, min shared loci threshold, and project name
min_loci_threshold <- opt$loci
project_name <- opt$name

# combine parents and offspring genotype files
combo_geno <- bind_rows(
  parents,
  offspring
)

# fix loci with '.' in names
names(combo_geno) <- gsub("\\.", "-", names(combo_geno))

## Compute allele frequencies from genotype data
# make note of the current order of loci in the data set
nc <- ncol(combo_geno)
loci <- names(combo_geno)[seq(2, nc, by = 2)]

#  reset the locus names
names(combo_geno)[seq(3, nc, by = 2)] <- str_c(names(combo_geno)[seq(2, nc, by = 2)], "1", sep = ".")
names(combo_geno)[seq(2, nc, by = 2)] <- str_c(names(combo_geno)[seq(2, nc, by = 2)], "2", sep = ".")

# then make some long format genotypes
long_genos <- combo_geno %>% 
  gather(key = "loc", value = "Allele", -ID) %>%
  separate(loc, into = c("Locus", "gene_copy"), sep = "\\.") %>%
  mutate(Allele = as.character(Allele)) %>%
  mutate(Allele = ifelse(Allele == "0", NA, Allele)) %>%
  rename(Indiv = ID)

# and now we can compute the allele frequencies
alle_freqs <- long_genos %>%
  count(Locus, Allele) %>%
  group_by(Locus) %>%
  mutate(Freq = n / sum(n),
         Chrom = "Unk",
         Pos = as.integer(factor(Locus, levels = loci))) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele))

# and use the included function to make a dataframe of frequencies
afreqs_ready <- reindex_markers(alle_freqs)

## Create a CKMR object
PO_ckmr <- create_ckmr(
  D = afreqs_ready,
  #kappa_matrix = kappas[c("MZ", "PO", "FS", "HS", "GP", "AN", "DFC", "HAN", "FC", "HFC", "DHFC", "SC", "HSC", "U"), ],
  #kappa_matrix = kappas[c("PO", "FS", "HS", "GP", "AN", "DFC", "HAN", "FC", "U"), ],
  kappa_matrix = kappas[c("PO", "FS", "HS", "FC", "U"), ],
  ge_mod_assumed = ge_model_microhap1,
  ge_mod_true = ge_model_microhap1,
  ge_mod_assumed_pars_list = list(miscall_rate = 0.005, dropout_rate = 0.005),
  ge_mod_true_pars_list = list(miscall_rate = 0.005, dropout_rate = 0.005)
)

# Split back into parent and offspring genotype dataframes
parent_genos_long <- long_genos %>%
  filter(Indiv %in% parent_ids)
offspring_genos_long <- long_genos %>%
  filter(Indiv %in% offspring_ids)


po_results <- pairwise_kin_logl_ratios(
  D1 = parent_genos_long,
  D2 = offspring_genos_long,
  CK = PO_ckmr,
  numer = "PO",
  denom = "U"
)

## Simulate genotype pairs and calculate log-probabilities
# We are interested in our predicted power for parent-offspring, fullsibs, halfsibs, and first cousins
Qs <- simulate_Qij(PO_ckmr, 
                   #calc_relats = c("PO", "FS", "HS", "FC", "U"),
                   #sim_relats = c("PO", "FS", "HS", "FC", "U"))
                   calc_relats = c("PO", "FS", "HS", "U"),
                   sim_relats = c("PO", "FS", "HS", "U"))

if (opt$auto == FALSE) {
  logl_threshold <- as.numeric(opt$logl)
  cat("Applying user-defined log-likelihood ratio threshold of", logl_threshold, "\n")
  cat("Applying minimum shared loci threshold of", min_loci_threshold, "\n")
} else if (opt$auto == TRUE) {
  # calculate FPR threshold
  threshFPR <- .01*(length(parents)*length(offspring_ids))^(-1)
  print(threshFPR)
  
  # calculate series of lambda_star values for FPR
  lambda_0_50 <- mc_sample_simple(Qs, 
                                   nu = "PO",
                                   de = "U", 
                                   lambda_stars = seq(0, 50, by = .1))
  
  print(lambda_0_50)
  
  filtLamb_star <- lambda_0_50 %>%
    filter(FPR > threshFPR) %>%
    arrange(FPR) %>%
    slice(1)
  
  print(filtLamb_star)
  
  # test if lambda_star found; if not default to user input
  if( nrow(filtLamb_star) == 1){
    logl_threshold <- filtLamb_star$Lambda_star
  }else{
    logl_threshold <- as.numeric(opt$logl)
    cat("Could not automatically calculate log-likelihood ratio threshold.")
    cat("Defaulting to user input of", logl_threshold, "\n")
    cat("Applying minimum shared loci threshold of", min_loci_threshold, "\n")
  }

}

po_results_filtered <- po_results %>%
  filter(logl_ratio >= logl_threshold) %>%
  filter(num_loc >= min_loci_threshold) %>%
  arrange(desc(logl_ratio))

total_genos_long_dedup <- long_genos %>%
  distinct(Indiv, Locus, gene_copy, .keep_all = TRUE)

mendelian_incompatibilities <- tag_mendelian_incompatibilities(po_results_filtered, total_genos_long_dedup)

sample_MI <- mendelian_incompatibilities %>%
  filter(!is.na(is_MI)) %>%
  group_by(D2_indiv, D1_indiv) %>%
  summarize(total_incompat = sum(is_MI, na.rm = TRUE), total_compat = sum(!is_MI, na.rm = TRUE), .groups = 'keep') %>%
  mutate(fraction_incompat = total_incompat / (total_incompat + total_compat))

po_results_filtered_with_MI <- po_results_filtered %>%
  left_join(sample_MI, by = c("D2_indiv","D1_indiv"))

# make output file name, combine with outdir path, and write tsv file
outfile = paste0(project_name, "_PO_results.tsv")
outfiledir = file.path(OUTDIR, outfile)
write_tsv(po_results_filtered_with_MI, file = outfiledir)

## Plots
## Simulate genotype pairs and calculate log-probabilities
# We are interested in our predicted power for parent-offspring, fullsibs, halfsibs, and first cousins
PO_U_logls <- extract_logls(Qs,
                            numer = c(PO = 1),
                            denom = c(U = 1))

#loglsPlot <- ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "FS", "HS", "FC", "U")),
loglsPlot <- ggplot(PO_U_logls %>% filter(true_relat %in% c("PO", "FS", "HS", "U")),
                    aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.25)
outplotdir = file.path(OUTDIR, "logLplot.png")
ggsave(outplotdir, plot=loglsPlot, dpi=300)

## Estimating False Negative and False Positive Rates
# estimate false positive rates when the true relationship is U, but we are looking for PO pairs
# borrowed some code from Cassie
FPRs <- mc_sample_simple(Qs,
                         nu = c("PO"),
                         FNRs = c(seq(0.001, 0.20, by = 0.001)))

# and generate the plots
FNR_FPR <- ggplot(FPRs, aes(x = FNR, y = FPR)) +
  geom_point() +
  geom_segment(aes(x = FNR, y = se, xend = FNR, yend = se)) +  # these are basically invisible because they are so small
  facet_wrap(~ numerator) +
  scale_y_continuous(trans = "log10") +
  xlim(0.005, 0.08) + theme(strip.text.x = element_text(size = 12, color = "Black", face = "bold")) +
  theme(axis.text.x = element_text(size=12, color= "Black", face = "bold")) +
  theme(axis.text.y = element_text(size = 12, color = "Black", face = "bold"))

outplotdir = file.path(OUTDIR, "FNR_FPRplot.png")
ggsave(outplotdir, plot=FNR_FPR, dpi=300)

quit()

################################################################################
################################## END SCRIPT ##################################
################################################################################
