#!/usr/bin/env Rscript

## Much of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

suppressPackageStartupMessages(library("tidyverse", quietly=TRUE))
suppressPackageStartupMessages(library("microhaplotopia", quietly=TRUE))
suppressPackageStartupMessages(library("vcfR", quietly=TRUE))
suppressPackageStartupMessages(library("optparse", quietly=TRUE))
suppressPackageStartupMessages(library("reshape2", quietly=TRUE))
suppressPackageStartupMessages(library("stringdist", quietly=TRUE))
suppressPackageStartupMessages(library("viridis", quietly=TRUE))

# for setting working directory when debugging in Rstudio
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# sexID model options - make string for help menu
modelOpts <- c("AFTC", "CDFW")
modelOptsHelp <- paste0("rds file containing logistic regression model for assigning sexID (default = AFTC; Options = ", paste0(modelOpts, collapse = ", "), ")")

## command line option parsing
option_list = list(
  make_option(
    c("-a", "--alleleBalance"), 
    type="double", 
    default=0.35, 
    help="minimum allele balance threshold (default = 0.35)", 
    metavar="alleleBalance"
  ),
  make_option(
    c("-c", "--cutoff"), 
    type="double", 
    default=0.1, 
    help="probability cutoff for sex ID calls with regression model (default = 0.1)", 
    metavar="cutoff"
  ),
  make_option(
    c("-d", "--hapDepth"), 
    type="integer", 
    default=4, 
    help="minimum haplotype sequencing depth (default = 4)", 
    metavar="hapDepth"
  ),
  make_option(
    c("-D", "--totDepth"), 
    type="integer", 
    default=10, 
    help="minimum locus sequencing depth (default = 10)", 
    metavar="totDepth"
  ),
  make_option(
    c("-f", "--finalOut"), 
    type="character", 
    default="haps_2col_final.csv", 
    help="name for final output .csv file (default = haps_2col_final.csv)", 
    metavar="finalOut"
  ),
  make_option(
    c("-g", "--greb1rosaOut"), 
    type="character", 
    default="greb1rosa_all_hapstr.txt", 
    help="output file for greb1rosa haplotype string (default = greb1rosa_all_hapstr.txt)", 
    metavar="greb1rosaOut"
  ),
  make_option(
    c("-G", "--grebInfo"), 
    type="character", 
    default="~/local/src/ca_chinook/example_files/greb1_roha_alleles_reordered_wr.txt", 
    #default="greb1_roha_alleles_reordered_wr.txt", 
    help="info needed to generate greb1rosa haplotype string (default = ~/local/src/ca_chinook/example_files/greb1_roha_alleles_reordered_wr.txt)", 
    metavar="grebInfo"
  ),
  make_option(
    c("-L", "--lfar"), 
    type="character", 
    default="FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds", 
    help="rds file of loci mapped to full genome (default = FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds)", 
    metavar="lfar"
  ),
  make_option(
    c("-m", "--microhaplot"), 
    type="character", 
    default="microhaplot", 
    help="microhaplot directory name (default = microhaplot)", 
    metavar="microhaplot"
  ),
  make_option(
    c("-M", "--sexidmodel"), 
    type="character", 
    default="AFTC", 
    help=modelOptsHelp, 
    metavar="sexidmodel"
  ),
  make_option(
    c("-o", "--output"), 
    type="character", 
    default="output", 
    help="output directory name (default = output)", 
    metavar="output"
  ),
  make_option(
    c("-O", "--reports"), 
    type="character", 
    default="reports", 
    help="reports directory name (default = reports)", 
    metavar="reports"
  ),
  make_option(
    c("-r", "--rosa"), 
    type="character", 
    default=NULL, 
    help="ROSA VCF file (No default; required)", 
    metavar="rosa"
  ),
  make_option(
    c("-R", "--run"), 
    type="character", 
    default=NULL, 
    #default="p133_gt027",
    help="Run number (No default; required)", 
    metavar="run"
  ),
  make_option(
    c("-s", "--sdy"), 
    type="character", 
    default="ordered-read-counts-table.csv", 
    help="sdy read count file (default = ordered-read-counts-table.csv)", 
    metavar="sdy"
  ),
  make_option(
    c("-S", "--snplicon"), 
    type="character", 
    default="FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds", 
    help="rds file of loci mapped to target fastas (default = FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds)", 
    metavar="snplicon"
  ),
  make_option(
    c("-y", "--propMale"), 
    type="double", 
    default=0.002, 
    help="minimum proportion of sdy reads per individual to call a male (default = 0.002)", 
    metavar="propMale"
  ),
  make_option(
    c("-Y", "--minReads"), 
    type="integer", 
    default=10000, 
    help="minimum number of reads per individual for making sex call (default = 10000)", 
    metavar="minReads"
  )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if model options specified correctly
if(!(opt$sexidmodel %in% modelOpts)){
	cat(paste0("Error: '--sexidmodel / -M' must be one of: ", paste(modelOpts, collapse = ", "), "\n\n"))
	#print_help(opt_parser)
	quit(status = 1)
}

if (is.null(opt$rosa)){
  print_help(opt_parser)
  stop("ROSA vcf file must be specified.", call.=FALSE)
}

# capture current working directory in WD
WD <- getwd()

################################################################################
## SET FILE AND PATH NAMES HERE
################################################################################

# name of folder containing microhaplot output
microhaplotDir <- opt$microhaplot

# sex ID file name
sdyReadCounts <- opt$sdy

#your quality and depth-filtered greb1rosa vcf file
vcfFile <- opt$rosa 
#vcfFile <- "CH_gt027_greb1_q20dp5.recode.vcf"

#location of logistic regression model for making sexID calls (distributed as .rds file)
modelFile <- paste0("~/local/src/ca_chinook/example_files/", opt$sexidmodel, "_model.rds")
#modelFile <- "AFTC_model.rds"

# name of folder where output will be written (will be created if doesn't exist)
out <- opt$output
outDir <- file.path(WD, out)
dir.create(outDir, showWarnings = FALSE)

# name of folder where reports will be written (will create if doesn't exist)
reports <- opt$reports
repDir <- file.path(WD, reports)
dir.create(repDir, showWarnings = FALSE)

# file with ROSA haplotype strings
rosaStringFile <- file.path(outDir, opt$greb1rosaOut)

# greb locus information file for making hapstring
grebInfo <- opt$grebInfo

sdyOutFile <- file.path(outDir, "sdy_calls.csv")

################################################################################
## SETTINGS FOR FILTERING
################################################################################

# sdy marker settings
cutoff <- opt$cutoff # probability cutoff for calling sex IDs with logistic regression model. 
propMale <- opt$propMale # minimum proportion of sdy reads per individual to call male
minReads <- opt$minReads # minimum number of total reads per individual to make sex call

# settings for retaining genotype calls
hapDepth <- opt$hapDepth # remove haplotypes with < specified hapDepth
totDepth <- opt$totDepth # remove genotypes with < specified totDepth
alleleBalance <- opt$alleleBalance # heterozygotes only. minimum allowable ratio of read depth for the two haplotypes

################################################################################
## turn ROSA VCF file into hapstring
################################################################################

# read and convert files
locs <- read.table(file=grebInfo, header=T, sep="\t", stringsAsFactors = F) %>% as_tibble()
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

write.table(ptypes, file=rosaStringFile, quote=F, sep="\t", row.names=F)

################################################################################
## read ROSA haplotype string file
################################################################################
rosaHapStr <- read.table(file=rosaStringFile, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>% as_tibble()

canon_rosa_geno_list <- c("EWWEEWEEEEEL", "EWWEEWEEEEEH", "ENNEENEEEEEE", "ENNEENEEEEEL", "ENNEENEEEEEH", "EHHEEHEEEEEL", "EHHEEHEEEEEH", "LNNLLNLLLLLL", "LNNLLNLLLLLH", "HHHHHHHHHHHL", "HHHHHHHHHHHH", "HNNHHNHHHHHL", "HNNHHNHHHHHH")

# match to canon_rosa_geno_list using Levenshtein distance and report values.
rosaHapStr <- rosaHapStr %>% mutate(match_idx = amatch(hapstr, canon_rosa_geno_list, maxDist = 2, method="lv"),
	hapstrMatch = canon_rosa_geno_list[match_idx]
	) %>% mutate(hapstr_dist = stringdist(hapstr, hapstrMatch, method="lv"))

## canonical ROSA phenotype classifications - match using the closest canonical genotype (max = 2 differences; Levenshtein distance)
rosaHapStr <- rosaHapStr %>% mutate(canonical_rosa_pheno = case_when(
	#Early/Early
	hapstrMatch == "EWWEEWEEEEEL" ~ "Winter",
	hapstrMatch == "EWWEEWEEEEEH" ~ "Winter",
	hapstrMatch == "ENNEENEEEEEE" ~ "Spring",
	hapstrMatch == "ENNEENEEEEEL" ~ "Spring",
	hapstrMatch == "ENNEENEEEEEH" ~ "Spring",
	hapstrMatch == "EHHEEHEEEEEL" ~ "Sp-Win",
	hapstrMatch == "EHHEEHEEEEEH" ~ "Sp-Win",

	#Late/Late
	hapstrMatch == "LNNLLNLLLLLL" ~ "Fall",
	hapstrMatch == "LNNLLNLLLLLH" ~ "Fall",

	#Early/Late
	hapstrMatch == "HHHHHHHHHHHL" ~ "Fall-Win",
	hapstrMatch == "HHHHHHHHHHHH" ~ "Fall-Win",
	hapstrMatch == "HNNHHNHHHHHL" ~ "Sp-Fall",
	hapstrMatch == "HNNHHNHHHHHH" ~ "Sp-Fall",

	# all other cases
	.default = "Unknown"
))

################################################################################
## read sex-ID marker file
################################################################################
sexy <- read_csv(file=sdyReadCounts, locale=locale(encoding="latin1")) %>% select(NMFS_DNA_ID, sdy_I183)

sdy <- bind_rows(sexy)

################################################################################
## getting microhaplotypes
################################################################################
#file1 = "FullPanel--target_fastas--target_fasta--rosa_microhap_snplicon.rds"
#file2 = "FullPanel--fullgex_remapped_to_thinned--Otsh_v1.0--lfar_wrap_vgll3six6.rds"

RDS_file1 <- readRDS(opt$snplicon)
RDS_file2 <- readRDS(opt$lfar)

mhp_RDS_file <- full_join(RDS_file1, RDS_file2)

hap <- mhp_RDS_file %>%
  rename(
    indiv.ID = id
  ) %>% 
  dplyr::filter(haplo != "haplo") %>%
  mutate(
    indiv.ID = str_remove(indiv.ID,"_*$") # Remove any trailing underscores from individual IDs
  ) %>%
  mutate(
	source = opt$run # add 'source' field with run number value so script stops crashing when looking for missing individuals
  )

#hap <- read_unfiltered_observed(file.path(WD, microhaplotDir))
locCount <- hap %>% group_by(source) %>% summarise(count = n_distinct(locus))

################################################################################
## handle sex-ID calls
################################################################################
# read model
modAFTC <- readRDS(modelFile)

#calculate reads per individual
readCount <- hap %>% group_by(indiv.ID) %>% summarise(sum_reads = sum(depth))
readCount <- rename(readCount, Indiv = indiv.ID)

#join readcounts to sdy data
sdy <- rename(sdy, Indiv = NMFS_DNA_ID)
sdy_out <- left_join(readCount, sdy, by = "Indiv")

# calculate proportion of reads that come from sdy marker
sdy_out <- sdy_out %>% mutate(sex_marker_read_prop = sdy_I183 / sum_reads)
sdy_out <- sdy_out %>% mutate(sdy_prop_sex = case_when(
		sum_reads >= minReads & sex_marker_read_prop >= propMale ~ "Male",
		sum_reads >= minReads & sex_marker_read_prop < propMale ~ "Female",
		.default = NULL
	)
)

# classify sex based on logistic regression model
sdy_out <- sdy_out %>% mutate(prob.male = predict(modAFTC, newdata = sdy_out, type = "response") )
sdy_out <- sdy_out %>% mutate(sdy_model_sex = case_when(
  sum_reads >= minReads & prob.male <= cutoff ~ "Female",
  sum_reads >= minReads & prob.male >= 1.0-cutoff ~ "Male",
  .default = NULL
))

# write out sex ID calls to file
write_csv(sdy_out, file = sdyOutFile)

# count proportion male and female
proportion_m <- sdy_out %>% summarise(prop_male = sum(sdy_model_sex == "Male", na.rm = TRUE) / n())
proportion_f <- sdy_out %>% summarise(prop_female = sum(sdy_model_sex == "Female", na.rm = TRUE) / n())
proportion_u <- 1.0 - (proportion_m$prop_male + proportion_f$prop_female) # calculate unknown sex
comb_sex <- bind_cols(proportion_m, proportion_f) %>% mutate(prop_unk = proportion_u)
cat("\nSex proportions:\n")
print(comb_sex)
cat("\n\n")

# make scatterplot showing sexID calls by number of reads
sdyPlot<-ggplot(data=sdy_out, aes(x=sum_reads, y=sex_marker_read_prop, color = sdy_model_sex)) +
	geom_point()
suppressMessages(ggsave("sdyPlot.png", path=repDir, dpi=600))

# make histogram showing sexID calls by proportion
sdyHisto<-ggplot(data=sdy_out, aes(x=sex_marker_read_prop, fill=sdy_model_sex)) + 
	geom_histogram(position="identity") + 
	scale_x_continuous(n.breaks = 10)
suppressMessages(ggsave("sdyHist.png", path=repDir, dpi=600))

readsPlot<-ggplot(data=sdy_out, aes(x=reorder(Indiv, sum_reads), y=sum_reads)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=2))
suppressMessages(ggsave("readsPlot.png", path=repDir, dpi=600))

################################################################################
## filtering microhaplotypes
################################################################################
# Filter based on depth and allele balance
hap_fil <- filter_raw_microhap_data(
  hap,
  haplotype_depth = hapDepth,
  total_depth = totDepth, 
  allele_balance = alleleBalance
)
locCount <- hap_fil %>% group_by(source) %>% summarise(count = n_distinct(locus))
#locCount

# Find and drop N/X alleles
nxa <- find_NXAlleles(hap_fil)
write_csv(nxa, file = file.path(repDir, "nxa.csv"))
nxaCount <- nxa %>% group_by(locus) %>% summarise(count = n()) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(nxaCount, file = file.path(repDir, "nxaCount.csv"))

hap_fil_nxa <- nxa %>% 
  select(group, indiv.ID, locus) %>% 
  distinct() %>% 
  anti_join(hap_fil, .)

# Find and drop extra alleles
xtralleles <- find_contaminated_samples(hap_fil_nxa)
write_csv(xtralleles, file = file.path(repDir, "extra_alleles.csv"))
cat(paste0("Samples and loci with extra alleles written to ", file.path(repDir, "extra_alleles.csv\n\n")))

# And am curious about indiv, run and locus effects
indiv <- xtralleles %>% group_by(indiv.ID) %>% summarise(count = n_distinct(locus)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(indiv, file = file.path(repDir, "extra_alleles_individuals.csv"))
cat(paste0("Count of loci per sample impacted by extra alleles written to ", file.path(repDir, "extra_alleles_individuals.csv\n\n")))

indivPlot<-ggplot(data=indiv, aes(x=reorder(indiv.ID, count), y=count)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=2))
suppressMessages(ggsave("extra_alleles_plot.png", path=repDir, dpi=600))

locus <- xtralleles %>% group_by(locus) %>% summarise(count = n_distinct(indiv.ID)) %>% arrange(desc(count)) #modified from Anthony's code to count the number of individuals per locus that were impacted by extra alleles
write_csv(locus, file = file.path(repDir, "extra_alleles_locus.csv"))
cat(paste0("Count of number of individuals with extra alleles per locus written to ", file.path(repDir, "extra_alleles_locus.csv\n\n")))

# remove genotypes with >2 alleles
hap_fil1 <- hap_fil_nxa %>% 
  anti_join(xtralleles)

## Summarize sequencing depth
loc_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "locus") %>% 
  arrange(., n_samples)
write_csv(loc_depth, file = file.path(repDir, "locus_depth_summary.csv"))
cat(paste0("Locus depth summary written to ", file.path(repDir, "locus_depth_summary.csv\n\n")))

ind_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "indiv.ID") %>% 
  arrange(., mean_depth)
write_csv(ind_depth, file = file.path(repDir, "individual_depth_summary.csv"))
cat(paste0("Individual depth written to ", file.path(repDir, "individual_depth_summary.csv\n\n")))

grps <- summarize_data(
  datafile = hap_fil1,
  group_var = "group") %>% 
  arrange(., mean_depth)

## Sequencing Depth Tables
ind_depth_select <- hap_fil1 %>% select(indiv.ID, locus, depth) # grab relevant columns

# table with total depth
ind_depth_table_total <- pivot_wider(
  ind_depth_select,
  names_from = locus,
  values_from = depth, values_fn = list(depth = sum),
)
cat(paste("Depth per locus per individual written to ", file.path(repDir, "totaldepth_per_locus_per_indiv.tsv\n\n")))
write_tsv(ind_depth_table_total,file = file.path(repDir, "totaldepth_per_locus_per_indiv.tsv"))

# plot total depth
totaldepth_long <- melt(ind_depth_table_total)
manual_breaks <- c(0, 10, 25, 50, 100, 500, 1000, Inf)
depth_labels <- c("0-9", "10-24","25-49", "50-99", "100-499", "500-999","1000+")
totaldepth_long$value_discrete <- cut(totaldepth_long$value, 
                                      breaks = manual_breaks, 
                                      labels = depth_labels,
                                      right = FALSE)
totaldepth_long$value_discrete <- factor(totaldepth_long$value_discrete)
totalDepthPlot <- ggplot(totaldepth_long, aes(x = indiv.ID, y = variable, fill = value_discrete)) +
  geom_tile() +
  scale_fill_viridis_d(option = "D") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=1),
        axis.text.y = element_text(size=4))
suppressMessages(ggsave("totalDepth.png", path=repDir, dpi=1200))

# table with haplotype depths
ind_depth_table_haplo <- pivot_wider(
  ind_depth_select,
  names_from = locus,
  values_from = depth, values_fn = function(x) paste(x, collapse = "|"),
)
cat(paste("Depth per haplotype per individual written to ", file.path(repDir, "haplodepth_per_locus_per_indiv.tsv\n\n")))
write_tsv(ind_depth_table_haplo,file = file.path(repDir, "haplodepth_per_locus_per_indiv.tsv"))

# find samples that were removed by filters
missing_samples <- find_missing_samples(hap, hap_fil1)

# remove duplicates if necessary
duplicate_samples <- find_duplicates(hap_fil1)
if (!is.null(duplicate_samples)){
  hap_fil1 <- resolve_duplicate_samples(hap_fil1, resolve = "drop")
}

# Add second allele for homozygotes
hap_final <- add_hom_second_allele(hap_fil1)

# Format for RUBIAS/CKMR
haps_2col <- mhap_transform(
  long_genos = hap_final,
  program = "rubias"
)

# append "_1" and "_2" to the locus names
suffs <- c(1,2)
locs <- colnames(haps_2col)[-1]
addnums <- as_tibble(cbind(locname = locs, suffix = suffs)) %>% 
  mutate(twocol = paste(locs, suffs, sep = "_")) %>% 
  pull(twocol)

haps_2col_final <- haps_2col
names(haps_2col_final) <- c("indiv", addnums)

# add back in missing individuals
haps_2col_final <- haps_2col_final %>% add_row(indiv=missing_samples$indiv.ID)
haps_2col_final <- rename(haps_2col_final, Indiv = indiv) # rename indiv to Indiv for combining tibbles

# add back any mising loci
columnOrder <- c("Indiv", "NC_037099.1:62937268-62937373_1", "NC_037099.1:62937268-62937373_2", "NC_037104.1:55923357-55923657_1", "NC_037104.1:55923357-55923657_2", "NC_037104.1:56552815-56552925_1", "NC_037104.1:56552815-56552925_2", "NC_037104.1:56552952-56553042_1", "NC_037104.1:56552952-56553042_2", "NC_037108.1:73543706-73544006_1", "NC_037108.1:73543706-73544006_2", "NC_037112.1:24542569-24542869_1", "NC_037112.1:24542569-24542869_2", "NC_037130.1:1062935-1063235_1", "NC_037130.1:1062935-1063235_2", "NC_037130.1:864908-865208_1", "NC_037130.1:864908-865208_2", "OkiOts_120255-113_1", "OkiOts_120255-113_2", "Ots_100884-287_1", "Ots_100884-287_2", "Ots_101119-381_1", "Ots_101119-381_2", "Ots_101704-143_1", "Ots_101704-143_2", "Ots_102213-210_1", "Ots_102213-210_2", "Ots_102414-395_1", "Ots_102414-395_2", "Ots_102457-132_1", "Ots_102457-132_2", "Ots_102801-308_1", "Ots_102801-308_2", "Ots_102867-609_1", "Ots_102867-609_2", "Ots_103041-52_1", "Ots_103041-52_2", "Ots_104063-132_1", "Ots_104063-132_2", "Ots_104569-86_1", "Ots_104569-86_2", "Ots_105105-613_1", "Ots_105105-613_2", "Ots_105132-200_1", "Ots_105132-200_2", "Ots_105401-325_1", "Ots_105401-325_2", "Ots_105407-117_1", "Ots_105407-117_2", "Ots_106499-70_1", "Ots_106499-70_2", "Ots_106747-239_1", "Ots_106747-239_2", "Ots_107074-284_1", "Ots_107074-284_2", "Ots_107285-93_1", "Ots_107285-93_2", "Ots_107806-821_1", "Ots_107806-821_2", "Ots_108007-208_1", "Ots_108007-208_2", "Ots_108390-329_1", "Ots_108390-329_2", "Ots_109693-392_1", "Ots_109693-392_2", "Ots_110064-383_1", "Ots_110064-383_2", "Ots_110495-380_1", "Ots_110495-380_2", "Ots_110551-64_1", "Ots_110551-64_2", "Ots_111312-435_1", "Ots_111312-435_2", "Ots_111666-408_1", "Ots_111666-408_2", "Ots_111681-657_1", "Ots_111681-657_2", "Ots_112301-43_1", "Ots_112301-43_2", "Ots_112419-131_1", "Ots_112419-131_2", "Ots_112820-284_1", "Ots_112820-284_2", "Ots_112876-371_1", "Ots_112876-371_2", "Ots_113242-216_1", "Ots_113242-216_2", "Ots_117043-255_1", "Ots_117043-255_2", "Ots_117242-136_1", "Ots_117242-136_2", "Ots_117432-409_1", "Ots_117432-409_2", "Ots_118175-479_1", "Ots_118175-479_2", "Ots_118205-61_1", "Ots_118205-61_2", "Ots_118938-325_1", "Ots_118938-325_2", "Ots_122414-56_1", "Ots_122414-56_2", "Ots_123048-521_1", "Ots_123048-521_2", "Ots_123921-111_1", "Ots_123921-111_2", "Ots_124774-477_1", "Ots_124774-477_2", "Ots_127236-62_1", "Ots_127236-62_2", "Ots_128302-57_1", "Ots_128302-57_2", "Ots_128693-461_1", "Ots_128693-461_2", "Ots_128757-61_1", "Ots_128757-61_2", "Ots_129144-472_1", "Ots_129144-472_2", "Ots_129170-683_1", "Ots_129170-683_2", "Ots_129458-451_1", "Ots_129458-451_2", "Ots_130720-99_1", "Ots_130720-99_2", "Ots_131460-584_1", "Ots_131460-584_2", "Ots_94857-232_1", "Ots_94857-232_2", "Ots_96222-525_1", "Ots_96222-525_2", "Ots_96500-180_1", "Ots_96500-180_2", "Ots_97077-179_1", "Ots_97077-179_2", "Ots_99550-204_1", "Ots_99550-204_2", "Ots_AldB1-122_1", "Ots_AldB1-122_2", "Ots_AldoB4-183_1", "Ots_AldoB4-183_2", "Ots_AsnRS-60_1", "Ots_AsnRS-60_2", "Ots_aspat-196_1", "Ots_aspat-196_2", "Ots_BMP-2-SNP1_1", "Ots_BMP-2-SNP1_2", "Ots_CD59-2_1", "Ots_CD59-2_2", "Ots_CD63_1", "Ots_CD63_2", "Ots_EP-529_1", "Ots_EP-529_2", "Ots_mybp-85_1", "Ots_mybp-85_2", "Ots_myoD-364_1", "Ots_myoD-364_2", "Ots_NAML12_1-SNP1_1", "Ots_NAML12_1-SNP1_2", "Ots_PGK-54_1", "Ots_PGK-54_2", "Ots_Prl2_1", "Ots_Prl2_2", "Ots_S71-336_1", "Ots_S71-336_2", "Ots_SClkF2R2-135_1", "Ots_SClkF2R2-135_2", "Ots_SWS1op-182_1", "Ots_SWS1op-182_2", "Ots_u07-07.161_1", "Ots_u07-07.161_2", "Ots_u07-49.290_1", "Ots_u07-49.290_2", "Ots_u4-92_1", "Ots_u4-92_2", "Ots_unk_526_1", "Ots_unk_526_2", "tag_id_1030_1", "tag_id_1030_2", "tag_id_1079_1", "tag_id_1079_2", "tag_id_1126_1", "tag_id_1126_2", "tag_id_1144_1", "tag_id_1144_2", "tag_id_1191_1", "tag_id_1191_2", "tag_id_120_1", "tag_id_120_2", "tag_id_1243_1", "tag_id_1243_2", "tag_id_1276_1", "tag_id_1276_2", "tag_id_1281_1", "tag_id_1281_2", "tag_id_1363_1", "tag_id_1363_2", "tag_id_1413_1", "tag_id_1413_2", "tag_id_1425_1", "tag_id_1425_2", "tag_id_1470_1", "tag_id_1470_2", "tag_id_1551_1", "tag_id_1551_2", "tag_id_1554_1", "tag_id_1554_2", "tag_id_1692_1", "tag_id_1692_2", "tag_id_1733_1", "tag_id_1733_2", "tag_id_186_1", "tag_id_186_2", "tag_id_1872_1", "tag_id_1872_2", "tag_id_2_1016_1", "tag_id_2_1016_2", "tag_id_2_1158_1", "tag_id_2_1158_2", "tag_id_2_123_1", "tag_id_2_123_2", "tag_id_2_1268_1", "tag_id_2_1268_2", "tag_id_2_136_1", "tag_id_2_136_2", "tag_id_2_1382_1", "tag_id_2_1382_2", "tag_id_2_1539_1", "tag_id_2_1539_2", "tag_id_2_1579_1", "tag_id_2_1579_2", "tag_id_2_1586_1", "tag_id_2_1586_2", "tag_id_2_1693_1", "tag_id_2_1693_2", "tag_id_2_188_1", "tag_id_2_188_2", "tag_id_2_1887_1", "tag_id_2_1887_2", "tag_id_2_20_1", "tag_id_2_20_2", "tag_id_2_206_1", "tag_id_2_206_2", "tag_id_2_2222_1", "tag_id_2_2222_2", "tag_id_2_234_1", "tag_id_2_234_2", "tag_id_2_2632_1", "tag_id_2_2632_2", "tag_id_2_2741_1", "tag_id_2_2741_2", "tag_id_2_2787_1", "tag_id_2_2787_2", "tag_id_2_284_1", "tag_id_2_284_2", "tag_id_2_3026_1", "tag_id_2_3026_2", "tag_id_2_3094_1", "tag_id_2_3094_2", "tag_id_2_311_1", "tag_id_2_311_2", "tag_id_2_321_1", "tag_id_2_321_2", "tag_id_2_332_1", "tag_id_2_332_2", "tag_id_2_3452_1", "tag_id_2_3452_2", "tag_id_2_3471_1", "tag_id_2_3471_2", "tag_id_2_40_1", "tag_id_2_40_2", "tag_id_2_414_1", "tag_id_2_414_2", "tag_id_2_419_1", "tag_id_2_419_2", "tag_id_2_487_1", "tag_id_2_487_2", "tag_id_2_502_1", "tag_id_2_502_2", "tag_id_2_58_1", "tag_id_2_58_2", "tag_id_2_633_1", "tag_id_2_633_2", "tag_id_2_661_1", "tag_id_2_661_2", "tag_id_2_694_1", "tag_id_2_694_2", "tag_id_2_700_1", "tag_id_2_700_2", "tag_id_2_705_1", "tag_id_2_705_2", "tag_id_2_749_1", "tag_id_2_749_2", "tag_id_2_786_1", "tag_id_2_786_2", "tag_id_2_855_1", "tag_id_2_855_2", "tag_id_2_859_1", "tag_id_2_859_2", "tag_id_2_9_1", "tag_id_2_9_2", "tag_id_2_911_1", "tag_id_2_911_2", "tag_id_2_935_1", "tag_id_2_935_2", "tag_id_2_939_1", "tag_id_2_939_2", "tag_id_2_953_1", "tag_id_2_953_2", "tag_id_2_978_1", "tag_id_2_978_2", "tag_id_2_98_1", "tag_id_2_98_2", "tag_id_235_1", "tag_id_235_2", "tag_id_251_1", "tag_id_251_2", "tag_id_275_1", "tag_id_275_2", "tag_id_278_1", "tag_id_278_2", "tag_id_282_1", "tag_id_282_2", "tag_id_3194_1", "tag_id_3194_2", "tag_id_32_1", "tag_id_32_2", "tag_id_3221_1", "tag_id_3221_2", "tag_id_381_1", "tag_id_381_2", "tag_id_384_1", "tag_id_384_2", "tag_id_3920_1", "tag_id_3920_2", "tag_id_423_1", "tag_id_423_2", "tag_id_425_1", "tag_id_425_2", "tag_id_427_1", "tag_id_427_2", "tag_id_430_1", "tag_id_430_2", "tag_id_481_1", "tag_id_481_2", "tag_id_4969_1", "tag_id_4969_2", "tag_id_542_1", "tag_id_542_2", "tag_id_5617_1", "tag_id_5617_2", "tag_id_5720_1", "tag_id_5720_2", "tag_id_600_1", "tag_id_600_2", "tag_id_603_1", "tag_id_603_2", "tag_id_650_1", "tag_id_650_2", "tag_id_664_1", "tag_id_664_2", "tag_id_669_1", "tag_id_669_2", "tag_id_684_1", "tag_id_684_2", "tag_id_695_1", "tag_id_695_2", "tag_id_70_1", "tag_id_70_2", "tag_id_716_1", "tag_id_716_2", "tag_id_744_1", "tag_id_744_2", "tag_id_757_1", "tag_id_757_2", "tag_id_773_1", "tag_id_773_2", "tag_id_787_1", "tag_id_787_2", "tag_id_819_1", "tag_id_819_2", "tag_id_826_1", "tag_id_826_2", "tag_id_871_1", "tag_id_871_2", "tag_id_945_1", "tag_id_945_2", "tag_id_999_1", "tag_id_999_2")

# insert missing locus columns and populate with 'NA' values
haps_2col_final <- haps_2col_final %>%
	bind_cols(setNames(rep(list(NA), length(columnOrder)), columnOrder)[!columnOrder %in% names(.)]) %>%
	select(all_of(columnOrder), everything())

# calculate and plot missing data per individual now that missing individuals and loci have been added back
himiss <- haps_2col_final %>% mutate(n_miss = rowSums(is.na(haps_2col_final)/2)) # count missing loci and add as n_miss field
himiss <- select(himiss, Indiv, n_miss) # select only Indiv and n_miss fields
hm<-ggplot(data=himiss, aes(x=reorder(Indiv, n_miss), y=n_miss)) +
	geom_bar(stat="identity")
suppressMessages(ggsave("missingData.png", path=repDir, dpi=600))

# calculate percent of microhaps genotyped per individual and add to the tibble as column 'percMicroHap'
percentLoci <- unlist(as_tibble(((ncol(haps_2col_final)-1)-rowSums(is.na(haps_2col_final)))/(ncol(haps_2col_final)-1))*100)
haps_2col_final <- haps_2col_final %>% mutate(percMicroHap = percentLoci)

## plot depth vs. percMicroHap
percMicro <- select(haps_2col_final, Indiv, percMicroHap) # get relevant columns
seqSucc <- left_join(readCount, percMicro, by = c("Indiv" = "Indiv")) # combine relevant columns

# make zoomed-in plot
cat("If one of the next few lines is a warning about 'rows containing missing values or values outside the scale range,' just ignore it.\n\n")
succPlotZoomed <- ggplot(seqSucc, aes(y=percMicroHap, x=sum_reads)) +
	# geom_smooth(method = "lm") + 
	geom_point() + xlim(c(0,80000)) + theme_minimal() + xlab("Total Reads") + ylab("Percent Success") + ggtitle("Success Rate vs. Total Reads (Zoomed)")
suppressMessages(ggsave("success_rate_v_total_reads_zoomed.png", path=repDir, dpi=600))

# plot all individuals
succPlot <- ggplot(seqSucc, aes(y=percMicroHap, x=sum_reads)) +
	# geom_smooth(method = "lm") + 
	geom_point() + theme_minimal() + xlab("Total Reads") + ylab("Percent Success") + ggtitle("Success Rate vs. Total Reads (All)")
suppressMessages(ggsave("success_rate_v_total_reads_all.png", path=repDir, dpi=600))


# add the sex and ROSA columns
haps_2col_final <- full_join(rosaHapStr,haps_2col_final, by="Indiv") # add ROSA hap string
haps_2col_final <- full_join(sdy_out,haps_2col_final, by="Indiv") # add ROSA hap string

# add extra allele counts
xtraAlleles <- indiv %>% rename(Indiv = indiv.ID)
haps_2col_final <- full_join(xtraAlleles,haps_2col_final, by="Indiv")
haps_2col_final <- haps_2col_final %>% mutate(count = replace_na(count, 0))

# remove extra columns added by sdy_out
haps_2col_final <- haps_2col_final %>% select(!c(sdy_I183, sum_reads, sex_marker_read_prop, sdy_prop_sex, prob.male))

haps_2col_final <- arrange(haps_2col_final, Indiv) # sort by Indiv column

haps_2col_final <- rename(haps_2col_final, indiv = Indiv) # rename indiv to Indiv for combining tibbles

## sort alleles alphabetically per locus
cat("\nSorting alleles alphabetically per locus per individual. This part takes longer to run than it probably should...\n\n")
n <- names(haps_2col_final) # get column names
r <- c("count", "indiv", "sdy_model_sex", "hapstr", "match_idx", "hapstrMatch", "hapstr_dist", "canonical_rosa_pheno", "percMicroHap") # vector of columns to remove
result <- setdiff(n, r) # remove columns
result2 <- str_sub(result, end = -3) # remove _1 and _2 from end of allele names
loci <- unique(result2) # reduce to only locus names

# calculate proportion of loci that had extra alleles
haps_2col_final <- haps_2col_final %>% mutate(perc_Xtra = (count/length(loci))*100) # calculate percent loci that had extra alleles
haps_2col_final <- arrange(haps_2col_final, indiv) # sort by indiv column
haps_2col_final <- haps_2col_final %>% select(!c(count)) # remove count column

# do the allele sort
for (locus in loci) {
  allele1 = paste0(locus, "_1")
  allele2 = paste0(locus, "_2")
  cols <- haps_2col_final %>% select(any_of(c(allele1, allele2)))
  
  cols2 <- cols |> 
    rowwise() |>
    mutate(sorted_vals = list(sort(c_across(everything())))) %>%
    mutate(c1_sorted = sorted_vals[1],
           c2_sorted = sorted_vals[2]
    ) %>%
    select(c1_sorted, c2_sorted) %>%
    rename(
      "{allele1}" := c1_sorted, # this is the dumbest way of accessing a variable to rename something I've ever seen.
      "{allele2}" := c2_sorted
    )
  #print(columns2, n=100)
  
  haps_2col_final[[allele1]] <- cols2[[allele1]]
  haps_2col_final[[allele2]] <- cols2[[allele2]]
  
}
cat("Done sorting.\n\n")

# prepare the final column order that we want
columnOrder <- columnOrder[-1] # drop first element of previous list
frontColumns <- c("indiv", "sdy_model_sex", "hapstr", "hapstr_dist", "canonical_rosa_pheno", "percMicroHap", "perc_Xtra") # get columns that we want appearing in first few columns of the output table
columnOrder <- c(frontColumns, columnOrder) # combine the two lists

# do the sorting and insert missing columns as 'NA' values
haps_2col_final <- haps_2col_final %>%
	select(all_of(columnOrder))

# write final genotype file
write_csv(haps_2col_final, file=file.path(outDir, opt$finalOut))

quit()

