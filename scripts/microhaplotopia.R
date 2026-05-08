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
    c("-M", "--modaftc"), 
    type="character", 
    default="~/local/src/ca_chinook/example_files/AFTC_model.rds", 
    #default="AFTC_model.rds", 
    help=".rds file containing logistic regression model for assigning sexID (default = ~/local/src/ca_chinook/example_files/AFTC_model.rds)", 
    metavar="modaftc"
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

#logistic regression model for making sexID calls (distributed as .rds file)
modelFile <- opt$modaftc 
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
	) %>% mutate(lv_dist = stringdist(hapstr, hapstrMatch, method="lv"))

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
sdy_out %>% summarise(prop_male = sum(sdy_model_sex == "Male", na.rm = TRUE) / n())
sdy_out %>% summarise(prop_female = sum(sdy_model_sex == "Female", na.rm = TRUE) / n())

# make scatterplot showing sexID calls by number of reads
sdyPlot<-ggplot(data=sdy_out, aes(x=sum_reads, y=sex_marker_read_prop, color = sdy_model_sex)) +
	geom_point()
ggsave("sdyPlot.png", path=repDir, dpi=600)

# make histogram showing sexID calls by proportion
sdyHisto<-ggplot(data=sdy_out, aes(x=sex_marker_read_prop, fill=sdy_model_sex)) + 
	geom_histogram(position="identity") + 
	scale_x_continuous(n.breaks = 10)
ggsave("sdyHist.png", path=repDir, dpi=600)

readsPlot<-ggplot(data=sdy_out, aes(x=reorder(Indiv, sum_reads), y=sum_reads)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=2))
ggsave("readsPlot.png", path=repDir, dpi=600)


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
locCount

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
print(paste("Samples/loci with extra alleles written to ", file.path(repDir, "extra_alleles.csv")))

# And am curious about indiv, run and locus effects
indiv <- xtralleles %>% group_by(indiv.ID) %>% summarise(count = n_distinct(locus)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(indiv, file = file.path(repDir, "extra_alleles_individuals.csv"))
print(paste("Count of loci per sample impacted by extra alleles written to ", file.path(repDir, "extra_alleles_individuals.csv")))

indivPlot<-ggplot(data=indiv, aes(x=reorder(indiv.ID, count), y=count)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=2))
ggsave("extra_alleles_plot.png", path=repDir, dpi=600)

locus <- xtralleles %>% group_by(locus) %>% summarise(count = n_distinct(indiv.ID)) %>% arrange(desc(count)) #modified from Anthony's code to count the number of individuals per locus that were impacted by extra alleles
write_csv(locus, file = file.path(repDir, "extra_alleles_locus.csv"))
print(paste("Count of number of individuals with extra alleles per locus written to ", file.path(repDir, "extra_alleles_locus.csv")))

# remove genotypes with >2 alleles
hap_fil1 <- hap_fil_nxa %>% 
  anti_join(xtralleles)

# these statistics can be examined before/after subsequent filtering steps
# Check for chinookie
#hap_fil1 %>% filter(str_detect(locus,"^OkiOts")) %>%  pull(haplo) %>% table()

## Summarize sequencing depth
loc_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "locus") %>% 
  arrange(., n_samples)
write_csv(loc_depth, file = file.path(repDir, "locus_depth_summary.csv"))
print(paste("Locus depth summary written to ", file.path(repDir, "locus_depth_summary.csv")))

ind_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "indiv.ID") %>% 
  arrange(., mean_depth)
write_csv(ind_depth, file = file.path(repDir, "individual_depth_summary.csv"))
print(paste("Individual depth written to ", file.path(repDir, "individual_depth_summary.csv")))

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
print(paste("Depth per locus per individual written to ", file.path(repDir, "totaldepth_per_locus_per_indiv.tsv")))
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
ggsave("totalDepth.png", path=repDir, dpi=1200)

# table with haplotype depths
ind_depth_table_haplo <- pivot_wider(
  ind_depth_select,
  names_from = locus,
  values_from = depth, values_fn = function(x) paste(x, collapse = "|"),
)
print(paste("Depth per haplotype per individual written to ", file.path(repDir, "haplodepth_per_locus_per_indiv.tsv")))
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

# calculate and plot missing data now that missing individuals have been added back
himiss <- haps_2col_final %>% mutate(n_miss = rowSums(is.na(haps_2col_final)/2)) # count missing loci and add as n_miss field
himiss <- select(himiss, Indiv, n_miss) # select only Indiv and n_miss fields
hm<-ggplot(data=himiss, aes(x=reorder(Indiv, n_miss), y=n_miss)) +
  geom_bar(stat="identity")
ggsave("missingData.png", path=repDir, dpi=600)

# calculate percent of microhaps genotyped per individual
percentLoci <- unlist(as_tibble(((ncol(haps_2col_final)-1)-rowSums(is.na(haps_2col_final)))/(ncol(haps_2col_final)-1))*100)
# add to tibble and move to first column after individual name
haps_2col_final <- haps_2col_final %>% mutate(percMicroHap = percentLoci)
haps_2col_final <- haps_2col_final %>% relocate(percMicroHap, .after = Indiv)

## plot depth vs. percMicroHap
percMicro <- select(haps_2col_final, Indiv, percMicroHap) # get relevant columns
seqSucc <- left_join(readCount, percMicro, by = c("Indiv" = "Indiv")) # combine relevant columns

# make zoomed-in plot
print("If the next line is a warning about missing values just ignore it.")
succPlotZoomed <- ggplot(seqSucc, aes(y=percMicroHap, x=sum_reads)) +
  # geom_smooth(method = "lm") + 
  geom_point() + xlim(c(0,80000)) + theme_minimal() + xlab("Total Reads") + ylab("Percent Success") + ggtitle("Success Rate vs. Total Reads (Zoomed)")
ggsave("success_rate_v_total_reads_zoomed.png", path=repDir, dpi=600)

# plot all individuals
succPlot <- ggplot(seqSucc, aes(y=percMicroHap, x=sum_reads)) +
  # geom_smooth(method = "lm") + 
  geom_point() + theme_minimal() + xlab("Total Reads") + ylab("Percent Success") + ggtitle("Success Rate vs. Total Reads (All)")
ggsave("success_rate_v_total_reads_all.png", path=repDir, dpi=600)



# add the sex and ROSA columns
haps_2col_final <- full_join(rosaHapStr,haps_2col_final, by="Indiv") # add ROSA hap string
haps_2col_final <- full_join(sdy_out,haps_2col_final, by="Indiv") # add ROSA hap string

# add extra allele counts
xtraAlleles <- indiv %>% rename(Indiv = indiv.ID)
haps_2col_final <- full_join(xtraAlleles,haps_2col_final, by="Indiv")
haps_2col_final <- haps_2col_final %>% mutate(count = replace_na(count, 0))

# remove extra columns added by sdy_out
haps_2col_final <- haps_2col_final %>% select(!c(sdy_I183))
haps_2col_final <- haps_2col_final %>% select(!c(sum_reads))
haps_2col_final <- haps_2col_final %>% select(!c(sex_marker_read_prop))
haps_2col_final <- haps_2col_final %>% select(!c(sdy_prop_sex))
haps_2col_final <- haps_2col_final %>% select(!c(prob.male))

haps_2col_final <- arrange(haps_2col_final, Indiv) # sort by Indiv column

haps_2col_final <- rename(haps_2col_final, indiv = Indiv) # rename indiv to Indiv for combining tibbles

## sort alleles alphabetically per locus
# I hate R, and especially tidyverse.
print("Sorting alleles alphabetically per locus per individual. This part takes longer to run than it probably should.")
n <- names(haps_2col_final) # get column names
r <- c("count", "indiv", "sdy_model_sex", "hapstr", "match_idx", "hapstrMatch", "lv_dist", "canonical_rosa_pheno", "percMicroHap") # vector of columns to remove
result <- setdiff(n, r) # remove columns
result2 <- str_sub(result, end = -3) # remove _1 and _2 from end of allele names
loci <- unique(result2) # reduce to only locus names

# calculate proportion of loci that had extra alleles
haps_2col_final <- haps_2col_final %>% mutate(perc_Xtra = (count/length(loci))*100) # calculate percent loci that had extra alleles
haps_2col_final <- haps_2col_final %>% relocate(perc_Xtra, .after = percMicroHap) # move location
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
print("Done sorting.")

# write final genotype file
write_csv(haps_2col_final, file=file.path(outDir, opt$finalOut))

quit()

