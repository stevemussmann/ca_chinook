#!/usr/bin/env Rscript

## Much of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

suppressPackageStartupMessages(library("tidyverse", quietly=TRUE))
suppressPackageStartupMessages(library("microhaplotopia", quietly=TRUE))
suppressPackageStartupMessages(library("vcfR", quietly=TRUE))
suppressPackageStartupMessages(library("optparse", quietly=TRUE))

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
#  ),
#  make_option(
#    c("-x", "--maxFemale"), 
#    type="integer", 
#    default=2, 
#    help="maximum sdy reads allowed to call female (default = 2)", 
#    metavar="maxFemale"
#  ),
#  make_option(
#    c("-X", "--maxUnknown"), 
#    type="integer", 
#    default=5, 
#    help="maximum number of sdy reads for unknown (default = 5)", 
#    metavar="maxUnknown"
  )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rosa)){
  print_help(opt_parser)
  stop("ROSA vcf file must be specified.", call.=FALSE)
}

# set working directory in Rstudio and capture path in WD
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
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
propMale <- opt$propMale # minimum proportion of sdy reads per individual to call male
minReads <- opt$minReads # minimum number of total reads per individual to make sex call
#maxFemale <- opt$maxFemale #maximum number of sdy reads allowed to call female
#maxUnknown <- opt$maxUnknown #maximum number of sdy reads for unknown. Values > maxFemale and <= maxUnknown will be called "Unk" sex. Values > maxUnknown will be called male.

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

## old way
#rosaHapStr <- rosaHapStr %>% mutate(rosa_pheno = case_when(
#	str_extract(hapstr, "^.{1}") == "E" ~ "sp",
#	str_extract(hapstr, "^.{1}") == "H" ~ "sp-fall",
#	str_extract(hapstr, "^.{1}") == "L" ~ "fall",
#	str_extract(hapstr, "^.{1}") == "?" ~ "unk"
#))

## canonical ROSA haplotype classifications
rosaHapStr <- rosaHapStr %>% mutate(canonical_rosa_pheno = case_when(
	#Early/Early
	hapstr == "EWWEEWEEEEEL" ~ "Winter",
	hapstr == "EWWEEWEEEEEH" ~ "Winter",
	hapstr == "ENNEENEEEEEE" ~ "Spring",
	hapstr == "ENNEENEEEEEL" ~ "Spring",
	hapstr == "ENNEENEEEEEH" ~ "Spring",
	hapstr == "EHHEEHEEEEEL" ~ "Sp-Win",
	hapstr == "EHHEEHEEEEEH" ~ "Sp-Win",

	#Late/Late
	hapstr == "LNNLLNLLLLLL" ~ "Fall",
	hapstr == "LNNLLNLLLLLH" ~ "Fall",

	#Early/Late
	hapstr == "HHHHHHHHHHHL" ~ "Fall-Win",
	hapstr == "HHHHHHHHHHHH" ~ "Fall-Win",
	hapstr == "HNNHHNHHHHHL" ~ "Sp-Fall",
	hapstr == "HNNHHNHHHHHH" ~ "Sp-Fall",

	# all other cases
	.default = "Unknown"
))

################################################################################
## read sex-ID marker file
################################################################################
sexy <- read_csv(file=sdyReadCounts, locale=locale(encoding="latin1")) %>% select(NMFS_DNA_ID, sdy_I183)

sdy <- bind_rows(sexy)

#sdy_out <- sdy %>% mutate(sdy_sex = case_when(
#    sdy_I183 <= maxFemale ~ "Female",
#    sdy_I183 > maxFemale & sdy_I183 <= maxUnknown ~ "Unk",
#    sdy_I183 > maxUnknown ~ "Male"
#  )) %>% 
#  relocate(sdy_sex, .after = sdy_I183)
#sdy_out <- rename(sdy_out, Indiv = NMFS_DNA_ID) # rename NMFS_DNA_ID to Indiv

# write sdy output
#write_csv(sdy_out, file = sdyOutFile)

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
#calculate reads per individual
readCount <- hap %>% group_by(indiv.ID) %>% summarise(sum_reads = sum(depth))
readCount <- rename(readCount, Indiv = indiv.ID)

#join readcounts to sdy data
sdy <- rename(sdy, Indiv = NMFS_DNA_ID)
sdy_out <- left_join(readCount, sdy, by = "Indiv")

#calculate proportion of reads that come from sdy marker
sdy_out <- sdy_out %>% mutate(prop = sdy_I183 / sum_reads)
sdy_out <- sdy_out %>% mutate(sdy_sex = case_when(
		sum_reads >= minReads & prop >= propMale ~ "Male",
		sum_reads >= minReads & prop < propMale ~ "Female",
		.default = NULL
	)
)

# write out sex ID calls to file
write_csv(sdy_out, file = sdyOutFile)

# count proportion male and female
sdy_out %>% summarise(prop_male = sum(sdy_sex == "Male", na.rm = TRUE) / n())
sdy_out %>% summarise(prop_female = sum(sdy_sex == "Female", na.rm = TRUE) / n())

# make scatterplot showing sexID calls by number of reads
sdyPlot<-ggplot(data=sdy_out, aes(x=sum_reads, y=prop, color = sdy_sex)) +
	geom_point()
ggsave("sdyPlot.png", path=repDir)

# make histogram showing sexID calls by proportion
sdyHisto<-ggplot(data=sdy_out, aes(x=prop, fill=sdy_sex)) + 
	geom_histogram(position="identity") + 
	scale_x_continuous(n.breaks = 10)
ggsave("sdyHist.png", path=repDir)


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

# Find and drop extra aleles
xtralleles <- find_contaminated_samples(hap_fil_nxa)
write_csv(xtralleles, file = file.path(repDir, "extra_alleles.csv"))
print(paste("Samples/loci with extra alleles written to ", file.path(repDir, "extra_alleles.csv")))

#And am curious about indiv, run and locus effects
indiv <- xtralleles %>% group_by(indiv.ID) %>% summarise(count = n_distinct(locus)) #modified from Anthony's code to count the number of loci impacted by extra alleles per individual
write_csv(indiv, file = file.path(repDir, "extra_alleles_individuals.csv"))
print(paste("Count of loci per sample impacted by extra alleles written to ", file.path(repDir, "extra_alleles_individuals.csv")))

source <- xtralleles %>% group_by(source) %>% summarise(count = n_distinct(indiv.ID)) #modified from Anthony's code to count the number of individuals per input file that were impacted by extra alleles
write_csv(source, file = file.path(repDir, "extra_alleles_source.csv"))
print(paste("Count of individuals per input file impacted by extra alleles written to ", file.path(repDir, "extra_alleles_source.csv")))

locus <- xtralleles %>% group_by(locus) %>% summarise(count = n_distinct(indiv.ID)) %>% arrange(desc(count)) #modified from Anthony's code to count the number of individuals per locus that were impacted by extra alleles
write_csv(locus, file = file.path(repDir, "extra_alleles_locus.csv"))
print(paste("Count of number of individuals with extra alleles per locus written to ", file.path(repDir, "extra_alleles_locus.csv")))

# remove genotypes with >2 alleles
hap_fil1 <- hap_fil_nxa %>% 
  anti_join(xtralleles)

# these statistics can be examined before/after subsequent filtering steps
# Check for chinookie
#hap_fil1 %>% filter(str_detect(locus,"^OkiOts")) %>%  pull(haplo) %>% table()

#Summarize
loc_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "locus") %>% 
  arrange(., n_samples)
write_csv(loc_depth, file = file.path(repDir, "locus_depth.csv"))
print(paste("Locus depth written to ", file.path(repDir, "locus_depth.csv")))


ind_depth <- summarize_data(
  datafile = hap_fil1,
  group_var = "indiv.ID") %>% 
  arrange(., mean_depth)
write_csv(ind_depth, file = file.path(repDir, "individual_depth.csv"))
print(paste("Locus depth written to ", file.path(repDir, "individual_depth.csv")))

grps <- summarize_data(
  datafile = hap_fil1,
  group_var = "group") %>% 
  arrange(., mean_depth)

# Plot missing data
himiss <- calculate_missing_data(hap_fil1)
hm<-ggplot(data=himiss, aes(x=reorder(indiv.ID, n_miss), y=n_miss)) +
  geom_bar(stat="identity")
ggsave("missingData.png", path=repDir)

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

# calculate percent of microhaps genotyped per individual
percentLoci <- unlist(as_tibble(((ncol(haps_2col_final)-1)-rowSums(is.na(haps_2col_final)))/(ncol(haps_2col_final)-1))*100)
# add to tibble and move to first column after individual name
haps_2col_final <- haps_2col_final %>% mutate(percMicroHap = percentLoci)
haps_2col_final <- haps_2col_final %>% relocate(percMicroHap, .after = Indiv)

# add the sex and ROSA columns
haps_2col_final <- full_join(rosaHapStr,haps_2col_final, by="Indiv") # add ROSA hap string
haps_2col_final <- full_join(sdy_out,haps_2col_final, by="Indiv") # add ROSA hap string
haps_2col_final <- haps_2col_final %>% select(!c(sdy_I183)) # remove extra columns added by sdy_out
haps_2col_final <- haps_2col_final %>% select(!c(sum_reads)) # remove extra columns added by sdy_out
haps_2col_final <- haps_2col_final %>% select(!c(prop)) # remove extra columns added by sdy_out

haps_2col_final <- arrange(haps_2col_final, Indiv) # sort by Indiv column

haps_2col_final <- rename(haps_2col_final, indiv = Indiv) # rename indiv to Indiv for combining tibbles

## sort alleles alphabetically per locus
# I hate R, and especially tidyverse.
print("Sorting alleles alphabetically per locus per individual. This part takes longer to run than it probably should.")
n <- names(haps_2col_final) # get column names
r <- c("indiv", "sdy_sex", "hapstr", "canonical_rosa_pheno", "percMicroHap") # vector of columns to remove
result <- setdiff(n, r) # remove columns
result2 <- str_sub(result, end = -3) # remove _1 and _2 from end of allele names
loci <- unique(result2) # reduce to only locus names

for (locus in loci) {
  allele1 = paste0(locus, "_1")
  allele2 = paste0(locus, "_2")
  cols <- haps_2col_final %>% select(any_of(c(allele1, allele2)))
  
  cols2 <- cols |> # what even is this '|>' nonsense?
    rowwise() |>
    mutate(sorted_vals = list(sort(c_across(everything())))) %>% # this %>% is nonsense too.
    mutate(c1_sorted = sorted_vals[1],
           c2_sorted = sorted_vals[2]
    ) %>%
    select(c1_sorted, c2_sorted) %>%
    rename(
      "{allele1}" := c1_sorted, # this is the dumbest way of accessing a variable to rename something I've ever seen.
      "{allele2}" := c2_sorted
    )
  #print(columns2, n=100)
  
  haps_2col_final[[allele1]] <- cols2[[allele1]] # and now I have to access the variable a different way.
  haps_2col_final[[allele2]] <- cols2[[allele2]] # why is it so hard to access a column name by variable?
  
}
print("Done sorting.")

# write final genotype file
write_csv(haps_2col_final, file=file.path(outDir, opt$finalOut))

quit()

