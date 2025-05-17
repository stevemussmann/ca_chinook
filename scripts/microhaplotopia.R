#!/usr/bin/env Rscript

## Much of the code in this script was modified from code provided by 
## Anthony Clemento at NOAA

library("tidyverse")
library("microhaplotopia")

# set working directory in Rstudio and capture path in WD
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
WD <- getwd()

################################################################################
## SET FILE AND PATH NAMES HERE
################################################################################

# sex ID file name
sdyReadCounts <- "ordered-read-counts-table.csv"

# file with ROSA haplotype strings (output of makeROSAstring.R)
rosaStringFile <- "greb1rosa_all_16MAY2025.txt"

# name of file with new column names
#panelInfo <- "Calif-Chinook-Amplicon-Panel-Information.csv"

# name of folder containing microhaplot output
microhaplotDir <- "microhaplot"

# name of folder where output will be written (will be created if doesn't exist)
out <- "output"

# name of folder where reports will be written (will create if doesn't exist)
reports <- "reports"

################################################################################
## SETTINGS FOR FILTERING
################################################################################

# sdy marker settings
maxFemale <- 2 #maximum number of sdy reads allowed to call female
maxUnknown <- 5 #maximum number of sdy reads for unknown. Values > maxFemale and <= maxUnknown will be called "Unk" sex. Values > maxUnknown will be called male.

# settings for retaining genotype calls
hapDepth <- 4 # remove haplotypes with < specified hapDepth
totDepth <- 8 # remove genotypes with < specified totDepth
alleleBalance <- 0.35 # heterozygotes only. minimum allowable ratio of read depth for the two haplotypes

################################################################################
## CREATE DIRECTORIES
################################################################################
# create directory for output files
outDir <- file.path(WD, out)
dir.create(outDir, showWarnings = FALSE)

# create directory for reports and summary files
repDir <- file.path(WD, reports)
dir.create(repDir, showWarnings = FALSE)

################################################################################
## read ROSA haplotype string file
################################################################################
rosaHapStr <- read.table(file=rosaStringFile, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>% as_tibble()

################################################################################
## read panel info (contains new column names)
################################################################################
#panelInfoStr <- read.table(file=panelInfo, header=TRUE, sep=",", stringsAsFactors = FALSE) %>% as_tibble()


################################################################################
## sex-ID marker
################################################################################
sexy <- read_csv(file=sdyReadCounts, locale=locale(encoding="latin1")) %>% select(NMFS_DNA_ID, sdy_I183)

sdy <- bind_rows(sexy)

sdy_out <- sdy %>% mutate(sdy_sex = case_when(
    sdy_I183 <= maxFemale ~ "Female",
    sdy_I183 > maxFemale & sdy_I183 <= maxUnknown ~ "Unk",
    sdy_I183 > maxUnknown ~ "Male"
  )) %>% 
  relocate(sdy_sex, .after = sdy_I183)
sdy_out <- rename(sdy_out, Indiv = NMFS_DNA_ID) # rename NMFS_DNA_ID to Indiv
# write sdy output
write_csv(sdy_out, file = file.path(outDir, "sdy_calls.csv"))
################################################################################


################################################################################
## getting microhaplotypes
################################################################################
hap <- read_unfiltered_observed(file.path(WD, microhaplotDir))
locCount <- hap %>% group_by(source) %>% summarise(count = n_distinct(locus))


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
#view(nxaCount)
write_csv(nxaCount, file = file.path(repDir, "nxaCount.csv"))

hap_fil_nxa <- nxa %>% 
  select(group, indiv.ID, locus) %>% 
  distinct() %>% 
  anti_join(hap_fil, .)

locCount <- hap_fil_nxa %>% group_by(source) %>% summarise(count = n_distinct(locus))
locCount

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
write_csv(source, file = file.path(repDir, "extra_alleles_locus.csv"))
print(paste("Count of number of individuals with extra alleles per locus written to ", file.path(repDir, "extra_alleles_locus.csv")))


# remove genotypes with >2 alleles
hap_fil1 <- hap_fil_nxa %>% 
  anti_join(xtralleles)
locCount <- hap_fil1 %>% group_by(source) %>% summarise(count = n_distinct(locus))
locCount

# these statistics can be examined before/after subsequent filtering steps
# Check for chinookie
hap_fil1 %>% filter(str_detect(locus,"^OkiOts")) %>%  pull(haplo) %>% table()

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
grps


# Plot missing data
himiss <- calculate_missing_data(hap_fil1)
hm<-ggplot(data=himiss, aes(x=reorder(indiv.ID, n_miss), y=n_miss)) +
  geom_bar(stat="identity")
hm
ggsave("missingData.png", path=repDir)


# find samples that were removed by filters
missing_samples <- find_missing_samples(hap, hap_fil1)
#table(missing_samples$group)

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

###CODE THAT DID NOT WORK FOR RENAMING COLUMNS
# get list of unique values in column names, without first column (indiv)
#loci <- unique(colnames(haps_2col)[-1])

# remove loci from panelInfoStr that are not in our tibble
#panelInfoStr <- panelInfoStr %>% filter(OtherName %in% loci)

#stuff <- haps_2col %>% rename_with(~ panelInfoStr$AmpliconName[which(panelInfoStr$OtherName == .x)], .cols = panelInfoStr$OtherName )

#stuff <- haps_2col %>% rename_at(vars(panelInfoStr$OtherName), ~ panelInfoStr$AmpliconName)


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
haps_2col_final <- haps_2col_final %>% select(!c(sdy_I183))
haps_2col_final <- arrange(haps_2col_final, Indiv) # sort by Indiv column

#View(haps_2col_final)

write_csv(haps_2col_final, file="haps_2col_final.csv")



## FILTERING TO DEAL WITH LATER
# Remove troublesome loci
loci2chuck <- c(
  "tag_id_1872",       #elevated error rate in WHOA!
  "tag_id_3194",       #elevated error rate in WHOA!
  "tag_id_1470",       #missing in ~half of individuals
  "OkiOts_120255-113"    #species ID
)

hap_fil2 <- filter_bad_loci(
  long_genos = hap_fil1,
  bad_loci = loci2chuck
)

# And  filter himissers 
hap_fil3 <- filter_missing_data(long_genos = hap_fil2, n_locs = 75)


#quit()