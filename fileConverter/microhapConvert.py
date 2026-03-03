#!/usr/bin/env python3

from comline import ComLine
from mhconvert import MHconvert
from microhap import Microhap
from stats import GTStats

import argparse
import json
import os
import pandas
import sys

def main():
	input = ComLine(sys.argv[1:])

	# make directory for output if doesn't exist
	convertedDir = "convertedFiles"
	if os.path.exists(convertedDir) == False:
		os.mkdir(convertedDir)

	# make list of file formats; grab relevant options from argparse object
	d = vars(input.args)
	convDict = dict() # all file formats
	snpDict = dict() # only snp file formats
	convList = ['ckmr', 'colony', 'csv', 'genepop', 'sequoia', 'snppit']
	snpList = ['sequoia', 'snppit']
	for key, value in d.items():
		if key in convList:
			convDict[key] = value
		if key in snpList:
			snpDict[key] = value

	mhFile = Microhap(input.args.infile, input.args.pmissloc, input.args.pmissind, input.args.mono) #initialize new file
	logfile = mhFile.getLog() # retrieve logfile name

	startIndsPerPop = mhFile.getCounts() # get counts of individuals per population at beginning of analysis
	
	# pull out special columns
	snppitCols = mhFile.removeSnppit() #removes optional columns for SNPPIT
	pops = mhFile.getPops() #remove populations column; variable 'pops' is a dict

	# parse file
	colonyData = mhFile.parseFile(input.args.colony, input.args.ckmr)

	try:
		locusdict = mhFile.getDict() # make locus dictionary
		alleleFreqs = mhFile.getFreqs() # get number occurrences of each allele per locus
	except KeyError as e:
		print("\nKeyError:", e, "was not found.")
		print("Make sure all locus columns in your input file end in _1 or _2.\n")
		raise SystemExit(1)

	# calculate beginning stats before removing any loci/individuals
	mLocStart = mhFile.calcMissingLocPCT()
	mIndStart = mhFile.calcMissingIndPCT()

	mLocStartList = pandas.Series(mLocStart).to_list()
	mIndStartList = pandas.Series(mIndStart).to_list()

	mLocStartStats = GTStats(mLocStartList)
	mIndStartStats = GTStats(mIndStartList)

	mLocStartStats.calcStats()
	mIndStartStats.calcStats()

	mLocStartStats.printStats(logfile, "loci", "before", True) # boolean value should be true for loci; false for individuals
	mIndStartStats.printStats(logfile, "individuals", "before", False)

	# remove 
	if input.args.removeloci:
		mhFile.removeLoci(input.args.removeloci) # remove blacklisted loci (if invoked)
	mhFile.runFilters() # run missing data filters

	# calculate ending stats after running all filters
	mLocEnd = mhFile.calcMissingLocPCT()
	mIndEnd = mhFile.calcMissingIndPCT()

	mLocEndList = pandas.Series(mLocEnd).to_list()
	mIndEndList = pandas.Series(mIndEnd).to_list()

	mLocEndStats = GTStats(mLocEndList)
	mIndEndStats = GTStats(mIndEndList)

	mLocEndStats.calcStats()
	mIndEndStats.calcStats()

	mLocEndStats.printStats(logfile, "loci", "after", True) # boolean value should be true for loci; false for individuals
	mIndEndStats.printStats(logfile, "individuals", "after", False)

	# calculate stats after running filters
	endIndsPerPop = mhFile.getFinalCounts(pops) # get count of remaining individuals after filtering
	
	# dump locus dictionary to text file (locusDictionary.json)
	if input.args.colony == True or input.args.genepop == True:
		jsonpath = os.path.join(convertedDir, "locusDictionary.json")
		print("\nPrinting locus dictionary used for COLONY and genepop format conversions to", str(jsonpath), "\n")
		with open(jsonpath, 'w') as jsonfile:
			json.dump(locusdict, jsonfile, indent='\t')

	# conversion process
	conversion = MHconvert(mhFile.df, input.args.infile, locusdict, colonyData, input.args.droperr, input.args.genoerr, input.args.pmale, input.args.pfemale, input.args.runname, input.args.inbreed, input.args.runlength, convertedDir, alleleFreqs, snppitCols, pops, input.args.snppitmap, snpDict)
	conversion.convert(convDict)

	# print starting and ending individuals per population
	print("Printing number of individuals per population before/after filtering...")
	print("Population\tPre-Filter\tPost-Filter")
	with open(logfile, 'a') as fh:
		fh.write("\nIndividuals per population before/after filtering:\n")
		fh.write("Population\tPre-Filter\tPost-Filter\n")
	popKeys = list(startIndsPerPop.keys())
	popKeys.sort()
	for pop in popKeys:
		final = 0
		if pop in endIndsPerPop:
			final = endIndsPerPop[pop]
		print(f"{pop}\t{startIndsPerPop[pop]}\t{final}")
		with open(logfile, 'a') as fh:
			fh.write(str(pop))
			fh.write("\t")
			fh.write(str(startIndsPerPop[pop]))
			fh.write("\t")
			fh.write(str(final))
			fh.write("\n")
			

main()

raise SystemExit
