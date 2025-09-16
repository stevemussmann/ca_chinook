#!/usr/bin/env python3

from comline import ComLine
from mhconvert import MHconvert
from microhap import Microhap

import argparse
import json
import os
import sys

def main():
	input = ComLine(sys.argv[1:])

	# make directory for output if doesn't exist
	convertedDir = "convertedFiles"
	if os.path.exists(convertedDir) == False:
		os.mkdir(convertedDir)

	# make list of file formats; grab relevant options from argparse object
	d = vars(input.args)
	convDict = dict()
	convList = ['colony', 'csv']
	for key, value in d.items():
		if key in convList:
			convDict[key] = value
	#print(convDict)

	mhFile = Microhap(input.args.infile, input.args.pmissloc, input.args.pmissind) #initialize new file
	colonyData = mhFile.parseFile(input.args.colony) # parse file
	locusdict = mhFile.getDict() # make locus dictionary

	if input.args.removeloci:
		mhFile.removeLoci(input.args.removeloci) # remove blacklisted loci (if invoked)

	mhFile.runFilters() # run missing data filters
	
	# dump locus dictionary to text file (locusDictionary.json)
	jsonpath = os.path.join(convertedDir, "locusDictionary.json")
	with open(jsonpath, 'w') as jsonfile:
		json.dump(locusdict, jsonfile, indent='\t')

	# conversion process
	conversion = MHconvert(mhFile.df, input.args.infile, locusdict, colonyData, input.args.droperr, input.args.genoerr, input.args.pmale, input.args.pfemale, input.args.runname, input.args.inbreed, input.args.runlength, convertedDir)
	conversion.convert(convDict)

main()

raise SystemExit
