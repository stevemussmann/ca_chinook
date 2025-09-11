#!/usr/bin/env python3

from comline import ComLine
from mhconvert import MHconvert
from microhap import Microhap

import argparse
#import os
import sys

def main():
	input = ComLine(sys.argv[1:])

	# make list of file formats; grab relevant options from argparse object
	d = vars(input.args)
	convDict = dict()
	convList = ['colony', 'csv']
	for key, value in d.items():
		if key in convList:
			convDict[key] = value
	#print(convDict)

	mhFile = Microhap(input.args.infile, input.args.pmissloc, input.args.pmissind)
	mhFile.parseFile()
	locusdict = mhFile.getDict()

	# conversion process
	conversion = MHconvert(mhFile.df, input.args.infile, locusdict)
	conversion.convert(convDict)

main()

raise SystemExit
