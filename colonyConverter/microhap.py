#import json
import pandas
#import os.path
#import distutils.util

from locusdict import LocusDict

class Microhap():
	'Class for operating on microhap genotype files'

	def __init__(self, infile, pmissLoc, pmissInd, mono):
		self.mhFile = infile #input file name
		self.df = pandas.DataFrame()
		self.pmissLoc = pmissLoc # allowable proportion of missing data locus
		self.pmissInd = pmissInd # allowable proportion of missing data individual
		self.colonyData = pandas.DataFrame()
		self.mono = mono # boolean to control monomorphic locus filter

	def getDict(self):
		ld = LocusDict(self.df)
		ldict = ld.getUnique()
		ld.countAlleles()
		return ldict

	def parseFile(self, colonyBool):
		print("Reading input .csv file.")
		print("")
		self.df = pandas.read_csv(self.mhFile, index_col=0, header=0)

		# remove unneeded columns
		toRemove = ["sdy_sex", "hapstr", "rosa_pheno", "percMicroHap"] # summary columns inserted by genotyping pipeline
		for col in toRemove:
			if col in self.df.columns:
				print("Removing", col, "column from input.")
				self.df.pop(col) # remove column

		## extract colony2 column; exit with error if it doesn't exist
		# only do this if colony2 conversion requested
		if colonyBool == True:
			try:
				self.colonyData = self.df.pop('colony2')
			except KeyError as e:
				print("\nERROR. The following column is missing from your input file:", e)
				print("The 'colony2' column should exist and contain information (status as potential male parent, female parent, or offspring) for all individuals.")
				print("Exiting program...\n")
				raise SystemExit(1)

		# validate that remaining columns all end in _1 or _2
		columnNames = list(self.df.columns)
		for columnName in columnNames:
			if not columnName.endswith(("_1", "_2")):
				print("\nERROR.")
				print("Make sure all locus columns in your input file end in _1 or _2.\n")
				raise SystemExit(1)


		return self.colonyData
		
	def removeLoci(self, blacklist):
		# remove blacklisted columns
		print("\nRemoving blacklisted loci (if present):")
		removeLoci = list()
		with open(blacklist, 'r') as fh:
			for line in fh:
				removeLoci.append(line.strip())
		for col in removeLoci:
			a1 = col + "_1"
			a2 = col + "_2"
			if a1 in self.df.columns:
				print(a1)
				self.df.pop(a1)
			if a2 in self.df.columns:
				print(a2)
				self.df.pop(a2)

	def runFilters(self):
		# filter individuals
		self.filterInds()
		
		# filter loci
		self.filterLoci()

		# remove monomorphic loci
		if self.mono == True:
			self.filterMono()
	
	def filterMono(self):
		removeLoci = list()
		findMono = LocusDict(self.df)
		mono = findMono.getUnique()
		#print(mono)
		#with open("mono.json", 'w') as jfile:
		#	json.dump(mono, jfile, indent='\t')
		for key, val in mono.items():
			if len(mono[key]) == 1:
				removeLoci.append(key)
		print("\nRemoving monomorphic loci:")
		if removeLoci:
			for col in removeLoci:
				a1 = col + "_1"
				a2 = col + "_2"
				if a1 in self.df.columns:
					print(a1)
					self.df.pop(a1)
				if a2 in self.df.columns:
					print(a2)
					self.df.pop(a2)
		else:
			print("None found!")

	def filterLoci(self):
		# get number of loci
		if len(self.df.columns)%2 == 0: # test if even number of columns
			nLoci = int(len(self.df.columns)/2) # calculate number of loci; ensure output is integer
			#print(nLoci)
		else:
			print("\nERROR: Uneven number of locus columns in input file.")
			print("Exiting program...\n")
			raise SystemExit

		# get number of individuals
		nInds = len(self.df)

		missLoc = self.df.isnull().sum(axis=0) # count missing genotypes per column (locus)
		missLocPCT = missLoc / nInds # divide missing data series by number of individuals
		
		## uncomment for debugging
		#pandas.set_option('display.max_rows', None)
		#print(missLoc)
		
		removeLocPCT = missLocPCT[missLocPCT > self.pmissLoc].index # get list of column indexes to remove
		#print(removeLocPCT)

		# print records that didn't pass missing data filter
		missRecords = missLocPCT.loc[missLocPCT.index.intersection(removeLocPCT)]

		print("\nMissing data proportion per removed locus:")
		#pandas.set_option("display.max_rows", None, "display.max_columns", None)
		print(missRecords.to_string(index=True))

		self.df.drop(removeLocPCT, axis=1, inplace=True)

	def filterInds(self):
		# get number of individuals
		nInds = len(self.df)
		
		# get number of loci
		if len(self.df.columns)%2 == 0: # test if even number of columns
			nAlleles = len(self.df.columns) # calculate number of alleles per ind
		else:
			print("\nERROR: Odd count of locus columns post-missing data (locus) filter.")
			print("Exiting program...\n")
			raise SystemExit

		missInd = self.df.isnull().sum(axis=1) # count missing genotypes per individual
		missIndPCT = missInd / nAlleles # divide missing data series by number of alleles
		removeIndPCT = missIndPCT[missIndPCT > self.pmissInd].index # get list of row indexes to remove
		
		# print records that didn't pass missing data filter
		missRecords = missIndPCT.loc[missIndPCT.index.intersection(removeIndPCT)]

		print("\nMissing data proportion per individual (indiv):")
		#pandas.set_option("display.max_rows", None, "display.max_columns", None)
		print(missRecords.to_string(index=True))
		
		self.df.drop(removeIndPCT, axis=0, inplace=True)
		
