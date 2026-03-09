import os
import pandas
import sys

from duplicates import Duplicates
from locusdict import LocusDict

class Microhap():
	'Class for operating on microhap genotype files'

	def __init__(self, infile, pmissLoc, pmissInd, mono, dup, t, k):
		self.mhFile = infile #input file name
		self.df = pandas.DataFrame()
		self.pmissLoc = pmissLoc # allowable proportion of missing data locus
		self.pmissInd = pmissInd # allowable proportion of missing data individual
		self.colonyData = pandas.DataFrame()
		self.mono = mono # boolean to control monomorphic locus filter
		self.dup = dup # boolean to control duplicate identification
		self.dupThresh = t # threshold for identifying duplicate individuals
		self.keepDups = k # method for keeping duplicates

		# deal with input file name to create log file name
		fn, ext = os.path.splitext(infile)
		self.log = fn + ".log"
		# remove old log file if still exists
		if os.path.isfile(self.log):
			os.remove(self.log)
		# write command used to launch program
		with open(self.log, 'a') as fh:
			fh.write("#microhapConvert.py was launched with command:\n#")
			comm = ' '.join(sys.argv)
			fh.write(comm)
			fh.write("\n\n")
		
		print("Reading input .csv file.")
		print("")
		self.df = pandas.read_csv(self.mhFile, index_col=0, header=0)


	def runFilters(self):
		# filter individuals
		self.filterInds()
		
		# filter loci
		self.filterLoci()

		# find duplicates
		if self.dup:
			dups = Duplicates(self.df, self.dupThresh, self.keepDups, self.log)
			dups.findDups()
			removeList = dups.removeDups() # get list of individuals to remove
			if removeList:
				self.dropRows(removeList) # drop duplicate individuals

		# remove monomorphic loci
		if self.mono == True:
			self.filterMono()


	def getLog(self):
		return self.log


	def getCounts(self):
		indsPerPop = self.df['Population ID'].value_counts().to_dict()
		return indsPerPop


	def getFinalCounts(self, pops):
		column_labels = ["Population ID"]
		popDF = pandas.DataFrame.from_dict(pops, orient='index', columns=column_labels)
		mergedDF = pandas.merge(self.df, popDF, how='left', left_index=True, right_index=True)
		indsPerPop = mergedDF['Population ID'].value_counts().to_dict()
		return indsPerPop


	def getDict(self):
		ld = LocusDict(self.df)
		ldict = ld.getUnique()
		#ld.countAlleles()
		return ldict


	def getFreqs(self):
		ld = LocusDict(self.df)
		freqs = ld.getFreqs()
		return freqs


	def parseFile(self, colonyBool, ckmrBool):
		# remove unneeded columns
		toRemove = ["sdy_sex", "hapstr", "rosa_pheno", "canonical_rosa_pheno", "percMicroHap"] # summary columns inserted by genotyping pipeline
		for col in toRemove:
			if col in self.df.columns:
				print("Removing", col, "column from input.")
				self.df.pop(col) # remove column

		## extract colony2 column if it exists
		if 'colony2' in self.df.columns:
			self.colonyData = self.df.pop('colony2')
		else:
			## exit with error if it doesn't exist
			# only do this if colony2 conversion requested
			if colonyBool == True or ckmrBool == True:
				print("\nERROR. The colony2 column is missing from your input file.")
				print("The 'colony2' column is needed when creating COLONY or CKMR outputs.")
				print("See documentation at https://github.com/stevemussmann/ca_chinook/tree/main/fileConverter")
				print("Exiting program...\n")
				raise SystemExit(1)

		# validate that remaining columns all end in _1 or _2
		columnNames = list(self.df.columns)
		#print(columnNames)
		for columnName in columnNames:
			if not columnName.endswith(("_1", "_2")):
				print("\nERROR.")
				print("Make sure all locus columns in your input file end in _1 or _2.\n")
				raise SystemExit(1)


		return self.colonyData


	def getPops(self):
		print("Extracting Population IDs...\n")
		try:
			pops = self.df.pop('Population ID').to_dict()
		except KeyError as e:
			print("\nERROR. The following column is missing from your input file:", e)
			print("Exiting program...\n")
			raise SystemExit(1)

		return pops


	def removeSnppit(self):
		print("Checking for presence of optional SNPPIT columns.")
		#list of all possible optional snppit columns
		optionalCols = ['POPCOLUMN_SEX', 'POPCOLUMN_REPRO_YEARS', 'POPCOLUMN_SPAWN_GROUP', 'OFFSPRINGCOLUMN_BORN_YEAR', 'OFFSPRINGCOLUMN_SAMPLE_YEAR', 'OFFSPRINGCOLUMN_AGE_AT_SAMPLING']

		remove = list() #will hold list of snppit columns that appear in pandas df
		snppitCols = pandas.DataFrame() #declare empty dataframe to be returned even if no optional columns were used.

		for col in optionalCols:
			if col in self.df.columns:
				remove.append(col) #add existing cols to remove list

		if remove:
			print("The following optional SNPPIT columns were detected in the input file:")
			for col in remove:
				print(col)
			print("")
			snppitCols = self.removeColumns(self.df, remove)
		else:
			print("No optional SNPPIT columns detected in input file.")
			print("")

		return snppitCols


	def removeColumns(self, df, removelist):
		junk = pandas.concat([df.pop(x) for x in removelist], axis=1)
		return junk

	
	def dropRows(self, removelist):
		self.df.drop(removelist, axis=0, inplace=True) # remove individuals

	
	def removeInds(self, blacklist):
		# remove blacklisted individuals
		print("\nRemoving blacklisted individuals (if present):")
		with open(self.log, 'a') as fh:
			fh.write("\nRemoving blacklisted individuals (if present):\n")
		removeInds = list()
		with open(blacklist, 'r') as fh:
			for line in fh:
				line = line.strip() # remove endline character
				print(line)
				removeInds.append(line)
				with open(self.log, 'a') as lfh:
					lfh.write(line)
					lfh.write("\n")
		if removeInds:
			self.dropRows(removeInds)


	def removeLoci(self, blacklist):
		# remove blacklisted columns
		print("\nRemoving blacklisted loci (if present):")
		with open(self.log, 'a') as fh:
			fh.write("\nRemoving blacklisted loci (if present):\n")
		removeLoci = list()
		with open(blacklist, 'r') as fh:
			for line in fh:
				removeLoci.append(line.strip())
		for col in removeLoci:
			a1 = col + "_1"
			a2 = col + "_2"
			if a1 in self.df.columns:
				print(a1) # print to stdout
				# write to log
				with open(self.log, 'a') as fh:
					fh.write(a1)
					fh.write("\n")
				self.df.pop(a1)
			if a2 in self.df.columns:
				print(a2) # print to stdout
				# write to log
				with open(self.log, 'a') as fh:
					fh.write(a1)
					fh.write("\n")
				self.df.pop(a2)


	def filterMono(self):
		removeLoci = list()
		findMono = LocusDict(self.df)
		mono = findMono.getUnique()
		findMono.countAlleles()
		#print(mono)
		#with open("mono.json", 'w') as jfile:
		#	json.dump(mono, jfile, indent='\t')
		for key, val in mono.items():
			if len(mono[key]) == 1:
				removeLoci.append(key)
		print("\nRemoving monomorphic loci:")
		with open(self.log, 'a') as fh:
			fh.write("\nRemoving monomorphic loci:\n")
		if removeLoci:
			for col in removeLoci:
				a1 = col + "_1"
				a2 = col + "_2"
				if a1 in self.df.columns:
					print(a1) # write to stdout
					# write to log
					with open(self.log, 'a') as fh:
						fh.write(a1)
						fh.write("\n")
					self.df.pop(a1)
				if a2 in self.df.columns:
					print(a2) # write to stdout
					# write to log
					with open(self.log, 'a') as fh:
						fh.write(a1)
						fh.write("\n")
					self.df.pop(a2)
			with open(self.log, 'a') as fh:
				fh.write("\n")
		else:
			print("None found!")
			with open(self.log, 'a') as fh:
				fh.write("None found!\n")


	def calcMissingLocPCT(self):
		# get number of individuals
		nInds = len(self.df)
		
		missLoc = self.df.isnull().sum(axis=0) # count missing genotypes per column (locus)
		missLocPCT = missLoc / nInds # divide missing data series by number of individuals

		return missLocPCT


	def filterLoci(self):
		# get number of loci
		if len(self.df.columns)%2 == 0: # test if even number of columns
			nLoci = int(len(self.df.columns)/2) # calculate number of loci; ensure output is integer
			#print(nLoci)
		else:
			print("\nERROR: Uneven number of locus columns in input file.")
			print("Exiting program...\n")
			raise SystemExit

		missLocPCT = self.calcMissingLocPCT()
		
		## uncomment for debugging
		#pandas.set_option('display.max_rows', None)
		#print(missLoc)
		
		removeLocPCT = missLocPCT[missLocPCT > self.pmissLoc].index # get list of column indexes to remove

		# print records that didn't pass missing data filter
		missRecords = missLocPCT.loc[missLocPCT.index.intersection(removeLocPCT)]

		print("\nMissing data proportion per removed locus:") # write to stdout
		# write to log
		with open(self.log, 'a') as fh:
			fh.write("\nMissing data proportion per removed locus:\n")
		#pandas.set_option("display.max_rows", None, "display.max_columns", None)
		print(missRecords.to_string(index=True)) # write to stdout
		# write to log
		with open(self.log, 'a') as fh:
			fh.write(missRecords.to_string(index=True))
			fh.write("\n")

		self.df.drop(removeLocPCT, axis=1, inplace=True)

	
	def calcMissingIndPCT(self):
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
		
		return missIndPCT
		
	
	def filterInds(self):
		missIndPCT = self.calcMissingIndPCT()

		removeIndPCT = missIndPCT[missIndPCT > self.pmissInd].index.to_list() # get list of row indexes to remove
		# print records that didn't pass missing data filter
		missRecords = missIndPCT.loc[missIndPCT.index.intersection(removeIndPCT)]

		print("\nMissing data proportion per individual (indiv):") # write to stdout
		# write to log
		with open(self.log, 'a') as fh:
			fh.write("\nMissing data proportion per individual (indiv):\n")
		#pandas.set_option("display.max_rows", None, "display.max_columns", None)
		print(missRecords.to_string(index=True))
		# write to log
		with open(self.log, 'a') as fh:
			fh.write(missRecords.to_string(index=True))
			fh.write("\n")
		
		self.dropRows(removeIndPCT)
