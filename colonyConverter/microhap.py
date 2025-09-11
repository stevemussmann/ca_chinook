import pandas
#import os.path
#import distutils.util

from locusdict import LocusDict

class Microhap():
	'Class for operating on microhap genotype files'

	def __init__(self, infile, pmissLoc, pmissInd):
		self.mhFile = infile #input file name
		self.df = pandas.DataFrame()
		self.pmissLoc = pmissLoc # allowable proportion of missing data locus
		self.pmissInd = pmissInd # allowable proportion of missing data individual

	def getDict(self):
		ld = LocusDict(self.df)
		ldict = ld.getUnique()
		return ldict

	def parseFile(self):
		print("Reading input .csv file.")
		print("")
		self.df = pandas.read_csv(self.mhFile, index_col=0, header=0)

		# filter loci
		self.filterLoci()

		# filter individuals
		self.filterInds()

	def filterLoci(self):
		# get number of loci
		if len(self.df.columns)%2 == 0: # test if even number of columns
			nLoci = int(len(self.df.columns)/2) # calculate number of loci; ensure output is integer
		else:
			print("\nERROR: Uneven number of locus columns in input file.")
			print("Exiting program...\n")
			raise SystemExit

		missLoc = self.df.isnull().sum(axis=0) # count missing genotypes per column (locus)
		missLocPCT = missLoc / nLoci # divide missing data series by number of loci
		removeLocPCT = missLocPCT[missLocPCT > self.pmissLoc].index # get list of column indexes to remove

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
		
		self.df.drop(removeIndPCT, axis=0, inplace=True)
		
