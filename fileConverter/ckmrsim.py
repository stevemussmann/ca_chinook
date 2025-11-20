import os
import pandas

class CKMR():
	'Class for outputting filtered pandas dataframe to CKMRsim (csv) format'

	def __init__(self, df, cDat):
		self.df = df
		self.cDat = cDat # colony data (potential male parent, female parent, offspring)
		#self.ldict = ldict

	def convert(self):
		offspring = list() # list to contain offspring genotypes
		parents = list() # list to contain parent genotypes

		# handle header line
		headerList = list()
		headerList.append("indiv")
		colNames = list(self.df.columns)
		headerList.extend(colNames)
		headerStr = "\t".join(headerList)

		# add header line to each
		offspring.append(headerStr)
		parents.append(headerStr)

		# handle offspring genotypes
		for (sampleName, row) in self.df.iterrows():
			if self.cDat[sampleName].casefold() == "offspring".casefold():
				sampleStr = self.parseInd(sampleName, row)
				offspring.append(sampleStr)
			else:
				sampleStr = self.parseInd(sampleName, row)
				parents.append(sampleStr)
		
		# combine offspring and parents lists into output list
		output = offspring + parents 

		return output

	def parseInd(self, sampleName, row):
		sampleList = list()
		locusList = list()
		sampleList.append(str(sampleName))

		for (locus, genotype) in row.items():
			if not pandas.isnull(genotype):
				locusList.append(str(genotype))
			else:
				locusList.append("NA")

		locusStr = "\t".join(locusList)
		sampleList.append(locusStr)
		sampleStr = "\t".join(sampleList)

		return sampleStr
