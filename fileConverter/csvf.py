import os
import pandas

class CSVfiltered():
	'Class for outputting filtered pandas dataframe to csv format'

	def __init__(self, df):
		self.df = df
		#self.ldict = ldict
		#self.cDat = cDat # colony data (potential male parent, female parent, offspring)

	def convert(self):
		output = list()
		# handle header line
		headerList = list()
		headerList.append("indiv")
		colNames = list(self.df.columns)
		headerList.extend(colNames)
		headerStr = ",".join(headerList)
		output.append(headerStr)

		# handle genotypes
		for (sampleName, row) in self.df.iterrows():
			sampleList = list()
			locusList = list()
			sampleList.append(str(sampleName))
			for (locus, genotype) in row.items():
				if not pandas.isnull(genotype):
					locusList.append(str(genotype))
				else:
					locusList.append("NA")
			locusStr = ",".join(locusList)
			sampleList.append(locusStr)
			sampleStr = ",".join(sampleList)
			output.append(sampleStr)

		return output
