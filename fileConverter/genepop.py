from popmap import Popmap

import os
import pandas

class Genepop():
	'Class for converting pandas dataframe to Genepop format'

	def __init__(self, df, popmap, convDir, ldict):
		self.df = df
		self.pops = popmap
		self.convertedDir = convDir
		self.ldict = ldict

	def convert(self):
		pm = Popmap(self.pops)
		mapDict = pm.parseMap()
		mapDict = dict(sorted(mapDict.items())) # sort dict

		# open file for writing population map
		popmapOut = os.path.join(self.convertedDir, "genepopmap.txt")
		fh = open(popmapOut, 'w')

		lineList = list()

		lineList.append('Title line:""')

		self.df.columns = self.df.columns.str.slice(stop=-2)
		uniq = pandas.unique(pandas.Series(self.df.columns.tolist())).tolist() # get unique column names
		for columnName in uniq:
			lineList.append(columnName)

		for (pop, num) in mapDict.items():
			lineList.append("Pop")
			for sampleName, row in self.df.iterrows():
				sampleList = list()
				if self.pops[sampleName] == pop:
					# write to popmap
					fh.write(sampleName)
					fh.write("\t")
					fh.write(pop)
					fh.write("\n")

					# append data to sampleList
					sampleList.append(sampleName)
					sampleList.append(",")
					sampleList.append("")
					counter=0 # counter to determine when to add spaces
					locusList = list()
					for (locus, genotype) in row.items():
						if pandas.isna(genotype) == True:
							locusList.append("000")
						else:
							locusList.append(self.ldict[locus][genotype])

						counter+=1 # increment counter
						# on even numbers
						if counter % 2 == 0:
							locusStr = ''.join(locusList) # make string for locus
							sampleList.append(locusStr) # append string to list for sample
							locusList = list() # make new locustList for next locus
					sampleStr = ' '.join(sampleList)
					lineList.append(sampleStr)


		fh.close()

		return lineList
