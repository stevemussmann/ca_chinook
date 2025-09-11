import collections
import os
import pandas

class LocusDict():
	'Class for making dict to translate microhap genotypes to integer data'

	def __init__(self, df):
		self.df = df
		self.recodeAlleles = collections.defaultdict(dict)
	
	def getUnique(self):
		colNames = list(self.df.columns)
		dupLoci = [item[:-2] for item in colNames]
		singleLoci = dupLoci[1::2]

		# for each locus
		for locus in singleLoci:
			name1 = locus + "_1"
			name2 = locus + "_2"

			# find all alleles for the locus (both columns of alleles combined)
			alleles = list(pandas.concat([self.df[name1], self.df[name2]]).dropna().unique())

			# enumerate all alleles
			tempdict = dict(enumerate(alleles))

			# swap values and add 1 to original keys to make individual locus dict
			swap = {value: str(key+1) for key, value in tempdict.items()}

			self.recodeAlleles[locus] = swap
			#print(swap)

		return self.recodeAlleles
		

