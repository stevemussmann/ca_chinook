import collections
import os
import pandas
import matplotlib.pyplot

class LocusDict():
	'Class for making dict to translate microhap genotypes to integer data'

	def __init__(self, df):
		self.df = df
		self.recodeAlleles = collections.defaultdict(dict)
		self.alleleCounts = collections.defaultdict(int)
	
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

	def getFreqs(self):
		alleleFreqs = collections.defaultdict(dict)

		colNames = list(self.df.columns)
		dupLoci = [item[:-2] for item in colNames]
		singleLoci = dupLoci[1::2]

		for locus in singleLoci:
			name1 = locus + "_1"
			name2 = locus + "_2"

			# get frequencies - maybe move to new function so not run multiple times
			freqs = dict(pandas.concat([self.df[name1], self.df[name2]]).dropna().value_counts())
			#print(freqs)

			alleleFreqs[locus] = freqs

		return alleleFreqs

	def countAlleles(self):
		for key in self.recodeAlleles.keys():
			self.alleleCounts[key] = len(self.recodeAlleles[key])
			#print(nAlleles)
		
		# use Counter to count occurrences of each value
		valCounts = collections.Counter(self.alleleCounts.values())
		#print(valCounts)	

		# make histogram of alleles per locus
		matplotlib.pyplot.bar(list(valCounts.keys()), list(valCounts.values())) # histogram
		matplotlib.pyplot.title("Allele Counts per Locus")
		matplotlib.pyplot.xlabel("Alleles per Locus")
		matplotlib.pyplot.ylabel("Observation Count")
		matplotlib.pyplot.savefig("histo.png")
