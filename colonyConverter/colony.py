import os
import pandas
import random

class Colony():
	'Class for converting pandas dataframe to colony format'

	def __init__(self, df, ldict, cDat, derr, gerr, pm, pf, runname):
		self.df = df
		self.ldict = ldict
		self.cDat = cDat # colony data (potential male parent, female parent, offspring)
		self.derr = derr # allelic dropout rate
		self.gerr = gerr # genotyping error rate
		self.pmale = pm # probability of father being present among candidates
		self.pfemale = pf # probability of mother being present among candidates
		self.runname = runname
		#self.convertedDir = convDir

	def convert(self):
		output = list()

		randseed = random.randint(1000, 9999) # 4-digit random number seed
		#offspring = len(self.df) # number of individuals in dataframe
		colonyCounts = self.cDat.str.lower().value_counts().to_dict() # counts of offspring and parents
		#print(colonyCounts)

		loci = nLoci = int(len(self.df.columns)/2) # number of loci in dataframe

		datasetnamelist = list()
		datasetnamelist.append("'")
		datasetnamelist.append(str(self.runname))
		datasetnamelist.append("'")
		datasetline = "".join(datasetnamelist)
		output.append(datasetline)
		output.append(datasetline)

		offspringline = str(colonyCounts["offspring"]) + "      ! Number of offspring in the sample"
		output.append(offspringline)

		lociline = str(loci) + "       ! Number of loci"
		output.append(lociline)

		randseedline = str(randseed) + "      ! Seed for random number generator"
		output.append(randseedline)

		output.append("0         ! 0/1=Not updating/updating allele frequency")
		output.append("2         ! 2/1=Dioecious/Monoecious species")
		output.append("0         ! 0/1=Inbreeding absent/present")
		output.append("0         ! 0/1=Diploid species/HaploDiploid species")
		output.append("0  0      ! 0/1=Polygamy/Monogamy for males & females")
		output.append("0         ! 0/1 = Clone inference = No/Yes")
		output.append("1         ! 0/1=Scale full sibship=No/Yes")
		output.append("0         ! 0/1/2/3/4=No/Weak/Medium/Strong sibship prior; 4=Optimal sibship prior for Ne")
		output.append("0         ! 0/1=Unknown/Known population allele frequency")
		output.append("1         ! Number of runs")
		output.append("2         ! 1/2/3/4 = Short/Medium/Long/VeryLong run")
		output.append("1         ! 0/1=Monitor method by Iterate#/Time in second")
		output.append("1         ! Monitor interval in Iterate# / in seconds")
		output.append("0         ! 0/1=DOS/Windows version")
		output.append("1         ! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method")
		output.append("2         ! 0/1/2/3=Low/Medium/High/VeryHigh precision")
		output.append("")

		# print string of locus names
		locusNames = self.getLocusNames()
		locusString = " ".join(locusNames)
		output.append(locusString)

		# print string of marker types. 0 = codominant, 1 = dominant
		mtString = self.prepValues(loci, 0)
		output.append(mtString)

		# print string of allelic dropout rates
		adrString = self.prepValues(loci, self.derr)
		output.append(adrString)

		# print string of genotyping error rates
		gerString = self.prepValues(loci, self.gerr)
		output.append(gerString)

		
		for (sampleName, row) in self.df.iterrows():
			if self.cDat[sampleName].casefold() == "offspring".casefold():
				# if popmap value = offspring
				sampleList = list()
				locusList = list()
				sampleList.append(str(sampleName))
				for (locus, genotype) in row.items():
					loc = locus[:-2]
					if not pandas.isnull(genotype):
						locusList.append(str(self.ldict[loc][genotype]))
					else:
						locusList.append("0")
				locusStr = " ".join(locusList)
				sampleList.append(locusStr)
				sampleStr = " ".join(sampleList)
				output.append(sampleStr)

		output.append("")

		# set probabilities that male and/or female parent included among candidates
		if ("male" not in colonyCounts) and ("female" not in colonyCounts):
			output.append("0.0  0.0     !prob. of dad/mum included in the candidates")
		else:
			templine = list() # list to build probabilities line
			if "male" in colonyCounts:
				templine.append(str(self.pmale))
			else:
				templine.append("0.0")

			if "female" in colonyCounts:
				templine.append(str(self.pfemale))
			else:
				templine.append("0.0")

			templine.append("      !prob. of dad/mum included in the candidates")
			mfCountString = " ".join(templine)

			output.append(mfCountString)
	
		# set number of male and/or female parents
		if ("male" not in colonyCounts) and ("female" not in colonyCounts):
			output.append("0  0         !numbers of candidate males & females")
		else:
			templine = list() # list to build line from
			if "male" in colonyCounts:
				templine.append(str(colonyCounts["male"]))
			else:
				templine.append("0")

			if "female" in colonyCounts:
				templine.append(str(colonyCounts["female"]))
			else:
				templine.append("0")

			templine.append("      !numbers of candidate males & females")
			mfCountString = " ".join(templine)

			output.append(mfCountString)

		output.append("")
		# genotypes of male parent candidates go here
		if "male" in colonyCounts:
			for (sampleName, row) in self.df.iterrows():
				if self.cDat[sampleName].casefold() == "male".casefold():
					sampleList = list()
					locusList = list()
					sampleList.append(str(sampleName))
					for (locus, genotype) in row.items():
						loc = locus[:-2]
						if not pandas.isnull(genotype):
							locusList.append(str(self.ldict[loc][genotype]))
						else:
							locusList.append("0")
					locusStr = " ".join(locusList)
					sampleList.append(locusStr)
					sampleStr = " ".join(sampleList)
					output.append(sampleStr)

		output.append("")
		# genotypes of female parent candidates go here
		if "female" in colonyCounts:
			for (sampleName, row) in self.df.iterrows():
				if self.cDat[sampleName].casefold() == "female".casefold():
					sampleList = list()
					locusList = list()
					sampleList.append(str(sampleName))
					for (locus, genotype) in row.items():
						loc = locus[:-2]
						if not pandas.isnull(genotype):
							locusList.append(str(self.ldict[loc][genotype]))
						else:
							locusList.append("0")
					locusStr = " ".join(locusList)
					sampleList.append(locusStr)
					sampleStr = " ".join(sampleList)
					output.append(sampleStr)

		output.append("")
		output.append("0  0        !#known father-offspring dyads, paternity exclusion threshold")
		output.append("")
		output.append("0  0        !#known mother-offspring dyads, maternity exclusion threshold")
		output.append("")
		output.append("0           !#known paternal sibship with unknown fathers")
		output.append("")
		output.append("0           !#known maternal sibship with unknown mothers")
		output.append("")
		output.append("0           !#known paternity exclusions")
		output.append("")
		output.append("0           !#known maternity exclusions")
		output.append("")
		output.append("0           !#known paternal sibship exclusions")
		output.append("")
		output.append("0           !#known maternal sibship exclusions")
		output.append("")

		return output

	def getLocusNames(self):
		colNames = list(self.df.columns) #get column names from pandas dataframe
		dupLoci = [item[:-2] for item in colNames] #strip allele identifiers from ends
		#print(dupLoci)

		singleLoci = dupLoci[1::2] #keep odd numbered list elements
		#print(singleLoci)

		return singleLoci

	def prepValues(self, nloci, val):
		valList = [str(val)] * nloci
		valString = " ".join(valList)
		return valString
		
