from colony import Colony
from csvf import CSVfiltered

import collections
import json
import os
import pandas

class MHconvert():
	'Class for converting pandas dataframes into various genotype files'

	def __init__(self, df, infile, ldict, cDat, derr, gerr, pm, pf, runname, inbreed, runlen, cdir, afreqs):
		self.df = df
		self.ldict = ldict
		self.infile = infile
		self.cDat = cDat # colony data; offspring, male parent, female parent, etc.
		self.derr = derr # allelic dropout rate
		self.gerr = gerr # genotyping error rate
		self.pmale = pm # probability of father being present among candidates
		self.pfemale = pf # probability of mother being present among candidates
		self.runname = runname
		self.inbreed = inbreed
		self.runlen = runlen
		self.suffix = {'colony': 'Dat', 'csv': 'csv'}
		self.convertedDir = cdir # directory to hold converted files
		self.snpdf = pandas.DataFrame() #dataframe to hold SNPs if SNP output option is used
		self.alleleFreqs = afreqs
		
		self.convSNP(self.df)

	def convert(self, d):
		output = list()
		for filetype, boolean in d.items():
			if boolean == True:
				print("Converting to", filetype, "format file.")
				output = self.convert_to(filetype)
				self.printOutput(output, self.infile, self.suffix[filetype])

	def conv_csv(self):
		#print("This function will print a filtered .csv file")
		csv = CSVfiltered(self.df)
		output = csv.convert()
		return output

	def conv_colony(self): 
		#print("This function will convert to colony format.")
		cy = Colony(self.df, self.ldict, self.cDat, self.derr, self.gerr, self.pmale, self.pfemale, self.runname, self.inbreed, self.runlen)
		output = cy.convert()
		return output
	
	def convert_to(self, name: str):
		conv = f"conv_{name}"
		output = list()
		if hasattr(self, conv) and callable(func := getattr(self, conv)):
			output = func()
		else:
			print("Function not found for converting", name, "format.")
			print("Exiting program...")
			print("")
			raise SystemExit(1)
		return output

	def printOutput(self, output, fileName, suffix):
		# if colony conversion, use Colony2.Dat as output name
		if suffix == "Dat":
			outName = os.path.join(self.convertedDir, "Colony2.Dat")
		else:
			# make new file name for writing
			fileName = fileName.replace(" ", "_") #replace spaces in original filename if they exist
			nameList = fileName.split('.')
			nameList.pop() #remove old extension
			nameList.append(suffix) #add new file extension
			outName = '.'.join(nameList)
			outName = os.path.join(self.convertedDir, outName)

		print("Writing to", outName)
		print("")

		fh = open(outName, 'w')

		for line in output:
			fh.write(line)
			fh.write("\n")

	def convSNP(self, df):
		superDict = dict() # will hold nested dicts created for each locus
		colNames = list(self.df.columns)
		dupLoci = [item[:-2] for item in colNames]
		singleLoci = dupLoci[1::2]

		for locus in singleLoci:
			locDict = collections.defaultdict(lambda: collections.defaultdict(int)) # nested dict of dicts with default int value so that I can add value without first needing to check if key already exists

			for key, val in self.alleleFreqs[locus].items():
				hapList = list(key)
				count=0
				for char in hapList:
					locDict[count][char] += int(val)
					count += 1

			superDict[locus] = locDict # superDict holds all individual locDicts

		## uncomment next 3 lines to verify nested dict format was created as desired
		#jsonpath = os.path.join(os.getcwd(), "snpFreqsDict.json")
		#with open(jsonpath, 'w') as jsonfile:
		#	json.dump(superDict, jsonfile, indent='\t')

		removeList = list() # hold loci that will be removed from dataframe because there are no biallelic SNPs
		keepDict = dict()
		for locus, d in superDict.items():
			potentialKeep = list()
			for position, d2 in d.items():
				# test if position is biallelic
				if len(d2) == 2:
					potentialKeep.append(position) # record biallelic positions

			# find loci with no biallelic positions and push to removeList
			if len(potentialKeep) == 0:
				removeList.append(locus)
			# if single biallelic position, push locus and position to keepDict
			elif len(potentialKeep) == 1:
				keepDict[locus] = potentialKeep[0]
			# if multiple biallelic positions, find position with greatest minor allele count
			else:
				mac = dict() # dict to hold all minor allele counts for locus
				for position in potentialKeep:
					# determine minor allele at each position
					minAllele = min(superDict[locus][position], key=superDict[locus][position].get)

					# record minor allele count for position in mac dict
					mac[position] = superDict[locus][position][minAllele]
				
				# get position that has minor allele with greatest minor allele count
				maxminAllele = max(mac, key=mac.get)
				keepDict[locus] = maxminAllele # add retained position to keepDict
				#print(mac)

		#print(removeList)
			
		#print(keepDict)

