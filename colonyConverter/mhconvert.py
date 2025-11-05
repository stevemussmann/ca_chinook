from colony import Colony
from csvf import CSVfiltered
from sequoia import Sequoia
from snppit import Snppit

import collections
import json
import os
import pandas
import warnings

# PerformanceWarning STFU
warnings.simplefilter(action='ignore', category=pandas.errors.PerformanceWarning)

class MHconvert():
	'Class for converting pandas dataframes into various genotype files'

	def __init__(self, df, infile, ldict, cDat, derr, gerr, pm, pf, runname, inbreed, runlen, cdir, afreqs, snppitCols, pops, snppitmap):
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
		self.suffix = {'colony': 'Dat', 'csv': 'csv', 'sequoia': 'sequoia', 'snppit': 'snppit'}
		self.convertedDir = cdir # directory to hold converted files
		self.snpdf = pandas.DataFrame() #dataframe to hold SNPs if SNP output option is used
		self.alleleFreqs = afreqs
		self.snppitCols = snppitCols
		self.pops = pops # dict of population information for all individuals
		self.snppitmap = snppitmap
		
		# add if statement here to only do SNP conversion if SNP file formats requested
		kd = self.findSNP(self.df) # identify SNP positions to keep for SNP format output files
		self.snpDF = self.convSNP(kd, self.df) # make dataframe of SNPs
		self.snpDF.to_excel('output.xlsx', index=True)


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

	def conv_sequoia(self):
		#print("This function will convert to sequoia format.")
		seq = Sequoia(self.snpDF, self.convertedDir)
		output = seq.convert(self.snppitCols)
		return output
	
	def conv_snppit(self):
		#print("This function will convert to SNPPIT format.")
		snppit = Snppit(self.snpDF, self.pops)
		output = snppit.convert(self.snppitmap, self.snppitCols)
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

	def convSNP(self, kd, df):
		print("Making new SNP dataframe.\n")
		newDF = df.copy(deep=True) # make deep copy of microhap dataframe
		snpDF = pandas.DataFrame() # make dataframe that will hold SNP calls
		for locus, pos in kd.items():
			name1 = locus + "_1"
			name2 = locus + "_2"

			newDF[name1] = newDF[name1].str.slice(pos, pos+1) # extract SNP allele 1
			newDF[name2] = newDF[name2].str.slice(pos, pos+1) # extract SNP allele 2

			snpDF[locus] = newDF[name1] + newDF[name2] # concatenate SNP alleles 1 and 2
			snpDF[locus] = snpDF[locus].apply(lambda x: ''.join(sorted(str(x))) if pandas.notna(x) else "0") # sort new string

			# comment above line and uncomment next two lines for alternate handling of NA values
			#snpDF[locus] = snpDF[locus].apply(lambda x: ''.join(sorted(str(x)))) # sort new string
		#snpDF.replace('ann', pandas.NA, inplace=True) # make sure NA values are treated properly
		#print(snpDF)

		return snpDF

	def findSNP(self, df):
		print("SNP file format output invoked - Finding SNPs with highest MAC values per locus.\n")
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

		if removeList:
			print("The following loci were removed from SNP file outputs because they contained no biallelic SNPs:")
			for l in removeList:
				print(str(l))
			print("\n")

		return keepDict

