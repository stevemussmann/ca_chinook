from colony import Colony
from csvf import CSVfiltered

import os

class MHconvert():
	'Class for converting pandas dataframes into various genotype files'

	def __init__(self, df, infile, ldict, cDat, derr, gerr, pm, pf, runname, inbreed, runlen, cdir):
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
		#if os.path.exists(self.convertedDir) == False:
		#	os.mkdir(self.convertedDir)

	def convert(self, d):
		output = list()
		for filetype, boolean in d.items():
			if boolean == True:
				print("Converting to", filetype, "format file.")
				output = self.convert_to(filetype)
				self.printOutput(output, self.infile, self.suffix[filetype])

	def conv_csv(self):
		print("This function will print a filtered .csv file")
		csv = CSVfiltered(self.df)
		output = csv.convert()
		return output

	def conv_colony(self): 
		print("This function will convert to colony format.")
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
