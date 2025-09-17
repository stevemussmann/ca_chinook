import argparse
import os.path
#import distutils.util

class ComLine():
	'Class for implementing command line options'

	def __init__(self, args):
		parser = argparse.ArgumentParser()
		parser._action_groups.pop()
		required = parser.add_argument_group('required arguments')
		filtering = parser.add_argument_group('filtering arguments')
		colony = parser.add_argument_group('colony arguments')
		conversion = parser.add_argument_group('conversion arguments')

		required.add_argument("-f", "--infile",
							dest='infile',
							required=True,
							help="Specify input file in .csv format (required)."
		)
		required.add_argument("-r", "--runname",
							dest='runname',
							required=True,
							help="Provide a unique name for the Colony run. Output files from Colony will receive this name (required)."
		)
		filtering.add_argument("-i", "--pmissind",
							dest='pmissind',
							type=float,
							default=0.3,
							help="Enter the maximum allowable proportion of missing data for an individual (default = 0.3)."
		)
		filtering.add_argument("-l", "--pmissloc",
							dest='pmissloc',
							type=float,
							default=0.3,
							help="Enter the maximum allowable proportion of missing data for a locus (default = 0.3)."
		)
		filtering.add_argument("-R", "--removeloci",
							dest='removeloci',
							help='Specify a list of loci to remove.'
		)
		filtering.add_argument("-m", "--mono",
							dest='mono',
							action='store_false',
							help="Remove monomorphic loci from final output (default = True)."
		)
		colony.add_argument("-d", "--droperr",
							dest='droperr',
							type=float,
							default=0.0005,
							help="Enter the assumed allelic dropout rate (default = 0.0005)."
		)
		colony.add_argument("-g", "--genoerr",
							dest='genoerr',
							type=float,
							default=0.0005,
							help="Enter the assumed genotyping error rate (default = 0.0005)."
		)
		colony.add_argument("-I", "--inbreed",
							dest='inbreed',
							type=int,
							default=0,
							choices={0,1},
							help="0 = inbreeding absent; 1 = inbreeding present (default = 0)."
		)
		colony.add_argument("-L", "--runlength",
							dest='runlength',
							type=int,
							default=2,
							choices={1,2,3,4},
							help="1/2/3/4 = Short/Medium/Long/VeryLong run (default = 2)."
		)
		colony.add_argument("-M", "--pmale",
							dest='pmale',
							type=float,
							default=0.5,
							help="Enter the assumed probability of father being among candidate parents (default = 0.5). Value is ignored if no candidate fathers provided in the dataset."
		)
		colony.add_argument("-F", "--pfemale",
							dest='pfemale',
							type=float,
							default=0.5,
							help="Enter the assumed probability of mother being among candidate parents (default = 0.5). Value is ignored if no candidate mothers provided in the dataset."
		)
		#optional.add_argument("-o", "--outfile",
		#					dest='outfile',
		#					default="default.txt",
		#					help="Specify output file name (default=default.txt)."
		#)
		conversion.add_argument("-c", "--csv",
							dest='csv',
							action='store_true',
							help="Write filtered csv format file."
		)
		conversion.add_argument("-C", "--colony",
							dest='colony',
							action='store_true',
							help="Write colony format file."
		)
		self.args = parser.parse_args()

		# check if certain values between 0.0 and 1.0
		if self.zeroOne(self.args.pmale) is False:
			print("ERROR: option -p must be between 0.0 and 1.0")
			raise SystemExit(1)
		
		if self.zeroOne(self.args.pfemale) is False:
			print("ERROR: option -P must be between 0.0 and 1.0")
			raise SystemExit(1)

		if self.zeroOne(self.args.droperr) is False:
			print("ERROR: option -d must be between 0.0 and 1.0")
			raise SystemExit(1)

		if self.zeroOne(self.args.genoerr) is False:
			print("ERROR: option -g must be between 0.0 and 1.0")
			raise SystemExit(1)

		if self.zeroOne(self.args.pmissind) is False:
			print("ERROR: option -i must be between 0.0 and 1.0")
			raise SystemExit(1)

		if self.zeroOne(self.args.pmissloc) is False:
			print("ERROR: option -l must be between 0.0 and 1.0")
			raise SystemExit(1)

		#check if at least one conversion option was used.
		if not [x for x in (self.args.colony, self.args.csv) if x is True]:
			print("")
			print("No format conversion options were selected.")
			print("You must choose at least one file format for output.")
			print("")
			raise SystemExit(1)

		# check if files exist
		self.exists( str(self.args.infile) )
		if self.args.removeloci:
			self.exists( str(self.args.removeloci) )

	def zeroOne(self, num):
		return 0.0 <= num <= 1.0

	def exists(self, filename):
		if( os.path.isfile(filename) != True ):
			print("")
			print(filename, "does not exist")
			print("Exiting program...")
			print("")
			raise SystemExit
