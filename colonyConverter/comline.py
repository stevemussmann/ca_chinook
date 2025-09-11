import argparse
import os.path
#import distutils.util

class ComLine():
	'Class for implementing command line options'

	def __init__(self, args):
		parser = argparse.ArgumentParser()
		parser._action_groups.pop()
		required = parser.add_argument_group('required arguments')
		optional = parser.add_argument_group('optional arguments')
		conversion = parser.add_argument_group('conversion arguments')

		required.add_argument("-f", "--infile",
							dest='infile',
							required=True,
							help="Specify input file (.csv format)."
		)
		optional.add_argument("-o", "--outfile",
							dest='outfile',
							default="default.txt",
							help="Specify output file name (default=default.txt)."
		)
		optional.add_argument("-l", "--pmissloc",
							dest='pmissloc',
							type=float,
							default=0.3,
							help="Enter the maximum allowable proportion of missing data for a locus (default = 0.3)."
		)
		optional.add_argument("-i", "--pmissind",
							dest='pmissind',
							type=float,
							default=0.3,
							help="Enter the maximum allowable proportion of missing data for an individual (default = 0.3)."
		)
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

		#check if at least one conversion option was used.
		if not [x for x in (self.args.colony, self.args.csv) if x is True]:
			print("")
			print("No format conversion options were selected.")
			print("You must choose at least one file format for output.")
			print("")
			raise SystemExit

		# check if files exist
		self.exists( self.args.infile )

	def exists(self, filename):
		if( os.path.isfile(filename) != True ):
			print("")
			print(filename, "does not exist")
			print("Exiting program...")
			print("")
			raise SystemExit
