#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
library("optparse")

option_list = list(
	make_option(
		c("-f", "--file"), 
		type="character", 
		default=NULL, 
		help="Illumina sample sheet file name", 
		metavar="file"
	)
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
	print_help(opt_parser)
	stop("Input sample sheet must be provided.", call.=FALSE)
}

# before running this script for the first time, modify the path in the line below so that it points to the location of sample-sheet-processing-functions.R on your computer
source("~/local/src/mega-simple-microhap-snakeflow/preprocess/sample-sheet-processing-functions.R")

tryCatch({
		create_samples_and_units(opt$file)
	}, error = function(e) {
		message("\nError while running preprocess.R: ", e$message)
		quit(status = 1)
	}
)

quit(status = 0)
