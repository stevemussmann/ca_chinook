#!/usr/bin/env Rscript

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

source("/home/mussmann/local/src/microhap/mega-simple-microhap-snakeflow/preprocess/sample-sheet-processing-functions.R")

create_samples_and_units(opt$file)

quit()
