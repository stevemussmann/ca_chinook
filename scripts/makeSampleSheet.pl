#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Time::Piece;
use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'c:g:hi:m:o:p:r:R:s:u:v', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $imp, $sum, $mac, $pan, $gtn, $usr, $chn, $rl1, $rl2, $rev, $out ) = &parsecom( \%opts );

# make sure machine type is valid
$mac = lc($mac);
if( $mac ne "nextseq" and $mac ne "miseq" ){
	die "Machine type must be 'nextseq' or 'miseq' only.\n\n";
}
print("Preparing sample sheet for $mac.\n\n");

if( $rev eq "TRUE" ){
	print("The -v option was used; reverse complementing i7 indexes.\n\n");
}

# parse project hash information
my %projHash;
my %prjs;

# starting number for CH numbers
my $num = $chn;

# read in inputs
my @impLines;
my @sumLines;
&filetoarray( $imp, \@impLines );
&filetoarray( $sum, \@sumLines );

# remove header from sample_imports
my $impHead;
if( $impLines[0] =~ /^Sample/ ){
	$impHead = shift( @impLines );
}else{
	die "ERROR: Header line missing from $imp\n\n";
}

# check if extra lines at beginning of project_summary
while( @sumLines and $sumLines[0] !~ /^Project/ ){
	shift( @sumLines );
}

# remove header from project_summary
my $sumHead;
if( $sumLines[0] =~ /^Project/ ){
	$sumHead = shift( @sumLines );
}else{
	die "Error: Header line missing from $sum\n\n";
}

# parse @sumLines
foreach my $line( @sumLines ){
	$line =~ s/ /_/g;
	my @temp = split( /,/, $line );
	$projHash{$temp[2]} = $temp[0];
	$prjs{$temp[0]}++;
}

# get date
my $date = localtime->mdy('/');

# get project description
my @prjArr = sort keys %prjs;
my $prjdsc = join( '+', @prjArr );

# make hash for checking unique indexes
my %hashIndex;

# print output
open( OUT, '>', $out ) or die "Can't open $out: $!\n\n";

print OUT "[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Investigator Name,$usr,,,,,,,,
Experiment Name,$gtn,,,,,,,,
Date,$date,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,$prjdsc,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
$rl1,,,,,,,,,
$rl2,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,\n";

print OUT "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n";

foreach my $line( @impLines ){
	# only operate on lines with full data records
	if( $line !~ /,\#N\/A/ ){
		my @indexTemp; # for holding indexes for line

		# format CH numbers and increment
		my $fmtnum = sprintf( "%06d", $num);
		$num++;

		# remove spaces from lines
		$line =~ s/ /_/g;
		my @temp = split( /,/, $line );

		# make sure negatives are uniquely named
		if( $temp[11] =~ /^Negative/ or $temp[11] =~ /^0_/ ){
			$temp[11] =~ s/_/-/g;
		}
		# print Sample_ID,Sample_Name,Sample_Plate,
		print OUT "CH_$fmtnum,CH-$fmtnum,$temp[11],";

		# make sure well position is formatted properly
		my $row;

		if( $temp[2] =~ s/^(.)//) {
			$row = $1;
		}
		my $fmtcol = sprintf( "%02d", $temp[2]);
		my $well = join( '', $row, $fmtcol );

		if( $rev eq "TRUE" ){
			$temp[7] = reverse($temp[7]);
			$temp[7] =~ tr/ACGTacgtNn/TGCAtgcaNn/;
		}

		push( @indexTemp, $temp[7] ); #i7 sequence

		# print Sample_Well,I7_Index_ID,index,I5_Index_ID,
		print OUT "$well,$temp[6],$temp[7],$temp[3],";

		my $i5index;
		if( $mac eq "nextseq" ){
			$i5index = $temp[5];
			push( @indexTemp, $temp[5] ) #nextseq i5 sequence
		}elsif( $mac eq "miseq" ){
			$i5index = $temp[4];
			push( @indexTemp, $temp[4] ) #miseq i5 sequence
		}else{
			die "Must specify nextseq or miseq with -m option.\n\n";
		}

		# print index2,Sample_Project,Description\n
		print OUT "$i5index,$projHash{$temp[6]},\"$pan\"\n";

		# add to index hash for checking uniqueness
		my $combIndex = join( '+', @indexTemp );
		$hashIndex{$combIndex}++;
	}
}

close OUT;

# check that all indexes are unique
print("Checking indexes for duplicates (i7index + i5index)...\n");
my $dupCount = 0;
foreach my $index( sort keys %hashIndex ){
	if( $hashIndex{$index} > 1 ){
		print "Index $index is repeated $hashIndex{$index} times.\n";
		$dupCount++;
	}
}
if( $dupCount > 0 ){
	print("\nThere were $dupCount duplicated indexes found.\n");
	print("Check inputs and rerun.\n\n");
	unlink $out or warn "Could not delete $out: $!\n"; #delete $out if duplicates exist
	die "No files written - exiting program.\n\n";
}else{
	print("\nNo duplicate indexes found.\n\n");
}

print("Output written to $out.\n\n");

# if writing files for miseq, also write version that can be loaded onto miseq
if( $mac eq "miseq" ){
	# read $out
	my @miseqLines;
	&filetoarray( $out, \@miseqLines );

	# replace line 21 (header line)
	$miseqLines[21] = "Sample_ID,Sample_Name,Description,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Sample_Plate";

	# make new filename for miseq run sheet file
	my @outTemp = split( /\./, $out );
	splice(@outTemp, -1, 0, "miseqrun" );
	my $newout = join( '.', @outTemp );

	# write miseq run sheet file
	open( OUT2, '>', $newout ) or die "Can't open $newout: $!\n\n";
	foreach my $line( @miseqLines ){
		print OUT2 "$line\n";
	}
	close OUT2;

	print "Sample sheet for miseq run written to $newout.\n\n";
	print "IMPORTANT: files specific to the MiSeq option were written.\n";
	print "Use $out for running the snakemake pipeline and use $newout for loading the run on the MiSeq.\n\n";
}

#print Dumper( \%prjs );
#print Dumper( \%hashIndex );

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nmakeSampleSheet.pl is a perl script developed by Steven Michael Mussmann\n\n";
  print "To report bugs send an email to mussmann\@uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -c | -g | -h | -i | -m | -o | -p | -r | -R | -s | -u | -v ]\n\n";
  print "\t-c:\tEnter starting number for CH- numbers (default = 1).\n\n";
  print "\t-g:\tEnter the GT number (required).\n\n";
  print "\t-h:\tDisplay this help message and exit.\n\n";
  print "\t-i:\tSpecify the name of the sample_imports file (default = sample_imports.csv).\n\n";
  print "\t-m:\tSpecify the machine type (required; input must be nextseq or miseq).\n\n";
  print "\t-o:\tSpecify output file (default = SampleSheet.csv).\n\n";
  print "\t-p:\tSpecify panels (default = \"ROSA, FullPanel\").\n\n";
  print "\t-r:\tSpecify read length for Read 1 (default = 76).\n\n";
  print "\t-R:\tSpecify read length for Read 2 (default = 76).\n\n";
  print "\t-s:\tSpecify the name of the project_summary file (default = project_summary.csv).\n\n";
  print "\t-u:\tSpecify the user's initials (required).\n\n";
  print "\t-v:\tUse this option if you need to reverse complement the i7 indexes (boolean; specify '-v' as an option to turn on).\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
	my( $params ) =  @_;
	my %opts = %$params;
  
	# set default values for command line arguments
	my $chn = $opts{c} || "1"; #specify beginning CH- number.
	my $gtn = $opts{g} || die "Must specify GT number.\n\n"; #specify nextseq or miseq.
	my $imp = $opts{i} || "sample_imports.csv"; #used to specify input file name.
	my $mac = $opts{m} || die "Must specify machine type (nextseq or miseq).\n\n"; #specify nextseq or miseq.
	my $out = $opts{o} || "SampleSheet.csv"; #used to specify output file name.
	my $rl1 = $opts{r} || "76"; #specify length of Read 1.
	my $rl2 = $opts{R} || "76"; #specify length of Read 2.
	my $pan = $opts{p} || "ROSA, FullPanel"; #specify panel options
	my $sum = $opts{s} || "project_summary.csv"; #used to specify input file name.
	my $usr = $opts{u} || die "Must specify user's initials\n\n"; #used to specify input file name.

	my $rev;
	if( $opts{v} ){
		$rev = "TRUE";
	}else{
		$rev = "FALSE";
	}

	return( $imp, $sum, $mac, $pan, $gtn, $usr, $chn, $rl1, $rl2, $rev, $out );

}

#####################################################################################################
# subroutine to put file into an array

sub filetoarray{

  my( $infile, $array ) = @_;

  
  # open the input file
  open( FILE, $infile ) or die "Can't open $infile: $!\n\n";

  # loop through input file, pushing lines onto array
  while( my $line = <FILE> ){
    chomp( $line );
    next if( $line =~ /^\s*$/ );
	#next if( $line =~ /^,*$/ );
    push( @$array, $line );
  }
  
  # close input file
  close FILE;

}

#####################################################################################################
