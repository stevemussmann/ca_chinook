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
getopts( 'c:g:hi:m:o:p:r:R:s:u:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $imp, $sum, $mac, $pan, $gtn, $usr, $chn, $rl1, $rl2, $out ) = &parsecom( \%opts );

# make sure machine type is valid
$mac = lc($mac);
if( $mac ne "nextseq" and $mac ne "miseq" ){
	die "Machine type must be 'nextseq' or 'miseq' only.\n\n";
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
shift( @impLines );
shift( @sumLines );

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
[Data],,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n";

foreach my $line( @impLines ){
	# format CH numbers and increment
	my $fmtnum = sprintf( "%06d", $num);
	$num++;

	# remove spaces from lines
	$line =~ s/ /_/g;
	my @temp = split( /,/, $line );

	# make sure negatives are uniquely named
	if( $temp[11] =~ /^Negative/ ){
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

	# print Sample_Well,I7_Index_ID,index,I5_Index_ID,
	print OUT "$well,$temp[6],$temp[7],$temp[3],";

	my $i5index;
	if( $mac eq "nextseq" ){
		$i5index = $temp[5];
	}elsif( $mac eq "miseq" ){
		$i5index = $temp[4];
	}else{
		die "Must specify nextseq or miseq with -m option.\n\n";
	}

	# print index2,Sample_Project,Description\n
	print OUT "$i5index,$projHash{$temp[6]},\"$pan\"\n";

}

print("\nOutput written to $out.\n\n");

#print Dumper( \%prjs );

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
  print "\t\t[ -c | -g | -h | -i | -m | -o | -p | -r | -R | -s | -u ]\n\n";
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

  return( $imp, $sum, $mac, $pan, $gtn, $usr, $chn, $rl1, $rl2, $out );

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
    next if($line =~ /^\s*$/);
    push( @$array, $line );
  }
  
  # close input file
  close FILE;

}

#####################################################################################################
