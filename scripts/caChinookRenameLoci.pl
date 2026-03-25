#!/usr/bin/perl

#####################################################################################################
## loci are renamed using information from:
## https://github.com/eriqande/california-chinook-microhaps/blob/main/inputs/Calif-Chinook-Amplicon-Panel-Information.csv
#####################################################################################################

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'hf:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# parse the command line
my( $file ) = &parsecom( \%opts );

my %hash = (
	'tag_id_1126' => 'Ots_mhap001_01_13123559',
	'Ots_myoD-364' => 'Ots_scon001_01_14555984',
	'Ots_129458-451' => 'Ots_scon002_01_18936635',
	'tag_id_2_911' => 'Ots_mhap002_01_25667029',
	'tag_id_32' => 'Ots_mhap003_01_38431263',
	'tag_id_2_705' => 'Ots_mhap004_01_39070917',
	'tag_id_278' => 'Ots_mhap005_01_42560204',
	'Ots_Prl2' => 'Ots_scon003_01_47781127',
	'Ots_107806-821' => 'Ots_scon004_01_58114002',
	'tag_id_2_419' => 'Ots_mhap006_01_58405385',
	'tag_id_1363' => 'Ots_mhap007_01_68692192',
	'Ots_111312-435' => 'Ots_scon005_02_20523765',
	'tag_id_695' => 'Ots_mhap008_02_25262616',
	'tag_id_2_855' => 'Ots_mhap009_02_26540846',
	'Ots_105407-117' => 'Ots_scon006_02_27383109',
	'Ots_128302-57' => 'Ots_scon007_02_31507195',
	'tag_id_2_1586' => 'Ots_mhap010_03_09269999',
	'tag_id_2_284' => 'Ots_mhap011_03_14901098',
	'tag_id_282' => 'Ots_mhap012_03_24780480',
	'Ots_94857-232' => 'Ots_scon008_03_27932737',
	'Ots_112301-43' => 'Ots_scon009_03_29061420',
	'tag_id_2_3094' => 'Ots_mhap013_03_39003492',
	'Ots_118938-325' => 'Ots_scon010_03_43608004',
	'tag_id_2_234' => 'Ots_mhap014_03_45945016',
	'tag_id_2_311' => 'Ots_mhap015_03_54344571',
	'NC_037099.1:62937268-62937373' => 'Ots_vgll001_03_57888689',
	'Ots_100884-287' => 'Ots_scon011_03_69168017',
	'Ots_107285-93' => 'Ots_scon012_04_12764504',
	'Ots_105401-325' => 'Ots_scon013_04_19552989',
	'Ots_117043-255' => 'Ots_scon014_04_20534278',
	'tag_id_2_1539' => 'Ots_mhap016_04_34492892',
	'Ots_96500-180' => 'Ots_scon015_04_51066873',
	'tag_id_2_2222' => 'Ots_mhap017_04_52951178',
	'Ots_aspat-196' => 'Ots_scon016_04_65335753',
	'Ots_111681-657' => 'Ots_scon017_05_15779620',
	'tag_id_2_3026' => 'Ots_mhap018_05_31573007',
	'OkiOts_120255-113' => 'Ots_coho001_05_32691399',
	'Ots_BMP-2-SNP1' => 'Ots_scon018_05_36663656',
	'tag_id_2_633' => 'Ots_mhap019_05_38661515',
	'tag_id_384' => 'Ots_mhap020_05_39928311',
	'tag_id_3920' => 'Ots_mhap021_05_43678882',
	'tag_id_3221' => 'Ots_mhap022_05_44270127',
	'tag_id_2_123' => 'Ots_mhap023_05_44336205',
	'Ots_u4-92' => 'Ots_scon019_05_54626401',
	'Ots_127236-62' => 'Ots_scon020_05_58077375',
	'Ots_SClkF2R2-135' => 'Ots_scon021_05_65590022',
	'tag_id_423' => 'Ots_mhap024_05_72474583',
	'tag_id_1554' => 'Ots_mhap025_05_73286856',
	'Ots_u07-07.161' => 'Ots_scon022_05_75986890',
	'tag_id_744' => 'Ots_mhap026_06_30779074',
	'tag_id_2_414' => 'Ots_mhap027_06_34566531',
	'tag_id_120' => 'Ots_mhap028_06_40413020',
	'tag_id_650' => 'Ots_mhap029_06_40662054',
	'tag_id_1030' => 'Ots_mhap030_06_44208207',
	'tag_id_2_694' => 'Ots_mhap031_06_55836851',
	'tag_id_2_1579' => 'Ots_mhap032_06_78042744',
	'tag_id_1551' => 'Ots_mhap033_07_02298130',
	'Ots_128757-61' => 'Ots_scon023_07_05760620',
	'tag_id_871' => 'Ots_mhap034_07_25250444',
	'Ots_131460-584' => 'Ots_scon024_07_29861576',
	'tag_id_819' => 'Ots_mhap035_07_39787654',
	'tag_id_1079' => 'Ots_mhap036_07_50103482',
	'tag_id_2_700' => 'Ots_mhap037_07_53954141',
	'Ots_97077-179' => 'Ots_scon025_07_55870548',
	'Ots_105105-613' => 'Ots_scon026_07_62479561',
	'NC_037104.1:56552952-56553042' => 'Ots_sixx001_08_26821599',
	'NC_037104.1:56552815-56552925' => 'Ots_sixx002_08_26821716',
	'NC_037104.1:55923357-55923657' => 'Ots_wrap001_08_27450181',
	'tag_id_5617' => 'Ots_mhap038_08_33665803',
	'tag_id_2_58' => 'Ots_mhap039_08_39424024',
	'tag_id_2_1016' => 'Ots_mhap040_08_48275061',
	'Ots_102213-210' => 'Ots_scon027_08_54400656',
	'Ots_101704-143' => 'Ots_scon028_08_59645920',
	'tag_id_757' => 'Ots_mhap041_08_65278571',
	'Ots_103041-52' => 'Ots_scon029_08_65607839',
	'tag_id_2_1268' => 'Ots_mhap042_08_87336385',
	'Ots_102457-132' => 'Ots_scon030_09_30976363',
	'tag_id_425' => 'Ots_mhap043_09_39126107',
	'tag_id_773' => 'Ots_mhap044_09_39729374',
	'tag_id_1413' => 'Ots_mhap045_09_53074331',
	'tag_id_2_136' => 'Ots_mhap046_09_65932978',
	'Ots_110064-383' => 'Ots_scon031_09_75568964',
	'tag_id_2_332' => 'Ots_mhap047_09_76555021',
	'Ots_111666-408' => 'Ots_scon032_09_79852211',
	'tag_id_2_206' => 'Ots_mhap048_10_22283306',
	'Ots_112820-284' => 'Ots_scon033_10_30521814',
	'Ots_102414-395' => 'Ots_scon034_10_34047049',
	'Ots_108007-208' => 'Ots_scon035_10_35062210',
	'tag_id_5720' => 'Ots_mhap049_10_36219596',
	'Ots_129144-472' => 'Ots_scon036_10_45303613',
	'tag_id_684' => 'Ots_mhap050_10_53198846',
	'tag_id_999' => 'Ots_mhap051_10_55361826',
	'Ots_110495-380' => 'Ots_scon037_10_56919994',
	'tag_id_2_749' => 'Ots_mhap052_10_58898610',
	'Ots_99550-204' => 'Ots_scon038_10_84187271',
	'Ots_128693-461' => 'Ots_scon039_11_26122246',
	'tag_id_945' => 'Ots_mhap053_11_42403691',
	'tag_id_603' => 'Ots_mhap054_12_08080631',
	'tag_id_1692' => 'Ots_mhap055_12_19688408',
	'Ots_102867-609' => 'Ots_scon040_12_28121453',
	'Ots_u07-49.290' => 'Ots_scon041_12_37084685',
	'tag_id_826' => 'Ots_mhap056_12_62690243',
	'NC_037108.1:73543706-73544006' => 'Ots_wrap002_12_70794116',
	'tag_id_235' => 'Ots_mhap057_13_38673256',
	'tag_id_2_321' => 'Ots_mhap058_13_39744969',
	'Ots_mybp-85' => 'Ots_scon042_13_40125546',
	'tag_id_275' => 'Ots_mhap059_13_46503417',
	'tag_id_1276' => 'Ots_mhap060_14_09424912',
	'tag_id_251' => 'Ots_mhap061_14_42651215',
	'Ots_104063-132' => 'Ots_scon043_15_03971044',
	'Ots_124774-477' => 'Ots_scon044_15_09463688',
	'Ots_102801-308' => 'Ots_scon045_15_16260319',
	'Ots_104569-86' => 'Ots_scon046_15_19381451',
	'tag_id_430' => 'Ots_mhap062_15_25679743',
	'tag_id_2_1693' => 'Ots_mhap063_15_31506472',
	'tag_id_381' => 'Ots_mhap064_15_33672857',
	'tag_id_2_1158' => 'Ots_mhap065_15_37602330',
	'tag_id_2_978' => 'Ots_mhap066_16_09752947',
	'tag_id_1191' => 'Ots_mhap067_16_10421524',
	'tag_id_542' => 'Ots_mhap068_16_19318682',
	'NC_037112.1:24542569-24542869' => 'Ots_wrap003_16_36409831',
	'Ots_117242-136' => 'Ots_scon047_16_41265618',
	'Ots_unk_526' => 'Ots_scon048_16_41583413',
	'tag_id_186' => 'Ots_mhap069_16_52637329',
	'Ots_SWS1op-182' => 'Ots_scon049_17_08781484',
	'Ots_129170-683' => 'Ots_scon050_17_13921807',
	'tag_id_600' => 'Ots_mhap070_18_15451213',
	'Ots_123921-111' => 'Ots_scon051_18_24158022',
	'Ots_S71-336' => 'Ots_scon052_18_32018111',
	'tag_id_2_1382' => 'Ots_mhap071_18_36182429',
	'tag_id_1243' => 'Ots_mhap072_19_00235863',
	'tag_id_1872' => 'Ots_mhap073_19_03686095',
	'tag_id_2_40' => 'Ots_mhap074_19_07050017',
	'Ots_CD59-2' => 'Ots_scon053_19_07690242',
	'tag_id_2_188' => 'Ots_mhap075_19_12460365',
	'tag_id_2_953' => 'Ots_mhap076_19_21708672',
	'tag_id_2_859' => 'Ots_mhap077_19_31688140',
	'tag_id_1144' => 'Ots_mhap078_20_03804523',
	'tag_id_1733' => 'Ots_mhap079_20_12339375',
	'tag_id_2_1887' => 'Ots_mhap080_20_25674590',
	'Ots_AsnRS-60' => 'Ots_scon054_20_45928282',
	'Ots_AldoB4-183' => 'Ots_scon055_21_06327095',
	'Ots_AldB1-122' => 'Ots_scon056_21_06330148',
	'Ots_105132-200' => 'Ots_scon057_21_10893989',
	'tag_id_2_98' => 'Ots_mhap081_21_11344019',
	'Ots_123048-521' => 'Ots_scon058_22_09128987',
	'Ots_CD63' => 'Ots_scon059_22_14067771',
	'Ots_110551-64' => 'Ots_scon060_22_15288100',
	'tag_id_2_939' => 'Ots_mhap082_22_35581634',
	'tag_id_2_935' => 'Ots_mhap083_24_18120751',
	'Ots_106747-239' => 'Ots_scon061_24_18987029',
	'tag_id_787' => 'Ots_mhap084_24_26153777',
	'Ots_107074-284' => 'Ots_scon062_25_09879096',
	'tag_id_2_502' => 'Ots_mhap085_25_12311919',
	'tag_id_2_2632' => 'Ots_mhap086_25_14239425',
	'tag_id_2_661' => 'Ots_mhap087_25_18268602',
	'tag_id_1425' => 'Ots_mhap088_25_40766038',
	'tag_id_4969' => 'Ots_mhap089_26_09047824',
	'tag_id_70' => 'Ots_mhap090_26_21362428',
	'tag_id_2_3452' => 'Ots_mhap091_26_24638247',
	'Ots_112876-371' => 'Ots_scon063_26_29909709',
	'tag_id_716' => 'Ots_mhap092_26_33509326',
	'tag_id_3194' => 'Ots_mhap093_27_16923862',
	'tag_id_2_487' => 'Ots_mhap094_28_07482658',
	'NC_037124.1:12267397-12267697' => 'Ots_rosa001_28_13452385',
	'winter-greb1l-diagnostic_09_12269030_151' => 'Ots_rosa002_28_13453826',
	'winter-greb1l-diagnostic_10_12269914_151' => 'Ots_rosa003_28_13454761',
	'NC_037124.1:12270118-12270418' => 'Ots_rosa004_28_13455044',
	'NC_037124.1:12272852-12273152' => 'Ots_rosa005_28_13457813',
	'winter-greb1l-diagnostic_13_12275053_151' => 'Ots_rosa006_28_13459923',
	'NC_037124.1:12277401-12277701_Tasha-SNP-1' => 'Ots_rosa007_28_13462389',
	'NC_037124.1:12279142-12279478' => 'Ots_rosa008_28_13464137',
	'NC_037124.1:12281207-12281551' => 'Ots_rosa009_28_13466193',
	'NC_037124.1:12310649-12310949_Tasha-SNP-2' => 'Ots_rosa010_28_13495518',
	'Ots_122414-56' => 'Ots_scon064_28_29123949',
	'tag_id_1470' => 'Ots_mhap095_29_04417406',
	'tag_id_2_3471' => 'Ots_mhap096_29_09123413',
	'Ots_117432-409' => 'Ots_scon065_29_09282633',
	'Ots_96222-525' => 'Ots_scon066_29_11441491',
	'Ots_108390-329' => 'Ots_scon067_29_15275376',
	'tag_id_481' => 'Ots_mhap097_29_21839563',
	'tag_id_664' => 'Ots_mhap098_30_11148896',
	'Ots_EP-529' => 'Ots_scon068_30_16694275',
	'Ots_130720-99' => 'Ots_scon069_30_18929423',
	'tag_id_2_2741' => 'Ots_mhap099_30_24390612',
	'Ots_PGK-54' => 'Ots_scon070_30_24794215',
	'tag_id_1281' => 'Ots_mhap100_30_47727644',
	'Ots_113242-216' => 'Ots_scon071_31_05296625',
	'Ots_101119-381' => 'Ots_scon072_31_13703251',
	'Ots_109693-392' => 'Ots_scon073_31_24541562',
	'Ots_118175-479' => 'Ots_scon074_32_08224124',
	'Ots_112419-131' => 'Ots_scon075_32_10523540',
	'tag_id_2_20' => 'Ots_mhap101_33_00388829',
	'Ots_NAML12_1-SNP1' => 'Ots_scon076_33_09788425',
	'tag_id_2_786' => 'Ots_mhap102_33_23745808',
	'Ots_118205-61' => 'Ots_scon077_33_31230848',
	'Ots_106499-70' => 'Ots_scon078_33_45793420',
	'tag_id_2_2787' => 'Ots_mhap103_33_48609049',
	'NC_037130.1:864908-865208' => 'Ots_lfar001_34_00954022',
	'NC_037130.1:1062935-1063235' => 'Ots_lfar002_34_01151770',
	'tag_id_669' => 'Ots_mhap104_34_02236036',
	'tag_id_427' => 'Ots_mhap105_34_02392610',
	'tag_id_2_9' => 'Ots_mhap106_34_09985671',
	'sdy_I183' => 'Ots_sexy001_NW024608692.1_00004816'
);

my @lines; # holds lines from input .csv file

# read in csv file
&filetoarray( $file, \@lines );

my $header = shift( @lines );

my @head = split( /,/, $header );
for( my $i=0; $i<@head; $i++ ){
	if( $head[$i] =~ /(_[12]{1}$)/ ){
		my $newstring = substr($head[$i], 0, -2);
		if( exists( $hash{$newstring} ) ){
			$head[$i] = join( '', $hash{$newstring}, $1 );

		}
	}
	#print $head[$i], "\n";
}

my $newheader = join( ',', @head );
unshift( @lines, $newheader );

my @fileArr = split( /\./, $file );
my $ext = pop( @fileArr );
push( @fileArr, "lociRenamed" );
push( @fileArr, $ext );
my $out = join( '.', @fileArr );

open( OUT, '>', $out ) or die "Can't open $out: $!\n\n";
foreach my $line( @lines ){
	print OUT $line, "\n";
}
close OUT;

#print Dumper( \@lines );

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\ncaChinookRenameLoci.pl Program Options:\n";
  print "\t\t[ -h | -f ]\n\n";
  print "\t-h:\tDisplay this help message.\n";
  print "\t\tThe program will die after the help message is displayed.\n\n";
  print "\t-f:\tSpecify the name of the input file.\n";
  print "\t\tThis should be the final haplotype .csv file output by the ca_chinook microhaplotopia.R script.\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
  my( $params ) =  @_;
  my %opts = %$params;
  
  # set default values for command line arguments
  my $file = $opts{f} || die "No input .csv file specified.\n\n"; #used to specify input file name.

  return( $file );

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
