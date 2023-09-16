#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use AMBER qw(getTopInfo ParseAmberTrj getOpts getConn GetAmberByteOffset);
use General qw(FileTester TrjSelections);
use FileFormats qw(createBGF createHeaders addHeader);
use BOX qw(MakeBox);
use File::Basename;

sub init;
sub numerically { ($a<=>$b); }
sub createFile;
sub getPrec;

die "usage: $0 topFile trjFile \"trajectory selection\" [savePrefix]\n"
    if (! @ARGV || $#ARGV < 3);

my ($topFile, $trjFile, $selection, $savePrefix) = @ARGV;
my ($OPTS, $SELECT, $DATA, $totAtms, $precision);
my ($printStr, $VEC, $link, $atmCounter, $trjName);

$|++;
print "Initializing...";
&init;
print "Done\nParsing AMBER topology file $topFile...";
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
print "Done\nGenerating Connectivities...";
$DATA->{BONDS} = getConn(\%{ $DATA->{"BONDLIST"} }, $DATA->{"ATOMS"});
print "Done\n";
$printStr = "Parsing AMBER trajectory $trjFile...";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
$precision = getPrec($SELECT);
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&createFile, $printStr, $DATA->{BONDS});

sub init {
    FileTester($topFile);
    FileTester($trjFile);

    $OPTS = &getOpts;
    $savePrefix = basename($trjFile) if (! defined($savePrefix));
    $savePrefix =~ s/\.\w+$//;
 
    if (-e $selection) {
	open SEL, $selection or die "ERROR: Cannot open file $selection: $!\n";
	$selection = "";
	while (<SEL>) {
	    chomp;
	    $selection .= "$_ ";
	}
	close SEL;
    }   
    $SELECT = TrjSelections($selection);
}

sub createFile {
    my ($ATOMS, $BBOX, $frameNum, $BONDS) = @_;
    my ($i, $BOX, $counter, @tmp);

    $frameNum = sprintf("%0" . $precision . "d", $frameNum);
    my ($bgfFile) = "${savePrefix}_${frameNum}.bgf";

    $BOX->{1}{DATA} = 90;
    $counter = 2;
    @tmp = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (0 .. $#tmp) {
	$BOX->{$counter}{DATA} = $BBOX->{($i + 2)}{DATA};
	$counter++;
    }
    $BOX = MakeBox($BOX);
    my ($HEADERS) = createHeaders($BOX, $bgfFile);
    addHeader($ATOMS, $HEADERS);
    createBGF($ATOMS, $BONDS, $bgfFile);
}

sub getPrec {
   my ($trjSelect) = @_;
   my ($prec, @trjList);

   $prec = 1;
   @trjList = sort numerically keys %{ $trjSelect };
   $prec = length($trjList[$#trjList]);

   return $prec;
}

