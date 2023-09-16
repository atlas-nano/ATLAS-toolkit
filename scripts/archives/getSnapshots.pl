#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use warnings;
no warnings "recursion";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Packages::AMBER qw(ParseAmberTrj getTopInfo getOpts getConn GetAmberByteOffset);
use Packages::General qw(FileTester IsInteger GetSoluteAtoms CoM CenterOnMol);
use Packages::FileFormats qw(createBGF createHeaders addHeader GetBondList);
use Packages::GetParms;
use Packages::ManipAtoms qw(ImageAtoms);

sub init;
sub showUsage;
sub readEngFile;
sub numerically {($a<=>$b);}
sub makeBox;
sub getAtoms;
sub createEngProfile;

my ($DATA, $totAtms, $trjFile, $tot, $savePrefix, $molList); 
my ($SELECT, $engFile, $PARMS, $pS, $BONDS, $SOLUTEATMS);

$|++;
&init;
$SELECT = readEngFile($engFile, $tot);
print "Creating energy profile ";
&createEngProfile($SELECT, $savePrefix);
print "Done\n";
$pS = "Parsing AMBER traj $trjFile...";
&GetAmberByteOffset($SELECT, $trjFile, $totAtms);
ParseAmberTrj($DATA->{ATOMS}, $trjFile, $SELECT, $totAtms, \&saveSnapshot, $pS, undef);

sub createEngProfile {
    my ($selections, $engFile) = @_;
    
    $engFile .= "_eng_profile.dat";

    print "$engFile...";
    open ENGFILE, "> $engFile" || die "ERROR: Cannot create $engFile: $!\n";
    for (sort numerically keys %{ $selections }) {
	printf ENGFILE "%8d%12.5f\n",$_, $selections->{$_};
    }
    close ENGFILE;
}

sub saveSnapshot {
    my ($ATOMS, $BOX, $currFrame, $junk) = @_;
    my ($HEADERS, $MOL, $CENTER, $i, $BBOX);
    my ($bgfFile) = "${savePrefix}_${currFrame}.bgf";

    $MOL = getAtoms($ATOMS, $SOLUTEATMS);
    $CENTER = CoM($MOL);
    CenterOnMol($ATOMS, $CENTER);

    $BBOX = makeBox($BOX);
    for $i (keys %{ $molList }) {
	$MOL = getAtoms($ATOMS, $molList->{$i});
	$CENTER = CoM($MOL);
	ImageAtoms($ATOMS, $molList->{$i}, $CENTER, $BBOX);
    }

    $HEADERS = createHeaders($BOX, $bgfFile);
    addHeader($ATOMS, $HEADERS);
    createBGF($ATOMS, $BONDS, $bgfFile);
}
 
sub readEngFile {
    my ($fileName, $numberToSave) = @_;
    my ($i, @tmp, %ENG, %RET, $tot);

    open ENGFILE, $fileName || die "ERROR: Cannot open $fileName: $!\n";
    while (<ENGFILE>) {
	chomp;
	if ($_ =~ /^\s+(\d+)\s+(\-\d+\.\d+)/) {
	    $ENG{$2} = $1;
	}
    }
    close ENGFILE || die "ERROR: Cannot close $fileName: $!\n";

    die "ERROR: $fileName is invalid!\n" if (! keys %ENG);
    @tmp = sort numerically keys %ENG;
    $numberToSave = scalar @tmp if (scalar(@tmp) < $numberToSave);

    for $i (1 .. $numberToSave) {
	$tot = shift @tmp;
	$RET{$ENG{$tot}} = $tot;
    }

    return \%RET;
}

sub getAtoms {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);

    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = $allAtoms->{$i};
    }

    return \%ATOMS;
}

sub makeBox {
    my ($box) = $_[0];
    my (%BOX, $i, @dim);

    @dim = ("X", "Y", "Z");
    for $i (0 .. 2) {
        $BOX{$dim[$i]}{lo} = 0;
        $BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{lo} = 0;
        $BOX{$dim[$i] . "COORD"}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i] . "COORD"}{len} = $box->{$i + 2}{DATA};
    }

    return \%BOX;
}

sub init {
    my (%OPTS, $usage, $parmFile, $amberOPTS, $topFile);
    
    getopt('pets',\%OPTS);
    ($parmFile, $engFile, $tot, $savePrefix) = ($OPTS{p}, $OPTS{e}, $OPTS{t}, $OPTS{s});
    $usage = &showUsage;
    for ($parmFile, $engFile, $tot) {
	die "$usage\n" if (! defined($_));
    }
    print "Initializing...";
    FileTester($parmFile);
    FileTester($engFile);
    die "ERROR: Expected positive integer for number of snapshots, got \"$tot\"!\n"
	if (! IsInteger($tot) || $tot < 1);

    print "getting parms...";
    $PARMS = Packages::GetParms->new();
    die "Error in Paramater file\n" if (! $PARMS->IsValidParams($parmFile));
    
    $topFile = $PARMS->{Files}{topology};
    $trjFile = $PARMS->{Files}{trajectory};
    FileTester($topFile);
    FileTester($trjFile);
    
    print "parsing $topFile...";
    $amberOPTS = &getOpts;
    ($DATA, $totAtms) = getTopInfo($topFile, $amberOPTS);
    $BONDS = getConn($DATA->{BONDLIST}, $DATA->{ATOMS});
    $molList = GetBondList($DATA->{ATOMS}, $BONDS);
    $SOLUTEATMS = GetSoluteAtoms($DATA->{ATOMS}, $molList);

    if (! defined($savePrefix)) {
	$savePrefix  =  "./" . basename($trjFile);
	$savePrefix =~ s/\.\w+//;
    }
    print "${totAtms} atoms...Done\n";
}

sub showUsage {
    my ($usage) = "usage: $0 -p parameter file -e energy file -t number of snapshots -s [saveprefix]\n" .
	"options:\n" .
	"-p parameter file: location of the parameter file. This file has locations of trajectory/topology files\n" . 
	"-e energy file: file with two columns: snapshot number and snapshot energy\n" .
	"-t number of snapshots: the number of snapshots to create. Must be > 0 and <= total snapshots\n" .
	"-s [saveprefix]: the prefix of the bgf files created. optional\n";
	
	return $usage;
}
