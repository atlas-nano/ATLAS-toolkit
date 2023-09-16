#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester TrjSelections);
use Math::Complex;  # The roots may be complex numbers.
use Math::Polynomial::Solve qw(cubic_roots GetLeastSquaresPlane);
use Math::Trig;
use AMBER qw(getTopInfo ParseAmberTrj GetAmberByteOffset);
use FileFormats qw (GetBGFFileInfo sortByRes);
use File::Basename;
use constant PI => atan2(1,1) * 4;

sub init;
sub calcDIHDR;
sub doAnal;
sub getLeastSquaresFit;
sub getMin;
sub numerically { ($a<=>$b); }
sub printHeader;

die "usage: $0 bgfFile parmFile amberTrj \"trajectory selection\" [saveName]\n"
    if (! @ARGV || $#ARGV < 3);

my ($bgfFile, $parmFile, $amberTrj, $selection, $dihdrFile) = @ARGV;
my ($HELIX, $SELECT, $OPTS, $totAtms, $printStr);

print "Initializing...";
&init;
print "Done\n";
open DIHEDRFILE, "> $dihdrFile" or die "ERROR: Cannot create $dihdrFile: $!\n";
&printHeader(\*DIHEDRFILE, $HELIX, $bgfFile);
$printStr = "Parsing AMBER trajectory $amberTrj...";
&GetAmberByteOffset($SELECT, $amberTrj, $totAtms);
ParseAmberTrj(undef, $amberTrj, $SELECT, $totAtms, \&doAnal, $printStr, \*DIHEDRFILE);
close DIHEDRFILE;
print "Created file $dihdrFile\n";

sub printHeader {
    my ($printHandle, $helixInfo, $atoms) = @_;
    my (@tmp, $i, $j, $RES, $nameA, $nameB);

    printf $printHandle "%-8s" , "#Frame";
    @tmp = sort numerically keys %{ $helixInfo };
    for $i (0 .. $#tmp -1) {
	$nameA = $helixInfo->{$tmp[$i]}{NAME};
	for $j (($i+1) .. $#tmp) {
	    $nameB = $helixInfo->{$tmp[$j]}{NAME};
	    printf $printHandle "%12s", "${nameA}->${nameB}";
	}
    }
    printf $printHandle "\n";
}

sub init {
    my ($isValid, $helixNum, $inStr, $i, $rec, @resInfo);
    my ($ATOMS, $BONDS, @tmp);

    $|++;
    FileTester($bgfFile);
    FileTester($parmFile);
    FileTester($amberTrj);

    if (! defined($dihdrFile)) {
	$dihdrFile = $amberTrj;
	$dihdrFile =~ s/\.\w+$//;
	$dihdrFile .= "_dihdr.dat";
    }
    $SELECT = TrjSelections($selection);

    print "Parsing BGF file $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    $totAtms = scalar keys %{ $ATOMS };
    $HELIX = sortByRes($ATOMS);
    for $i (keys %{ $HELIX }) {
	@tmp = keys %{ $HELIX->{$i}{ATOMS} };
	$HELIX->{$i}{NAME} = $ATOMS->{ $tmp[0] }{CHAIN};
    }
    $isValid =  0;
    open PARMFILE, $parmFile || die "ERROR: Cannot open parameter file $parmFile: $!\n";
    while (<PARMFILE>) {
	chomp;
	if ($_ =~ /^Helix (\d): (.+)/) {
	    $helixNum = $1;
	    $inStr = $2;
	    while ($inStr =~ /(\d+)/g) {
		push @{ $OPTS->{HELIX}{$helixNum}{COM} }, $1;
		$isValid = 1;
	    }
	}
    }
    close PARMFILE;
    die "ERROR: Parameter file $parmFile is invalid!\n" if (! $isValid);
}

sub calcDIHDR {
    my ($HELIX) = $_[0];
    my ($angle, $a1, $a2, $b1, $b2, $c1, $c2);
    my ($i, $j, @tmp, $returnStr, $DATA);

    @tmp = keys %{ $HELIX };

    for $i (0 .. ($#tmp - 1)) {
	for $j (($i + 1) .. $#tmp) {

	    $a1 = $HELIX->{$tmp[$i]}{NORMAL}{a};
	    $b1 = $HELIX->{$tmp[$i]}{NORMAL}{b};
	    $c1 = $HELIX->{$tmp[$i]}{NORMAL}{c};
   
	    $a2 = $HELIX->{$tmp[$j]}{NORMAL}{a};
	    $b2 = $HELIX->{$tmp[$j]}{NORMAL}{b};
	    $c2 = $HELIX->{$tmp[$j]}{NORMAL}{c};

	    $angle = ($a1*$a2 + $b1*$b2 + $c1*$c2);
	    $angle /= (sqrt($a1**2 + $b1**2 + $c1**2) * sqrt($a2**2 + $b2**2 + $c2**2));
	    $DATA->{$tmp[$i]}{$tmp[$j]} = (acos($angle) * 180/PI);
	}
    }

    for $i (sort numerically keys %{ $DATA }) {
	for $j (sort numerically keys %{ $DATA->{$i} }) {
	    $returnStr .= sprintf("%12.1f", $DATA->{$i}{$j});
	}
    }

    return $returnStr;
}

sub doAnal {
    my ($ATOMS, $BOX, $frameNum, $fileHandle) = @_;
    my ($helixNum, $i, $dim, $PLANEATMS, $angle, $CoG, $j, $res);

    $j = 0;

    for $helixNum (keys %{ $HELIX } ) {
	$j = 0;
	$CoG = ();
	for $i (@{ $OPTS->{HELIX}{$helixNum}{COM} }) {
	    $j++;
	    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
		$CoG->{$dim} += $ATOMS->{$i}{$dim};
		$PLANEATMS->{$helixNum}{ATOMS}{$j}{$dim} = $ATOMS->{$i}{$dim};
	    }
	}
	for $dim ("XCOORD", "YCOORD", "ZCOORD") { # calculate the center of geometry
	    $CoG->{$dim} /= $j;
	}
	GetLeastSquaresPlane($PLANEATMS->{$helixNum}, $CoG);
    }

    $angle = calcDIHDR($PLANEATMS);
    printf $fileHandle "%8d$angle\n", $frameNum;
}

sub getMin {
    my ($vals) = $_[0];
    my ($min, $i);

    for $i (@{ $vals }) {
	$min = $i if (! defined($min) or $min > $i);
    }
    return $min;
}
