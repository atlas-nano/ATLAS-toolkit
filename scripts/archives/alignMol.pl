#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::General qw(FileTester CoP Rotate DotProduct);
use Packages::FileFormats qw(GetPDBFileInfo sortByRes createPDB);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Math::Trig qw(acos);

sub init;
sub alignMol;
sub createPDB;
sub getBaseCoG;

my ($pdbFile, $bases, $axis, $saveFile);
my ($ATOMS, $bCENTER, $RES);

$|++;
$bases = &init;
print "Parsing PDB file $pdbFile...";
$ATOMS = GetPDBFileInfo($pdbFile);
$RES = sortByRes($ATOMS);
print "Done\nAligning molecule to $axis axis...";
$bCENTER = getBaseCoG($ATOMS, $bases, $RES);
alignMol($ATOMS, $bCENTER, $axis);
print "Done\nCreating PDB file $saveFile...";
createPDB($ATOMS, undef, $saveFile);
print "Done\n";

sub alignMol {
    my ($atoms, $basesCenter, $axis) = @_;
    my (@dim, $angle, $DISP, $i, $j, $r, $AXISPOS, @angles);
    my ($phi, $theta);

    my ($CENTER) = CoP($atoms);

    @dim = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (@dim) {
	$DISP->{$i} = $basesCenter->{$i} - $CENTER->{$i};
    }

    for $j (keys %{ $atoms }) {
	for $i (@dim) {
	    $atoms->{$j}{$i} -= $CENTER->{$i};
	}
    }
    
    # spherical coordinates
    $r = sqrt($DISP->{XCOORD}**2 + $DISP->{YCOORD}**2 + $DISP->{ZCOORD}**2);
    $phi = acos($DISP->{ZCOORD}/$r);
    $theta = atan2($DISP->{YCOORD},$DISP->{XCOORD});

    if ($axis eq "ZCOORD") {
    } elsif ($axis eq "YCOORD") {
    } else {
    }

    Rotate($atoms, \@angles, 3);

    for $r (keys %{ $atoms }) {
	for $i (@dim) {
	    $atoms->{$r}{$i} += $CENTER->{$i};
	}
    }
}

sub getBaseCoG {
    my ($atoms, $baseInfo, $resInfo) = @_;
    my ($i, $j, $counter, $BASES, $currRes, $baseAtoms);

    $counter = 1;
    for $i (1 .. 2) {
	$currRes = $baseInfo->{$i};
	for $j (keys %{ $resInfo->{$currRes}{ATOMS} }) {
	    $baseAtoms->{$counter} = $ATOMS->{$j};
	    $counter++;
	}
	$BASES = CoP($baseAtoms);
    }

    return $BASES;
}

sub init {
    my (%OPTS, %BASES, $bp);

    getopt('psba', \%OPTS);
    ($pdbFile, $bp, $axis, $saveFile) = ($OPTS{p},$OPTS{b},$OPTS{a},$OPTS{s});
    for ($pdbFile, $bp, $axis) {
	die "usage: $0 -p pdb files -b \"base1:base2\" -a (x|y|z axis) -s [saveName]\n"
	    if (! defined($_));
    }

    print "Initializing...";
    FileTester($pdbFile);
    if ($bp =~ /(\d+):(\d+)/) {
	$BASES{1} = $1;
	$BASES{2} = $2;
    } else {
	die "ERROR: Basepair must be specified as Res1:Res2! Got \"$bp\"\n";
    }
    if ($axis =~ /(x|y|z)/i) {
	$axis = uc($1);
    } else {
	die "ERROR: Expected x|y|z for axis. Got \"$axis\"\n";
    }
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//; $saveFile .= "_" . lc($axis) . "align.pdb";
    }
    $axis .= "COORD";
    return \%BASES;
}
