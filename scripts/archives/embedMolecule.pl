#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester CoP CoM CombineMols);
use Packages::FileFormats qw(GetBGFFileInfo createBGF createHeaders addHeader addBoxToHeader);
use Packages::CERIUS2 qw(parseCerius2FF LoadFFs);
use Packages::BOX qw(GetBox);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub checkAtomTypes;
sub embeddMols;

my ($DATA, $FFILES, $ATOMS, $BONDS, $HEADERS, $BOX, $PARMS);
my ($i, $soluBGF, $centerType, $centerFunc, $solvBGF, $saveName);

$|++;
$FFILES = &init;
$PARMS = LoadFFs($FFILES);
print "Parsing solvent/membrane bgf $solvBGF...";
($DATA->{SOLVENT}{ATOMS}, $DATA->{SOLVENT}{BONDS}, $DATA->{SOLVENT}{HEADERS}) = GetBGFFileInfo($solvBGF, 1);
print "Done\nParsing solute bgf $soluBGF...";
($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLUTE}{BONDS}, $DATA->{SOLUTE}{HEADERS}) = GetBGFFileInfo($soluBGF, 1);
print "Done\n";
checkAtomTypes($DATA, $PARMS);
for $i ("SOLUTE", "SOLVENT") {
    print "Computing ${centerType} for $i...\r";
    $DATA->{$i}{CENTER} = $centerFunc->($DATA->{$i}{ATOMS});
}
print "Computing ${centerType}s...Completed\n";
print "Embedding systems...";
embeddMols($DATA);
($ATOMS, $BONDS) = CombineMols($DATA->{SOLUTE}{ATOMS}, $DATA->{SOLVENT}{ATOMS},
			       $DATA->{SOLUTE}{BONDS}, $DATA->{SOLVENT}{BONDS});
$BOX = GetBox($ATOMS, $PARMS, undef);
$HEADERS = createHeaders(undef, $saveName);
addBoxToHeader($HEADERS, $BOX);
addHeader($ATOMS, $HEADERS);
print "Done\nCreating BGF file $saveName...";
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub init {
    my (%OPTS, @tmp, @FFILES, $i, $FF);
    
    getopt('mswfc',\%OPTS);
    ($solvBGF, $soluBGF, $FF, $centerType, $saveName) = ($OPTS{m}, $OPTS{s}, $OPTS{f}, $OPTS{c}, $OPTS{w});
    
    die "usage: $0 -m solvent/membrane bgf -s solute bgf -f \"forcefield1 forcefield2...\"\n" . 
	"\t-d [x|y|z] constant dimension -c [com|cog] centering type -w [saveName.bgf]\n"
	if (! defined($solvBGF) || ! defined($soluBGF) || ! defined($FF));
    
    print "Initializing...";
    FileTester($solvBGF);
    FileTester($soluBGF);

    $centerType = "COM" if (! defined($centerType));
    $centerType = uc($centerType);

    if ($centerType eq "COM") {
	$centerFunc = \&CoM;
    } else {
	$centerFunc = \&CoP;
    }

    if (! defined($saveName)) {
	$saveName = basename($soluBGF);
	$saveName =~ s/\.\w+$/_embed\.bgf/;
    }

    if ($FF =~ /\s+/) {
        @tmp = split /\s+/, $FF;
    } else {
        @tmp = ($FF);
    }
  
    for $i (@tmp) {
        if (-e $i &&  -r $i && -T $i) {
            push @FFILES, $i;
	}
    }

    die "ERROR: No valid CERIUS2 forcefield found!\n" if (! @FFILES);
    print "Done\n";
    return \@FFILES;
}

sub checkAtomTypes {
    my ($mols, $parms) = @_;
    my ($i, $j, $ffType);

    for $i ("SOLUTE", "SOLVENT") {
	for $j (keys %{ $mols->{$i}{ATOMS} }) {
	    $ffType = $mols->{$i}{ATOMS}{$j}{FFTYPE};
	    die "ERROR: Force field type $ffType not found in forcefield(s). Aborting\n" 
		if (! exists($parms->{ATOMTYPES}{$ffType}));
	    $mols->{$i}{ATOMS}{$j}{MASS} = $parms->{ATOMTYPES}{$ffType}{MASS};
	}
    }
}

sub embeddMols {
    my ($mols) = $_[0];
    my ($i, $j, @dims);

    @dims = keys %{ $mols->{SOLVENT}{CENTER} };
    for $i (@dims) {
	$mols->{SOLUTE}{OFFSET}{$i} = $mols->{SOLUTE}{CENTER}{$i} - $mols->{SOLVENT}{CENTER}{$i};
    }

    for $i (keys %{ $mols->{SOLUTE}{ATOMS} }) {
	for $j (@dims) {
	    $mols->{SOLUTE}{ATOMS}{$i}{$j} -= $mols->{SOLUTE}{OFFSET}{$j};
	}
    }
}
