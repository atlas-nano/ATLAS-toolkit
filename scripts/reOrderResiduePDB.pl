#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetPDBFileInfo createPDB sortByRes);
use General qw(FileTester GetSelections);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub reOrderResNumbers;
sub numerically {($a<=>$b); }
sub checkRes;

my ($pdbFile, $saveFile, $RESLIST);
my ($ATOMS, $BONDS, $HEADERS, $RES);
$|++;
&init;
print "Parsing PDB file $pdbFile...";
($ATOMS, $BONDS, $HEADERS) = GetPDBFileInfo($pdbFile);
print "Done\nReordering residue numbers...";
$RES = sortByRes($ATOMS);
&checkRes($RES, $RESLIST);
$ATOMS = &reOrderResNumbers($ATOMS, $RES, $RESLIST);
print "Done\nCreating PDB file $saveFile...";
&createPDB($ATOMS, undef, $saveFile);
print "Done\n";

sub checkRes {
    my ($atomsRes, $resList) = @_;
    my ($totAtomsRes, $totResList);

    $totAtomsRes = scalar(keys %{ $atomsRes });
    $totResList = $#{ $resList } + 1;
    die "ERROR: Number of residues in pdb file ($totAtomsRes) is not same as residue list (" . 
	"$totResList). Aborting..\n" if ($totAtomsRes != $totResList);
}

sub reOrderResNumbers {
    my ($atoms, $res, $RESLIST) = @_;
    my ($i, $j, $newAtoms, $count, $resNum);

    $resNum = 1;
    for $i (@{ $RESLIST }) {
	for $j (keys %{ $res->{$i}{ATOMS} }) {
	    $atoms->{$j}{RESNUM} = $resNum;
	}
	$resNum++;
    }

    $res = ();
    $res = sortByRes($atoms);
    $count = 1;
    for $i (sort numerically keys %{ $res }) {
	for $j (sort numerically keys %{ $res->{$i}{ATOMS} }) {
	    %{ $newAtoms->{$count} } = %{ $atoms->{$j} };
	    $count++;
	}
    }

    return $newAtoms;
}

sub init {
    my (%OPTS, $select, $atms, $resList, $i);

    getopt('prs',\%OPTS);
    die "usage: $0 -p pdb file -r \"residue list\" -s (save file name)\n" 
	if (! exists($OPTS{p}) or ! exists($OPTS{r}));
    ($pdbFile, $resList, $saveFile) = ($OPTS{p}, $OPTS{r}, $OPTS{s});
    print "Initializing...";
    FileTester($pdbFile);
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_resSorted.pdb";
    }

    while ($resList =~ /(\d+)/g) {
	push @{ $RESLIST }, $1;
    }
    print "Done\n";
}
