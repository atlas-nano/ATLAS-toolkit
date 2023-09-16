#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF sortByRes GetBGFAtoms);
use General qw(GetSelections FileTester);
use ManipAtoms qw(GetAtmList);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub getAvgAtomCharge;
sub getResData;
sub numerically { ($a<=>$b); }

my ($saveFile, $bgfFile, $selection);
my ($ATOMS, $BONDS, $HEADERS, $SELECT, $SELECTATMS, $RES);
my ($BGF, $CONS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSelecting atoms...";
$SELECT = GetSelections($selection, 0);
$SELECTATMS = GetAtmList($SELECT, $ATOMS);
$RES = sortByRes($ATOMS, $SELECTATMS);
print "Done\nUpdating residue numbers to match molecules...";
($BGF, $CONS) = getResData($ATOMS, $BONDS, $RES);
&getAvgAtomCharge($ATOMS, $RES, $BGF);
print "Done\nCreating BGF file $saveFile...";
addHeader($BGF, $HEADERS);
createBGF($BGF, $CONS, $saveFile);
print "Done\n";

sub getResData {
    my ($atoms, $bonds, $resInfo) = @_;
    my (@tmp, $resAtoms, $resBonds, $selection, $i, $resTmp);

    @tmp = keys %{ $resInfo };

    for $i (keys %{ $resInfo->{ $tmp[0] }{ATOMS} }) {
	$selection->{$i} = 1;
    }
    ($resTmp, $resBonds, undef) = GetBGFAtoms($selection, $atoms, $bonds);
    for $i (keys %{ $resTmp }) {
	%{ $resAtoms->{$i} } = %{ $resTmp->{$i} };
	$resAtoms->{$i}{CHARGE} = 0;
    }
    return ($resAtoms, $resBonds);
}


sub getAvgAtomCharge {
    my ($atoms, $resList, $resData) = @_;
    my ($i, $j, @atomIndices, $startIndex, $currIndex, $count, $tot);

    $count = 0;
    for $i (keys %{ $resList }) {
	$count++;
	@atomIndices = sort numerically keys %{ $resList->{$i}{ATOMS} };
	$startIndex = $atomIndices[0];
	for $j (@atomIndices) {
	    $currIndex = $j - $startIndex + 1;
	    $resData->{$currIndex}{CHARGE} += $atoms->{$j}{CHARGE};
	    $tot += $atoms->{$j}{CHARGE};
	}
    }

    for $i (keys %{ $resData }) {
	$resData->{$i}{CHARGE} /= $count;
    }
}

sub init {
    my ($atomSelect, %OPTS);

    getopt('bsa',\%OPTS);
    for ("b", "a") {
	die "usage: $0 -b bgf file -a atom selection -s [save name]\n" if (! defined($OPTS{$_}));
    }

    print "Initializing...";
    ($bgfFile, $saveFile, $atomSelect) = ($OPTS{b}, $OPTS{s}, $OPTS{a});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }

    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }

    print "Done\n";
}
