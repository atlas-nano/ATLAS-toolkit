#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester GetSelections);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use ManipAtoms qw(GetAtmList);

sub init;
sub numerically { ($a<=>$b); }
sub getOffsets;
sub reorderAtoms;
sub updateResNum;

my ($bgfFile, $saveFile, $selection);
my ($ATOMS, $BONDS, $HEADERS, $SELECT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECT = GetAtmList($SELECT, $ATOMS);
die "ERROR: Not a valid atom selection!\n" if (! keys %{ $SELECT });
print "Done\nReordering atoms...";
&reorderAtoms($ATOMS, $BONDS, $SELECT);
&updateResNum($ATOMS);
print "Done\nCreating BGF $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateResNum {
    my ($atoms) = $_[0];
    my (@atmList, $i, $RES, $resCount, $resName);

    @atmList = sort numerically keys %{ $atoms };
    $resCount = $RES->{NUM} = $atoms->{$atmList[0]}{RESNUM};
    for $i (@atmList) {
	if ($atoms->{$i}{RESNAME} ne $RES->{NAME} or $atoms->{$i}{RESNUM} != $RES->{NUM}) {
	    $resCount++;
	    $RES->{NAME} = $atoms->{$i}{RESNAME};
	    $RES->{NUM} = $atoms->{$i}{RESNUM};
	}
	$atoms->{$i}{RESNUM} = $resCount;

    }
}

sub reorderAtoms {
    my ($atoms, $bonds, $atmSelection) = @_;
    my (@atmList, $atom, $i, $oldBonds, $atmOffset); 
    my ($oldAtoms, @revAtmList, $j, $resOffset, %MAP);

    @atmList = sort numerically %{ $atmSelection };
    for $i (@atmList) {
	%{ $oldAtoms->{$i} } = %{ $atoms->{$i} };
	@{ $oldBonds->{$i} } = @{ $bonds->{$i} } if ($bonds->{$i});
    }
    
    ($atmOffset, $resOffset) = getOffsets($atoms, \@atmList);

    @revAtmList = reverse @atmList;
    for $i (0 .. $#atmList) {
	$MAP{$atmList[$i]} = $revAtmList[$i];
    }
    for $i (0 .. $#revAtmList) {
	$atom = $revAtmList[$i];
	$atoms->{$atom} = $oldAtoms->{ $atmList[$i] };
	#$atoms->{$atom}{RESNUM} = $resOffset - $atoms->{$atom}{RESNUM};
	@{ $bonds->{$atom} } = @{ $oldBonds->{ $atmList[$i] } };
	for $j (0 .. $#{ $bonds->{$atom} }) {
	    $bonds->{$atom}[$j] = $MAP{$bonds->{$atom}[$j]};
	}
    }
}

sub getOffsets {
    my ($atoms, $atmList) = @_;
    my ($atmOffset, $firstAtm, $lastAtm);
    my ($resOffset, $firstRes, $lastRes);

    $firstAtm = $atmList->[0];
    $lastAtm = $atmList->[$#{ $atmList }];
    $atmOffset = $lastAtm - $firstAtm - 1;
    $firstRes = $atoms->{$firstAtm}{RESNUM};
    $lastRes = $atoms->{$lastAtm}{RESNUM};
    $resOffset = $lastRes-$firstRes-1;
    
    return ($atmOffset, $resOffset);
}

sub init {
    my (%OPTS, $atomSelect);
    getopt('bas',\%OPTS);
    die "usage: $0 -b bgffile -a [atom selection] -s [save name]\n" if (! exists($OPTS{b}));
    print "Initialzing...";
    ($bgfFile, $saveFile, $selection) = ($OPTS{b}, $OPTS{s}, $OPTS{a});
    FileTester($bgfFile);
    $atomSelect = "*" if (! defined($selection));
    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }
 
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_reordered.bgf";
    }
    print "Done\n";
}
