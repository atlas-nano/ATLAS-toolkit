#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use General qw(FileTester);
use ManipAtoms qw(GetMols GetAtmList SelectAtoms BuildAtomSelectionString);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub updateResNum;
sub numerically { ($a<=>$b); }

my ($saveFile, $bgfFile, $selection);
my ($ATOMS, $BONDS, $HEADERS, $SELECTIONS, $MOLS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
$MOLS = GetMols($ATOMS,$BONDS);
print "Done\nSelecting atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
print "Done\nUpdating residue numbers to match molecules...";
&updateResNum($ATOMS, $SELECTIONS, $MOLS);
print "Done\nCreating BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateResNum {
	my ($atoms, $selection, $mols) = @_;
	my ($i, $j, $valid);

	for $i (keys %{ $mols }) {
		$valid = 0;
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			$valid = 1 if (exists($selection->{$j}));
		}
		next if (! $valid);
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			$atoms->{$j}{RESNUM} = $i;
		}
	}
}


sub updateResNumN {
	my ($atoms, $selectAtoms) = @_;
	my ($i, $j, $resNum, $firstAtm);

	for $i (sort numerically keys %{ $selectAtoms }) {
		next if ($atoms->{$i}{UPDATED});
		if (! defined($resNum)) {
			$firstAtm = $i;
			if ($i == 1) {
				$resNum = 0;
			} else {
				$resNum = $atoms->{$i -1}{RESNUM};
			}
		}
		$resNum++;
		$atoms->{$i}{RESNUM} = $resNum;
		$atoms->{$i}{UPDATED} = 1;
		next if (! exists($atoms->{$i}{MOLECULE}));
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			$atoms->{$j}{RESNUM} = $resNum;
			$atoms->{$j}{UPDATED} = 1;
		}
	}

	for $i  (keys %{ $atoms }) {
		next if ($atoms->{$i}{UPDATED} or $i <= $firstAtm);
		$atoms->{$i}{RESNUM} += $resNum;
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

	$selection = BuildAtomSelectionString($atomSelect);

	print "Done\n";
}
