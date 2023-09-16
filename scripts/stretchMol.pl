#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester CoM);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);

sub usage;
sub stretchMols;

my ($bgfFile, $saveName, $selection, $molOpt, $stretchFactor);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
print "Done\nStretching molecules by $stretchFactor...";
&stretchMols($ATOMS, $SELECTIONS, $stretchFactor);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub stretchMols {
	my ($atoms, $sels, $sf) = @_;
	my ($i, $j, $k, $molCoM, $cMol);

	for $i (keys %{ $sels }) {
		next if (exists($atoms->{$i}{USED}));
		$cMol = ();
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			$cMol->{$j} = $atoms->{$j};
		}
		$atoms->{$i}{USED} = 1;
		$molCoM = CoM($cMol);
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			$atoms->{$j}{USED} = 1;
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$atoms->{$j}{$k} = $molCoM->{$k} + ($atoms->{$j}{$k} - $molCoM->{$k}) * $sf;
			}
		}
	}

}

sub init {
	my (%OPTS, $atomSel);
	getopt('bofs',\%OPTS);
	($bgfFile, $saveName, $atomSel, $stretchFactor) = ($OPTS{b},$OPTS{s},$OPTS{o}, $OPTS{f});
	for ($bgfFile, $atomSel, $stretchFactor) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	die "ERROR: Expected positive integer or decimal for stretch_factor! Got \'$stretchFactor\'"
		if ($stretchFactor !~ /^\d+\.?\d*$/);
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o options -f stretch_factor
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
