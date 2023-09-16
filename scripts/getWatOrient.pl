#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);

sub usage;

my ($bgfFile, $selection, $dim);
my ($ATOMS, $BONDS, $SELECTIONS, );

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
print "Done\nCalculating OH bond vector...";
$DATA = calcBondVector($ATOMS, $SELECTIONS, $dim);
print "Done\n";

sub calcBondVector {
}

sub init {
	my (%OPTS, $oSel);
	getopt('bda',\%OPTS);
	($bgfFile, $oSel, $dim) = ($OPTS{b},$OPTS{a},$OPTS{d});
	&usage if (! defined($bgfFile));
	print "Initializing...";
	FileTester($bgfFile);
	$oSel = "fftype eq 'OW'" if (! defined($oSel));
	$selection = BuildAtomSelectionString($oSel);
	$dim = "z" if (! defined($dim) or $dim !~ /x|y|z/i);
	$dim =~ /(x|y|z)/;
	$dim = uc $1;
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -d (axis_dimension) -a (oxygen_selection)
Arguments:
  bgf_file: name of bgf_file
  axis_dimension: plane to calculate OH vector to. x|y|z. Default z
  oxygen_selection:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
	Default is "fftype eq 'OW'"
DATA
die "\n";

}
