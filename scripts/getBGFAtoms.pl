#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(ParseStructFile addHeader createBGF GetBGFAtoms insertHeaderRemark);
use General qw(FileTester HasCell GetFileTypeStr);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use BOX qw(GetBox);

sub usage;

my ($structFile, $saveName, $selection, $molOpt, $atomSel);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF, $CONS, $tmp, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $structFile...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($structFile,1);
&GetMols($ATOMS, $BONDS);
$BOX = &GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\nCreating BGF file $saveName...";
&insertHeaderRemark($HEADERS, "REMARK selected atoms according to $atomSel");
&addHeader($BGF,$HEADERS);
&createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub init {
	my (%OPTS);
	getopt('boms',\%OPTS);
	($structFile, $saveName, $atomSel, $molOpt) = ($OPTS{b},$OPTS{s},$OPTS{o}, $OPTS{m});
	for ($structFile, $atomSel) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$structFile =~ /^\s*(\S+)/;
		$saveName = basename($1);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	$molOpt = 0 if (! defined($molOpt) or $molOpt !~ /yes|1/i);
	$molOpt = 1 if ($molOpt =~ /yes|1/i);
}

sub usage {
	my ($fTypeStr) = &GetFileTypeStr;
	print STDOUT <<DATA;
usage: $0 -b struct_file -s save_name -o options -m (entire mol = no)
Arguments:
  struct_file: name of structure file
$fTypeStr  
  save_name: name of file to save
  entire mol: select all atoms belonging to selected molecule even if not selected.
		default = no
  options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
