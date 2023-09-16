#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);

sub usage;

my ($bgfFile, $saveName, $selection, $molOpt);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF, $CONS, $tmp, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\nReimaging...";

print "Done\nCreating BGF file $saveName...";
&addHeader($BGF,$HEADERS);
&createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub init {
    my (%OPTS, $atomSel);
    getopt('boms',\%OPTS);
    ($bgfFile, $saveName, $atomSel, $molOpt) = ($OPTS{b},$OPTS{s},$OPTS{o}, $OPTS{m});
    for ($bgfFile, $atomSel) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    $selection = BuildAtomSelectionString($atomSel);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
    $molOpt = 0 if (! defined($molOpt) or $molOpt !~ /yes|1/i);
    $molOpt = 1 if ($molOpt =~ /yes|1/i);
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o options -m (entire mol = no)
Arguments:
  bgf_file: name of bgf_file
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
