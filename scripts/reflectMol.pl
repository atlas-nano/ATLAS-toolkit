#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);

sub usage;
sub reflectAtoms;

my ($bgfFile, $saveName, $selection, $molOpt, $line);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
print "Done\nReflecting " . scalar(keys %{ $SELECTIONS }) . " atoms about the line $line->{axis} = $line->{val}...";
&reflectAtoms($ATOMS, $SELECTIONS, $line);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub reflectAtoms {
    my ($atoms, $select, $lineOpts) = @_;
    my ($i, $j, $dim, $val);

    $dim = uc $lineOpts->{axis} . "COORD";
    for $i (keys %{ $select } ) {
	$atoms->{$i}{$dim} = $lineOpts->{val}-$atoms->{$i}{$dim};
    }
}

sub init {
    my (%OPTS, $atomSel,$lineEqn);
    getopt('bloms',\%OPTS);
    ($bgfFile, $lineEqn, $saveName, $atomSel, $molOpt) = ($OPTS{b},$OPTS{l},$OPTS{s},$OPTS{o}, $OPTS{m});
    for ($bgfFile, $lineEqn) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    if ($lineEqn =~ /(x|y|z) = (\-?\d+\.?\d*)/i) {
	$line->{axis} = lc $1;
	$line->{val} = $2;
    } else {
	die "ERROR: Expected x|y|z = xx.xx for line. Got \"$lineEqn\"!\n";
    }
    $atomSel = "index>0" if !(defined($atomSel));
    $selection = BuildAtomSelectionString($atomSel);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$//;
	$saveName .= ".reflect.bgf";
    }
    $molOpt = 0 if (! defined($molOpt) or $molOpt !~ /yes|1/i);
    $molOpt = 1 if ($molOpt =~ /yes|1/i);
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -l line -o (atom selection) -m (entire mol = no)
Arguments:
  bgf_file: name of bgf_file
  line: x|y|z = xx.xx
  save_name: name of file to save
  entire mol: select all atoms belonging to selected molecule even if not selected.
	default = no
  atom selection:
    any valid bgf field expression. E.g. resname eq 'WAT' will select
    all the "WAT" residues while index > 10 will select all indices > 10.
    combine multiple expressions to make complicated selections: e.g.
    (xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
