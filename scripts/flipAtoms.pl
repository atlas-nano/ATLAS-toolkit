#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester CrossProduct Normalize);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);

sub usage;
sub flipAtoms;

my ($bgfFile, $saveName, $selection, $DIM);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
print "Flipping atoms...";
&flipAtoms($ATOMS, $BONDS, $SELECTIONS, $DIM);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub flipAtoms {
    my ($atoms, $bonds, $select, $dim) = @_;
    my ($i, $j, $ref, $parent, $offset);

    for $i (sort { $a<=>$b } keys %{ $select }) {
	%{ $ref } = ({XCOORD=>0,YCOORD=>0,ZCOORD=>0});
	if ($#{ $bonds->{$i} } > -1) {
	    $parent = $atoms->{ $bonds->{$i}[0] };
	    for $j (qw /XCOORD YCOORD ZCOORD/) {
		$ref->{$j} = $parent->{$j};
	    }
	}
	for $j (keys %{ $dim }) {
	    $offset = 2 * ($atoms->{$i}{$j} - $ref->{$j});
	    $atoms->{$i}{$j} -= $offset;
	}
    }
}

sub init {
    my (%OPTS, $atomSel, $dimStr);

    getopt('bads',\%OPTS);
    ($bgfFile, $saveName, $atomSel, $dimStr) = ($OPTS{b},$OPTS{s},$OPTS{a}, $OPTS{d});
    for ($bgfFile, $atomSel, $dimStr) {
	&usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    $selection = BuildAtomSelectionString($atomSel);
    while ($dimStr =~ /(X|Y|Z)/gi) {
	$DIM->{uc $1. "COORD"} = 1;
    }
    die "ERROR: Expected x|y|z for flip dimension. Got \"$dimStr\"!\n"
    if (! defined($DIM));
    if (! defined($saveName)) {
	$saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -a atom_selection -d dimension 
Arguments:
  bgf_file: name of bgf_file
  atom_selection:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
  dimension: any combination of x,y,z
  save_name: name of file to save
DATA
die "\n";

}
