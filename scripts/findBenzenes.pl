#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);

sub usage;
sub findBenzenes;
sub findRing;

my ($bgfFile, $saveName, $selection, $fftype);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
&findBenzenes($ATOMS, $BONDS, $SELECTIONS);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub findBenzenes {
    my ($atoms, $bonds, $select) = @_;
    my ($i, $j, $blist, $count);

    for $i (keys %{ $select }) {
	$select->{$i} = 0;
    }

    $count = 0;
    for $i (keys %{ $select }) {
	next if ($select->{$i} > 0);
	$count++;
	$blist = ();
	&find_ring($i, \%{ $blist }, $atoms, $bonds, $select);
	for $j (keys %{ $blist }) {
	    $select->{$j} = $count if(exists($select->{$j}));
	    $atoms->{$j}{RESNUM} = 444 + $count;
        }
    }

    print "Found $count benzenes...";

}

sub find_ring {
    my ($curr, $blist, $atoms, $bonds, $select) = @_;
    my ($i);

    for $i (@{ $bonds->{$curr} }) {
	if(exists($select->{$i})) {
	    next if(exists($blist->{$i}));
	    $blist->{$i} = 1;
	    &find_ring($i, \%{ $blist }, $atoms, $bonds, $select);
	}elsif($atoms->{$i}{FFTYPE} =~ /^H/) {
	    $blist->{$i} = 1;
	}
    }
}

sub init {
    my (%OPTS, $atomSel);
    getopt('bfs',\%OPTS);
    ($bgfFile, $saveName, $atomSel, $fftype) = ($OPTS{b},$OPTS{s},$OPTS{f});
    &usage if(! defined($bgfFile));
    print "Initializing...";
    FileTester($bgfFile);
    $atomSel = "fftype eq 'C_R'" if (!defined($atomSel));
    $selection = BuildAtomSelectionString($atomSel);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -f (benzen C fftype)
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
