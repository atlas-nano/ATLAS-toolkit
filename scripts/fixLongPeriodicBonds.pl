#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);

sub usage;
sub fixLongBonds;

my ($bgfFile, $saveName, $sel1, $sel2);
my ($ATOMS, $BONDS, $SELECT1, $SELECT2, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = GetBox($ATOMS,undef,$HEADERS);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECT1 = SelectAtoms($sel1, $ATOMS);
$SELECT2 = SelectAtoms($sel2, $ATOMS);
print "Done\nFixing bonds...";
&fixLongBonds($ATOMS, $BONDS, $SELECT1, $SELECT2, $BOX);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub fixLongBonds {
	my ($atoms, $bonds, $atom1, $atom2, $box) = @_;
	my ($i, $j, $k, $del, $del2);

	for $k ("X", "Y", "Z") {
		$BOX->{"${k}COORD"} = $BOX->{$k};
	}
	for $i (keys %{ $atom1 }) {
		for $j (@{ $bonds->{$i} }) {
			next if (! exists($atom2->{$j}));
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$del = $atoms->{$i}{$k}-$atoms->{$j}{$k};
				$del += $BOX->{$k}{len} if (abs($del+$BOX->{$k}{len})<abs($del)); 
				$del -= $BOX->{$k}{len} if (abs($del-$BOX->{$k}{len})<abs($del)); 
				$atoms->{$j}{$k} = $atoms->{$i}{$k} - $del; 
			}
		}
	}

}

sub init {
	my (%OPTS, $atoms1, $atoms2);
	getopt('bsij',\%OPTS);
	($bgfFile, $saveName, $atoms1, $atoms2) = ($OPTS{b},$OPTS{s},$OPTS{i},$OPTS{j});
	for ($bgfFile, $atoms2, $atoms1) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	$sel1 = BuildAtomSelectionString($atoms1);
	$sel2 = BuildAtomSelectionString($atoms2);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s (save_name) -i atoms2 -j atoms2
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  atom selection options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
