#!/usr/bin/perl -w

use constant PI => 3.14159265358979;
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester Rotate CoM);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString GetAtmData AddMolsToSelection);

sub usage;
sub rotateMols;
sub getRandomNumber;

my ($bgfFile, $saveName, $selection);
my ($ATOMS, $BONDS, $SELECTIONS, $HEADERS, $MOLS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $ATOMS });
print "Done\nRandomely rotating molecules...";
&rotateMols($ATOMS, $MOLS);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub rotateMols {
	my ($atoms, $mols) = @_;
	my ($i, $j, $k, $angles, $currMol, $molcom);

	for $i (keys %{ $mols }) {
		$currMol = GetAtmData($atoms, $mols->{$i}{MEMBERS});
		$molcom = CoM($currMol);
		for $j (keys %{ $currMol }) {
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$currMol->{$j}{$k} -= $molcom->{$k};
			}
		}
	    $angles = getRandomNumber(360);
		Rotate($currMol, $angles, 3);
		for $j (keys %{ $currMol }) {
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$currMol->{$j}{$k} += $molcom->{$k};
			}
		}
	}
}

sub getRandomNumber(@) {
    my ($max) = $_[0];
    my (@output_array, $rand_angle, $counter);

    for $counter (0 .. 2) {
        $rand_angle = (rand $max) * PI/180;
        $output_array[$counter] = sprintf("%.2f", $rand_angle);
    }

    return \@output_array;

}

sub init {
	my (%OPTS, $atomSel);
	getopt('bms',\%OPTS);
	($bgfFile, $saveName, $atomSel) = ($OPTS{b},$OPTS{s},$OPTS{o});
	&usage if (! defined($bgfFile));
	print "Initializing...";
	$atomSel = "index>0" if (!defined($atomSel));
	FileTester($bgfFile);
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o options
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
