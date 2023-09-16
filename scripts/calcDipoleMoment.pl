#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(ParseStructFile GetBGFAtoms);
use General qw(FileTester GetStats CoM GetFileTypeStr);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use Bulk qw(CalcDipole);
use Storable qw(dclone);

sub usage;
sub getAvgDipoleMoment;
sub printStats;

my ($structFile, $saveName, $selection, $molOpt);
my ($ATOMS, $BONDS, $SELECTIONS, $MOLS, $dM, $i);

$|++;
&init;
print "Getting atom information from $structFile...";
($ATOMS, $BONDS) = ParseStructFile($structFile,0);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS);
($ATOMS, $BONDS, undef) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
$MOLS = GetMols($ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $BONDS });
print "Done\n";
$dM = getAvgDipoleMoment($ATOMS, $MOLS);
&printStats($dM);

sub printStats {
	my ($dM) = $_[0];

	printf "====================================\nAvg Dipole Moment\n====================================\n";
	for $i ("T", "X", "Y", "Z") {
		printf "%-10s %10.5f +/- %10.5f\n",$i, $dM->{$i}{AVG}, $dM->{$i}{STDEV};
	}
}

sub getAvgDipoleMoment {
	my ($atoms, $mols) = @_;
	my ($i, $j, $k, $tmp, $currMol, $data, $dM, $com);


	for $i (keys %{ $mols }) {
		$currMol = ();
		for $j (keys %{ $mols->{$i}{MEMBERS} }) {
			$currMol->{$j} = dclone($atoms->{$j});
		}
		$com = CoM($currMol);
		for $j (keys %{ $currMol }) {
			for $k ("XCOORD", "YCOORD", "ZCOORD") {
				$currMol->{$j}{$k} -= $com->{$k};
			}
		}
		$tmp->{1} = $mols->{$i};
		$k = CalcDipole($currMol, $tmp, undef, 1);
		for $j ("X", "Y", "Z", "T") {
			push @{ $data->{$j} }, $k->{$j};
		}
	}

	for $i ("X", "Y", "Z", "T") {
		$dM->{$i} = GetStats($data->{$i});
	}

	return $dM;
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bo',\%OPTS);
	($structFile, $atomSel) = ($OPTS{b},$OPTS{o});
	&usage if (! defined($structFile));

	print "Initializing...";
	$atomSel = "index>0" if (! defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);
}

sub usage {
	my ($fTypeStr) = GetFileTypeStr;
	print STDOUT <<DATA;
usage: $0 -b struct_file -o (options)
Arguments:
  struct_file: name of structure file
$fTypeStr  
  options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
