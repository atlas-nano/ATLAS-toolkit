#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo GetBGFAtoms);
use General qw(FileTester GetStats CoM);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use Bulk qw(CalcQuadrupole);
use Storable qw(dclone);

sub usage;
sub getAvgQuadrupoleMoment;
sub printStats;

my ($bgfFile, $saveName, $selection, $molOpt);
my ($ATOMS, $BONDS, $SELECTIONS, $MOLS, $dM, $i, @qEle);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile,0);
&GetMols($ATOMS, $BONDS);
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS);
&AddMolsToSelection($SELECTIONS, $ATOMS);
($ATOMS, $BONDS, undef) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
$MOLS = GetMols($ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $BONDS });
print "Done\n";
$dM = getAvgQuadrupoleMoment($ATOMS, $MOLS);
&printStats($dM);

sub printStats {
	my ($dM) = $_[0];

	printf "====================================\nQuadrupole Moment (x10^26 e.s.u)\n====================================\n";
	for $i (@qEle) {
		printf "%-10s %10.5f +/- %10.5f\n",$i, $dM->{$i}{AVG}, $dM->{$i}{STDEV};
	}
}

sub getAvgQuadrupoleMoment {
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
				$mols->{$i}{R2} += $currMol->{$j}{$k}**2;
			}
		}
		$tmp->{1} = $mols->{$i};
		$k = CalcQuadrupole($currMol, $tmp, undef, 1);
		for $j (@qEle) {
			push @{ $data->{$j} }, $k->{$j};
		}
	}

	for $i (@qEle) {
		$dM->{$i} = GetStats($data->{$i});
	}

	return $dM;
}

sub init {
	my (%OPTS, $atomSel);
	getopt('bo',\%OPTS);
	($bgfFile, $atomSel) = ($OPTS{b},$OPTS{o});
	&usage if (! defined($bgfFile));

	print "Initializing...";
	FileTester($bgfFile);
	$atomSel = "index>0" if (! defined($atomSel));
	$selection = BuildAtomSelectionString($atomSel);
	@qEle=("AVG", "ANISO", "Qxx", "Qyy", "Qzz","XX","XY","YY","XZ","YZ","ZZ");
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -o (options)
Arguments:
  bgf_file: name of bgf_file
  options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
