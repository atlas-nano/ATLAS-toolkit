#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms insertHeaderRemark);
use General qw(FileTester HasCell GetBondLength);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(SelectAtoms);
use BOX qw(GetBox);

sub usage;

my ($bgfFile, $saveName, $selection, $tol);
my ($ATOMS, $BONDS, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$BOX = &GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
print "Done\nRemoving atoms within ${tol}A of other atoms...";
($ATOMS, $BONDS) = removeOverlappingAtoms($ATOMS, $BONDS, $BOX, $tol);
print "Done\nCreating BGF file $saveName...";
&insertHeaderRemark($HEADERS, "REMARK remove overlapping atoms within $tol");
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub removeOverlappingAtoms {
	my ($atoms, $bonds, $box, $tol) = @_;
	my ($i, $j, $count, $bl, $nAtoms, $nBonds, $nList);

	%{ $nList } = keys %{ $atoms };

	for $i (keys %{ $nList }) {
		next if (! exists($nList->{$i}));
		for $j (keys %{ $nList }) {
			next if (! exists($nList->{$j}) or $i == $j);
			$bl = GetBondLength($atoms->{$i}, $atoms->{$j}, $BOX,0);
			if ($bl <= $tol) { #overlapping atom
				delete $nList->{$j};
				$count++;
			}
		}
	}
	($nAtoms, $nBonds, undef) = GetBGFAtoms($nList, $atoms, $bonds);
	print "removed $count atoms...";
	return ($nAtoms, $nBonds);
}

sub init {
	my (%OPTS);
	getopt('bts',\%OPTS);
	($bgfFile, $saveName, $tol) = ($OPTS{b},$OPTS{s},$OPTS{t});
	&usage if (! defined($bgfFile));
	
	print "Initializing...";
	FileTester($bgfFile);
	$tol=0.5 if (! defined($tol) or $tol !~ /\d+\.?\d*/);
	$tol = $1 if ($tol =~ /(\d+\.?\d*)/);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -t (tolerance=0.5A)
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  tolerance: minimum distance between atoms to consider as nonoverlapping. 
		default = 0.5A
  options:
	any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
