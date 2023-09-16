#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester IsInteger);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);

sub usage;
sub checkAtoms;
sub createBond;

my ($bgfFile, $saveName, $atom1, $atom2);
my ($ATOMS, $BONDS, $HEADERS, $SEL1, $SEL2);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
$SEL1 = SelectAtoms($atom1, $ATOMS);
$SEL2 = SelectAtoms($atom2, $ATOMS);
die "ERROR: No valid atoms selected for atom1\n"
	if (! keys %{ $SEL1 });
die "ERROR: No valid atoms selected for atom2\n"
	if (! keys %{ $SEL2 });
print "Done\nDeleting Bonds...";
&deleteBond(\%{ $BONDS }, $SEL1, $SEL2);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub deleteBond {
	my ($bonds, $seli, $selj) = @_;
	my ($i, $j, $ba);

	foreach $i (keys %{ $seli }) {
		$j = 0;
		while ($j <= $#{ $bonds->{$i} }) {
			$ba = $bonds->{$i}[$j];
			if (! exists($selj->{$ba})) {
				$j++;
			} else {
				splice @{ $bonds->{$i} }, $j, 1;
			}
		}
	}
	foreach $i (keys %{ $selj }) {
		$j = 0;
		while ($j <= $#{ $bonds->{$i} }) {
			$ba = $bonds->{$i}[$j];
			if (! exists($seli->{$ba})) {
				$j++;
			} else {
				splice @{ $bonds->{$i} }, $j, 1;
			}
		}
	}
}

sub init {
	my (%OPTS, $select, $atomSel1, $atomSel2);
	getopt('bijs',\%OPTS);
	($bgfFile, $atomSel1, $atomSel2, $saveName) = ($OPTS{b},$OPTS{i},$OPTS{j},$OPTS{s});
	for ($bgfFile, $atomSel1, $atomSel2) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	$atom1 = BuildAtomSelectionString($atomSel1);
	$atom2 = BuildAtomSelectionString($atomSel2);
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -i atom1 -j atom2 -s (save_name)
DATA

die "\n";

}
