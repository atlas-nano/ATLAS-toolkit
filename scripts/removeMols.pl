#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(FileTester ShuffleArray);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);

sub init;
sub removeMols;

my ($bgfFile, $saveName, $selection, $randOpt);
my ($ATOMS, $BONDS, $HEADERS, $SELECT, $atomCount, $molCount);
my ($BGF, $CONS, $tmp, $MOLS);

$|++;
&init;
print "Done\nGetting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$MOLS = GetMols($ATOMS, $BONDS);
print "Done\nGetting selected atoms...";
$SELECT = SelectAtoms($selection, $ATOMS);
die "ERROR: No valid atoms selected!\n" if (! keys %{ $SELECT });
print "Done\nRemoving " . scalar(keys %{ $SELECT }) . " atoms...";
&removeMols($ATOMS, $SELECT, $atomCount, $molCount, $randOpt);
($BGF, $CONS, $tmp) = GetBGFAtoms($ATOMS, $ATOMS, $BONDS);
print "Done\nCreating BGF file $saveName...";
&addHeader($BGF,$HEADERS);
&createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub removeMols {
	my ($atoms, $atomList, $totAtoms, $totMols, $randomize) = @_;
	my ($remainAtoms, $i, $aCount, $mCount, $valid, $j, $index, $list);

	$mCount = $aCount = 0;
	@{ $list } = reverse sort {($a<=>$b)} keys %{ $atomList };
	&ShuffleArray($list) if (defined($randomize));
	MAIN: for $i (@{ $list }) {
		$valid = 0;
		next if (! exists($atoms->{$i}));
		for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
			if (exists($atoms->{$j})) {
				delete $atoms->{$j};
				delete $BONDS->{$j};
				$aCount++;
				$valid = 1;
				last MAIN if (defined($totAtoms) and $aCount == $totAtoms);
			}
		}
		$mCount++ if ($valid);
		last if (defined($totMols) and $mCount == $totMols);
	}
	
	die "ERROR: No atoms remain after deletion!\n"
		if (! keys %{ $atoms });

	print "removed $aCount atoms ($mCount mols)...";
}

sub init {
	my (%OPTS, $atomSel);
	getopt('basnmr',\%OPTS);
	
	die "usage: $0 -b bgf file -a atom selection -n (num atoms) -m (num mols) -s (savename) -r (randomize=no)\n"
		if (! exists($OPTS{b}) or ! exists($OPTS{a}));
	print "Initializing...";
	($bgfFile, $atomSel, $saveName, $atomCount, $molCount, $randOpt) = ($OPTS{b}, $OPTS{a}, $OPTS{s}, $OPTS{n}, $OPTS{m}, $OPTS{r});
	FileTester($bgfFile);
	$selection = BuildAtomSelectionString($atomSel);
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$/_mod\.bgf/;
	}
	undef($atomCount) if (defined($atomCount) and $atomCount !~ /^\d+$/);
	undef($molCount) if (defined($molCount) and $molCount !~ /^\d+$/);
	undef($randOpt) if (defined($randOpt) and $randOpt !~ /^1|yes/i);
	if (defined($randOpt)) {
		print "will randomize...";
		$randOpt = 1;
	}

}

