#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use MolData;
use General qw(PrintProgress GetTime ShuffleArray);

my ($fileType, $fileName, $num, $selStr, $molOpt, $saveName, $randomize);
my ($molStructure, $delList);

$|++;
&init;
print "Gettting data from $fileType file $fileName...";
$molStructure->read($fileName, $fileType);
print "Done\nSelecting atoms according to \"$selStr\"...";
$delList = $molStructure->findAtoms($selStr, $num, $molOpt, $randomize);
die "ERROR: No valid atom found from selection\n" if (! $delList);
print "Done\n";
&remove($molStructure, $delList, $molOpt);
print "Creating $fileType file $saveName...";
$molStructure->write($saveName, $fileType);
print "Done\n";

sub remove {
	my ($structure, $list, $removeMol) = @_;
	my ($i, $tot, $pStr, $start, $count, $strLen, $sList);

	$start = time();
	@{ $sList } = keys %{ $list };
	$tot = scalar(@{ $sList });
	&ShuffleArray($sList);
	$pStr = "Removing $tot atoms...";
	$pStr = "Removing $tot mols..." if($removeMol);
	print "${pStr}Calculating time remaining...\r";
	for $i (@{ $sList }) {
		if (! $removeMol) {
			$structure->deleteAtom($list->{$i});
		} else {
			$structure->deleteMol($list->{$i});
		}	
		$count++;
		$strLen = PrintProgress($count, $tot, $start, $pStr);
	}
	$tot = GetTime(time() - $start);
	printf "${pStr}%-${strLen}s\n", "${tot}s elapsed...Done";	
}

sub init {
	my (%OPTS);
	
	getopt('fsatmnsr', \%OPTS);
	for ("f", "a") {
		die "usage: $0 -f structure file -a atom selection \n\t-t (file type = bgf/mol2/pdb) " . 
				"-m (entire molecule=no) -n (number of atoms/mols) -r (randomize=0) -s (save name)\n"
			if (! exists($OPTS{$_}));
	}

	print "Initializing...";
	$molStructure =  MolData->new();
	($fileName, $fileType, $saveName, $molOpt, $num, $selStr, $randomize) = 
		($OPTS{f}, $OPTS{t}, $OPTS{s}, $OPTS{m}, $OPTS{n}, $OPTS{a}, $OPTS{r});
	$molStructure->testFile($fileName);
	if (! defined($fileType)) {
		$fileType = "bgf";
		if ($fileName =~ /\.(\w+)$/) {
			$fileType = lc $1;
		}
	}
	$saveName = $molStructure->getFileName($fileName) if (! defined($saveName));
	$molOpt = 0 if (! defined($molOpt) || $molOpt !~ /^1|yes$/i);
	$molOpt = 1 if ($molOpt =~ /^1|yes$/i);
	$randomize = 0 if (! defined($randomize) || $randomize !~ /1|yes/i);
	$randomize = 1 if ($randomize =~ /1|yes/i);
	$num = -1 if (! defined($num) || $num !~ /^\d+/);
	$num = -1 if (! $num);
	print "Done\n";
}
