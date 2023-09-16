#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use MolData;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use Superimpose;
use General qw(LoadElements FindElement);

sub numerically { ($a<=>$b); }

my ($fileType, $fileName, $selStr, $saveFile);
my ($molStructure, $mutantRes, $atomList, $resList, $ELEMENTS);

$|++;
&init;
print "Gettting data from $fileType->{struct} file $fileName->{struct}...";
$molStructure->read($fileName->{struct}, $fileType->{struct});
print "Done\nSelecting atoms according to \"$selStr\"...";
$atomList = $molStructure->findAtoms($selStr);
die "ERROR: No valid atom found from selection\n" if (! $atomList);
$resList = getResList($atomList);
print "Done\nGetting mutant residue data from $fileType->{mutant} file $fileName->{mutant}...";
$mutantRes = setMutantData($fileName, $fileType, $resList);
print "Done\n";
&mutateResidues($molStructure, $mutantRes, $resList);
print "Creating $fileType->{save} file $fileName->{save}...";
$molStructure->write($fileName->{save}, $fileType->{save});
print "Done\n";

sub mutateResidues {
	my ($atoms, $mAtoms, $mList) = @_;
	my ($i, $mCount, $res, $atm, $resName, $tmp, $cAtoms); 
	my ($mutateName, $mutateStr, $printStr, $rBonds);
	
	$mCount = scalar(keys %{ $mList });

	for $i (sort numerically keys %{ $mList }) {
		$res = $atoms->shash->{resid}{$i}; #pointer to residue
		@{ $tmp } = keys %{ $res };
		$atm = $res->{ $tmp->[0] }; #first atom in residue
		$resName = $atm->resname; #residue name
		$mutateName = $mAtoms->{$i}->shash->{index}{1}{1}->resname if (!defined($mutateName));
		$printStr = "Mutating residue # $i ($resName) to $mutateName...";
		$mutateStr .= " # $i ($resName)";
		print "${printStr}\r";
		#first align mutant residue based on backbone atoms
		$cAtoms = alignRes($atoms, $i, $mAtoms->{$i});
		#delete old residue
		$atoms->deleteRes($i, \%{ $rBonds });
		#insert new mutant residue
		$atoms->insertRes($i, $cAtoms, $rBonds);
	}
	printf "%-" . (length($printStr) + 10) . "s\n", "Mutating residues: Sucessfully mutated $mutateStr -> $mutateName";
}

sub alignRes {
	my ($atoms, $resID, $mutantResAtoms) = @_;
	my ($i, $j, $d, $m, $tmp, $aList, $used, $mRes, $count);

	%{ $mRes } = ();
	@{ $tmp } = keys %{ $mutantResAtoms->shash->{resid} };
	$count=1;
	for $i (values %{ $mutantResAtoms->shash->{resid}{ $tmp->[0] } } ) { #loop over all atoms in the mutant residue
		undef $m;
		$m = findResAtom($i->atmname, $atoms, $resID, \%{ $used }) if($i->atmname !~ /^\s*H/i); #search for corresponding backbone atom in original residue

		#now we get the cordinates of the backbone atoms for the original residue and the mutant and align the mutant
		for $j ("x", "y", "z") {
			$d = uc $j . "COORD";
			$aList->{full}{ $i->index }{$d} = $i->$j; #this stores the atom coordinates in the old style hash for use with Superimpose
			if(defined($m)) { #backbone atoms that are common to ref and mutant residues
				$aList->{mut}{$count}{$d} = $i->$j; #mutant
				$aList->{ref}{$count}{$d} = $m->$j; #reference		
			}
		}
		#store fftype (needed by CoM)
		$aList->{full}{ $i->index }{FFTYPE} = $i->fftype;

		if(defined($m)) {
			$aList->{mut}{$count}{FFTYPE} = $i->fftype; #mutant
			$aList->{ref}{$count}{FFTYPE} = $m->fftype; #reference
			$count++;
		}
	}
	die "ERROR: Need at list 5 atoms for alignment!\n" if($count < 5); #need at least 5 atoms/3 intersecting planes

	&Kearsley::SuperimposeAtoms($aList->{ref}, $aList->{mut}, $aList->{full}); #this aligns the mutant atoms
	#&Quarternion::SuperimposeAtoms($aList->{ref}, $aList->{mut}, $aList->{full}); #this aligns the mutant atoms
	#&Kabash::SuperimposeAtoms($aList->{ref}, $aList->{mut}, $aList->{full}); #this aligns the mutant atoms

	#now update the mutant residue atom positions
	for $i (keys %{ $aList->{full} }) { #aList holds the updated mutant atom positions
		$mRes->{$i} = $mutantResAtoms->shash->{index}{$i}{$i}; #mRes holds the mutant atoms (original positions)
		$mRes->{$i}->("resid", $resID); #update resid
		for $j ("x", "y", "z") {
			$d = uc $j . "COORD";
			$mRes->{$i}->($j, $aList->{full}{$i}{$d}); #update atom position of mutant residue atom from aList
		}
	}

	return $mRes;
}


sub findResAtom {
	my ($atmName, $atoms, $resID, $used) = @_;
	my ($mList, $sStr, $i, $j, $m);

	$sStr = "'atmname' eq '$atmName' and 'resid' == $resID";
	$mList = $atoms->findAtoms($sStr);
	for $i (values %{ $mList }) {
		if(! exists($used->{ $i->index })) {
			$m = $i;
			$used->{ $i->index } = 1;
			last;
		}
	}

	return $m;
}

sub getResList {
	my ($atoms) = $_[0];
	my (%RES, $i);

	for $i (keys %{ $atoms }) {
		$RES{ $atoms->{$i}->resid } = 1
	}

	return \%RES;
}

sub setMutantData {
	my ($fle, $fTy, $res) = @_;
	my ($tmpMut, $nmutRes, $i, $mut);

	$tmpMut =  MolData->new();
	$tmpMut->read($fle->{mutant}, $fTy->{mutant});
	$nmutRes = scalar(keys %{ $tmpMut->shash->{resid}});
	die "ERROR: Found $nmutRes residues in mutant residue file! Expected 1... Aborting\n" if($nmutRes > 1);
	#now we read in the mutant residue file multiple times in order to not have overwriting issues

	for $i (keys %{ $res }) {
		$mut->{$i} = MolData->new();
		$mut->{$i}->read($fle->{mutant}, $fTy->{mutant});
	}

	return $mut;
}

sub init {
	my (%OPTS, $i);
	getopt('bsrm',\%OPTS);
	
	$molStructure = MolData->new();

	for $i ("b", "r") {
		die "usage: $0 -b bgf file -r residue bgf file -m (mutation residue(s) default all) -s (savename)\n"
			if (! exists($OPTS{$i}));
	}
	print "Initializing...";
	for $i ("b", "r") {
		$molStructure->testFile($OPTS{$i});
	}
	($fileName->{struct},$fileName->{save},$fileName->{mutant}, $selStr) = ($OPTS{b},$OPTS{s},$OPTS{r},$OPTS{m});
	$selStr = "'index'>0" if (! defined($selStr));

	#get file types
	for $i ("struct", "mutant", "save") {
		$fileType->{$i} = getFileType($fileName->{$i});
	}

	$ELEMENTS = &LoadElements;
	print "Done\n";
}

sub getFileType {
	my ($fName) = $_[0];
	my ($fType);

	$fType = "bgf";
	if ($fName =~ /\.(\w+)$/) {
		$fType = lc $1
	}
	return $fType;
}