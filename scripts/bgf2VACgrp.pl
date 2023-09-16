#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
sub createGrpFile;
sub getGrps;
sub numerically { ($a<=>$b); }
sub getConstraints;
sub getRange;

my ($field, $bgfFile, $saveName, $shake, $opts, $selection, $sortField);
my ($ATOMS, $BONDS, $GRPS, $zeroConstraints,$SELECTION);

$|++;
&init;

print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSelecting Relevant Atoms...";
&GetMols($ATOMS, $BONDS);
$SELECTION = SelectAtoms($selection, $ATOMS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTION });
&AddMolsToSelection($SELECTION, $ATOMS) if ($field =~ /MOL/i);
print "Done\nSelecting atoms by $field...";
$GRPS = getGrps($ATOMS, $field, $SELECTION, $sortField);
print "Done\nCreating Group file ${saveName}...";
&createGrpFile($GRPS, $saveName, $shake, $opts);
print "Done\n";

sub getGrps {
	my ($atoms, $field, $select, $sortKey) = @_;
	my ($gHASH, $i, $hKey, $tmp, $gMap, @list);

	for $i (keys %{ $select }) {
		$tmp->{ $atoms->{$i}{$field} } = 1 if ($field !~ /MOL/i);
		$tmp->{ ${ $atoms->{$i}{$field} } } = 1 if ($field =~ /MOL/i);
	}
	if ($field =~ /CHAIN|FFTYPE|RESNAME/) {
		@list = sort {$a cmp $b} keys %{ $tmp };
	} else {
		@list = sort numerically keys %{ $tmp };
	}
	for $i (0 .. $#list) {
		$gMap->{$list[$i]} = $i;
	}
	for $i (keys %{ $select }) {
		$hKey = $gMap->{ $atoms->{$i}{$field} } if ($field !~ /MOL/i);
		$hKey = $gMap->{ ${ $atoms->{$i}{$field}  }} if ($field =~ /MOL/i);
		if (defined($sortKey)) {
			$gHASH->{ $atoms->{$i}{$sortKey} }{$hKey}{$i} = $atoms->{$i}{$sortKey};
		} else {
			$gHASH->{1}{$hKey}{$i} = 1;
		}

	}

	return $gHASH;
}

sub createGrpFile {
	my ($data, $saveFile, $shakeOpt, $otherOpts) = @_;
	my ($i, @tmp, @tmp1, $j, $k, $constraints, $grpStr, $count);

	for $i (keys %{ $data }) {
		$count += scalar(keys %{ $data->{$i} });
	}
	open OUTDATA, "> $saveFile" || die "ERROR: Cannot write to ${saveFile}: $!\n";
	print OUTDATA "Total Groups: $count\n";

	$count = 1;
	@tmp1 = sort { $a cmp $b } keys %{ $data };
	for $k (@tmp1) {
		@tmp = sort numerically keys %{ $data->{$k} };
		for $i (@tmp) {
			$j = scalar(keys %{ $data->{$k}{$i} });
			print OUTDATA "Group $count Atoms $j\n";
			$grpStr = getRange($data->{$k}{$i});
			print OUTDATA "$grpStr\n";
			$count++;
		}
	}
	if ($shakeOpt) { 
			$constraints = getConstraints($data, $otherOpts, $zeroConstraints);
	} else { 
			$constraints = "";
	}
	print OUTDATA $constraints;
	close OUTDATA;
}

sub getRange {
	my ($grpAtoms) = @_;
	my ($sorted_list, $i, $j, $retStr, $range);

	@{ $sorted_list } = sort numerically keys %{ $grpAtoms };
	return $sorted_list->[0] if ($#{ $sorted_list } == 0);
	$i = 0;
	$j = 1;
	while (($i+$j) <= $#{ $sorted_list }) {
		if (($sorted_list->[$i+$j] - $sorted_list->[$i]) <= $j) {
			$j++;
			next;
		}
		if (($sorted_list->[$i]-$sorted_list->[$i+$j-1]) == 0) {
			push @{ $range }, "$sorted_list->[$i] ";
		} else {
			push @{ $range }, "$sorted_list->[$i] - " . $sorted_list->[$i+$j-1]. " ";
		}
		$i += $j;
		$j = 1;
	}
	if (($i-$j) == 1) {
		push @{ $range }, "$sorted_list->[$i]";
	} else {
		push @{ $range }, "$sorted_list->[$i] - " . $sorted_list->[$i+$j-1];
	}
	foreach $i (0 .. $#{ $range }) {
		$retStr .= $range->[$i];
		$retStr .= "\n" if ((($i+1) % 10) == 0);
	}
	return $retStr;
}

sub getConstraints {
	my ($data, $oOpts, $setZero) = @_;
	my ($i, $j, $k, $mols, @tmp, @tmp1, $cStr, $rStr, $lStr, $count, $eType, $atom, $tot);

	$cStr = "Constraints\n";
	$rStr = "RotationalSymmetryNumber\n";
	$lStr = "LinearMoleculeFlag\n";
	@tmp1 = sort { $a cmp $b } keys %{ $data };
	for $k (@tmp1) {
		@tmp = sort numerically keys %{ $data->{$k} };
		for $i (@tmp) {
			$count = 0;
			$tot = 0;
			$mols = ();
			$mols = GetMols($ATOMS, $BONDS, $data->{$k}{$i});
			for $j (keys %{ $mols }) {
				for $atom (keys %{ $mols->{$j}{MEMBERS} }) {
					$count++ if ($ATOMS->{$atom}{FFTYPE} =~ /^H/i);
				}
				$tot += $mols->{$j}{MOLSIZE};
				$count++ if ($mols->{$j}{MOLSIZE} == 3);
			}
			$cStr .= "$count " if (! $setZero);
			$cStr .= "0 " if ($setZero);
			if ((scalar(keys %{ $mols }) > 0) and (($tot/scalar(keys %{ $mols })) == 3)) {
				$rStr .= "2 "; 
			} else {
				$rStr .= "1 ";
			}
			$lStr .= "0 ";
		}
	}
	$rStr = $lStr = "" if(! $oOpts);
	return "$cStr\n$rStr\n$lStr\n";
}

sub init {
	my (%OPTS,$atomSelect);

	getopt('fbscozag',\%OPTS);

	($field, $bgfFile, $saveName, $shake, $opts, $zeroConstraints, $atomSelect, $sortField) = 
		($OPTS{f},$OPTS{b},$OPTS{s}, $OPTS{c}, $OPTS{o}, $OPTS{z}, $OPTS{a}, $OPTS{g});

	die "usage:$0 -b bgf file -f (field=chain|resid|resname|mol|fftype) -g (sort-field) -c (has shake=no) -o (other opts=yes if shake) -z (zero constraints) -s (savename) -a [atoms = all]\n" 
		if (! defined($OPTS{b}));
	print "Initializing...";

	FileTester($bgfFile);
	$atomSelect = "index>0" if (!defined($atomSelect));
	$selection = BuildAtomSelectionString($atomSelect);
	$field = "CHAIN" if (! defined($field));
	if ($field =~ /(chain|resname|res|molsize|mol|fftype)/i) {
		$field = uc $1;
		$field .= "ECULEID" if ($field eq "MOL");
		$field .= "NUM" if ($field eq "RES");
	}
	if (defined($sortField) and $sortField =~ /(chain|resname|res|molsize|mol|fftype|index)/i) {
		$sortField = uc $1;
		$sortField .= "ECULEID" if ($sortField eq "MOL");
		$sortField .= "NUM" if ($sortField eq "RES");
	} else {
		undef $sortField;
	}

	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$//;
		$saveName .= ".grps";
	}
	$shake = 0 if (! defined($shake) or $shake !~ /1|yes/i);
	$shake = 1 if (defined($shake) and $shake =~ /1|yes/i);
	$opts = 1 if (! defined($opts) or $opts !~ /0|no/i);
	$opts = 0 if (defined($opts) and $opts =~ /0|no/i);
	$zeroConstraints = 0 if (! defined($zeroConstraints) or $zeroConstraints !~ /1|yes/i);
	$zeroConstraints = 1 if ($zeroConstraints =~ /1|yes/i);
	print "Done\n";
}
