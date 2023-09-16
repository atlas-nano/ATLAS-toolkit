#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Storable qw(dclone);
use FileFormats qw(ParseStructFile addHeader insertHeaderRemark
				createBGF createPDB createMOL2 GetMOL2Header);
use General qw(FileTester ShuffleArray HasCell GetFileTypeStr);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols AddMolsToSelection);
use Getopt::Std qw(getopt);
use BOX qw(GetBox);

my ($FILES, $FIELDS, $selection, $SELECT, $saveType, $fieldStr, $molOpt, $nrand, $dList);
my ($ATOMS, $BONDS, $CONS, $HEADERS, $BOX, $SELECTATMS, $saveFunc);

$|++;
&init;
print "Parsing Structure file $FILES->{STRUCTURE}...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($FILES->{STRUCTURE},1);
&GetMOL2Header($ATOMS, \@{ $HEADERS }) if (exists($ATOMS->{HEADER}));
$CONS = dclone($BONDS);
&removeBonds($ATOMS, $CONS, $dList, $BOX) if (defined($dList));
&GetMols($ATOMS, $CONS);
$BOX = &GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
&updateNumbonds($ATOMS, $BONDS);
if ($saveType =~ /mol2/) {
	$HEADERS = $ATOMS->{HEADER};
	delete $ATOMS->{HEADER};
}
&parseFieldList($ATOMS, $FIELDS);
print "Done\nParsing atom/residue selection...";
$SELECTATMS = SelectAtoms($selection, $ATOMS);
die "ERROR: Not a valid atom selection!\n" if (! keys %{ $SELECTATMS });
&selectRandom($SELECTATMS, $nrand) if ($nrand);
&AddMolsToSelection($SELECTATMS, $ATOMS) if ($molOpt);
print "Done\nUpdating fields..";
&updateAtomFields($ATOMS, $SELECTATMS, $FIELDS);
print "Done\nCreating $FILES->{SAVE}...";
&insertHeaderRemark($HEADERS, "REMARK $FILES->{STRUCTURE} modified $fieldStr") 
		if ($saveType eq "bgf");
&addHeader($ATOMS, $HEADERS) if ($saveType eq "bgf");
$saveFunc->($ATOMS, $BONDS, $FILES->{SAVE});
print "Done\n";

sub removeBonds {
	my ($atoms, $bonds, $dL, $box) = @_;
	my ($a1, $a2, $i, $j);

	$a1 = SelectAtoms($dL->{1}, $atoms, $box);
	die "ERROR: No atoms matched selection $dL->{1}\n"
		if(! keys %{ $a1 });
	$a2 = SelectAtoms($dL->{2}, $atoms, $box);
	die "ERROR: No atoms matched selection $dL->{2}\n"
		if(! keys %{ $a2 });
	for $i (keys %{ $a1 }) {
		$j = 0;
		while ($j <= $#{ $bonds->{$i} }) {
			if(exists($a2->{$bonds->{$i}[$j]})) {
				splice @{ $bonds->{$i} }, $j, 1; 
			} else {
				$j++;
			}
		}
	}
	for $i (keys %{ $a2 }) {
		$j = 0;
		while ($j <= $#{ $bonds->{$i} }) {
			if(exists($a1->{$bonds->{$i}[$j]})) {
				splice @{ $bonds->{$i} }, $j, 1; 
			} else {
				$j++;
			}
		}
	}
}

sub updateNumbonds {
	my ($atoms, $bonds) = @_;
	my ($i);

	for $i (keys %{ $bonds }) {
		$atoms->{$i}{NUMBONDS} = 0;
		$atoms->{$i}{NUMBONDS} = scalar(@{ $bonds->{$i} }) if ($bonds->{$i});
	}
}

sub selectRandom {
	my ($select, $num) = @_;
	my ($list, $tot, $i);

	$tot = scalar(keys %{ $select });
	return if ($tot < $num);
	@{ $list } = keys %{ $select };
	&ShuffleArray($list);
	for $i ($num .. $tot) {
		delete $select->{$list->[$i]};
	}
}

sub updateAtomFields {
	my ($atoms, $atomSelection, $fields) = @_;
	my ($i, $j, $tmp);

	for $i (keys %{ $atomSelection }) {
		for $j (keys %{ $fields }) {
			if ($fields->{$j}{MOD} eq ".") {
				$atoms->{$i}{$j} .= $fields->{$j}{VAL};
			} elsif ($fields->{$j}{MOD} eq "+") {
				$atoms->{$i}{$j} += $fields->{$j}{VAL};
			} elsif ($fields->{$j}{MOD} eq "-") {
				$atoms->{$i}{$j} -= $fields->{$j}{VAL};
			} else {
				$atoms->{$i}{$j} = $fields->{$j}{VAL};
			}
		}
	}
}

sub parseFieldList {
	my ($atoms, $fields) = @_;
	my ($i, $currAtm, $fieldList, $j);

	for $i (keys %{ $atoms }) {
		$currAtm = \%{ $atoms->{$i} };
		for $j (keys %{ $fields }) {
			delete $fields->{$j} if (! exists($currAtm->{$j}));
			$fieldList .= "$j ";
		}
		last;
	}

	die "ERROR: No valid fields found while search for ${fieldList} in bgf file\n"
		if (! keys %{ $fields });
}

sub init {
	my (%OPTS, $usage, $atomSelect, $tmp, $rec, $delStr);

	getopt('safwtmrd',\%OPTS);

	$usage = &showUsage;
	for ("s","a","f") {
		die "$usage\n" if (! defined($OPTS{$_}));
	}
	print "Initializing...";
	($FILES->{STRUCTURE}, $atomSelect, $fieldStr, $FILES->{SAVE}, $saveType, $molOpt,  $nrand,   $delStr) = 
	($OPTS{s},            $OPTS{a},    $OPTS{f},  $OPTS{w},       $OPTS{t},  $OPTS{m}, $OPTS{r}, $OPTS{d});

	$atomSelect = "*" if (! defined($atomSelect));
	$selection = BuildAtomSelectionString($atomSelect);

	while ($fieldStr =~ /(\S+)/g) {
		$tmp = $1;
		if ($tmp =~ /(\w+):(\+|\-|\.)?(.*)/) {
			if (defined($2)) {
				$rec = (
						{
							"MOD" => $2,
							"VAL" => $3,
						}
						);
			} else {
				$rec = (
						{
							"MOD" => "",
							"VAL" => $3,
						}
						);
			}
			$FIELDS->{$1} = $rec;
		}
	}

	$saveType = "bgf" if (! defined($saveType));
	$saveType = lc($saveType);
	$saveFunc = \&createBGF;
	if ($saveType =~ /mol2/) {
		$saveFunc = \&createMOL2;
	} elsif ($saveType =~ /pdb/) {
		$saveFunc = \&createPDB;
	}
	if (! defined($FILES->{SAVE})) {
		$FILES->{STRUCTURE} =~ /^\s*(\S+)/;
		$FILES->{SAVE} = $1;
		$FILES->{SAVE} =~ s/\.\w+$//;
		$FILES->{SAVE} .= "_mod.${saveType}";
	}

	$molOpt = 0 if (! defined($molOpt) or $molOpt !~ /1|yes/i);
	$molOpt = 1 if ($molOpt =~ /1|yes/i);

	$nrand = 0 if (! defined($nrand) or $nrand !~ /^\d+$/);
	$nrand = $1 if ($nrand =~ /(\d+)/);

	die "ERROR: Invalid string found while searching \"FIELD\".\n" .
		"Expected field:(+|-|.)new_val. Got \"${fieldStr}\"!\n" if (! defined($FIELDS));

	if(defined($delStr) and $delStr =~ /\:\:/) {
		$tmp = ();
		@{ $tmp } = split /::/,$delStr;
		$dList->{1} = BuildAtomSelectionString($tmp->[0]);
		$dList->{2} = BuildAtomSelectionString($tmp->[1]);
	}
	print "Done\n";
}

sub showUsage {
	my ($fTypeStr) = GetFileTypeStr;
	my $rStr = <<DATA;
usage: $0 -s structure file -a atom selection -f field -w [save name] -m [mol option] -r [random option] -t [save type] -d [delete_bonds]
Options:
	-s structure file: Location of structure file
$fTypeStr	
	-a atom selection:
		any valid bgf field expression. E.g. resname eq 'WAT' will select
		all the "WAT" residues while index > 10 will select all indices > 10.
		combine multiple expressions to make complicated selections: e.g.
		(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
	-f field: field to adjust. Expected field:[+|-|.]new_val. Enclose multiple in quotes
	-d [delete_bonds]: (Optional) Temporarily delete bonds between atom. 
		Expected 'atom1_selection::atom2_selection'. Default in none
	-w [save name]: (Optional) Name to save the file as
	-m [mol option]: (Optional) Add atoms of selected molecule to selection. Default 0
	-r [random]: (Optional) only change these many randomely selected atoms/mols
	-t [save type]: (Optional) File formats (bgf (default), mol2 or pdb) to save file as
DATA

	return $rStr;
}
