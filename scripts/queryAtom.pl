#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(ParseStructFile);
use General qw(HasCell);
use Getopt::Std qw(getopt);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString AddMolsToSelection);
use BOX qw(GetBox);

sub usage;
sub testField;
sub printAtoms;
sub hasCell;
sub updateNumbonds;
sub numerically { ($a<=>$b) }

my ($bgfFile, $selection, $field, $sfield, $molOpt, $dList);
my ($ATOMS, $BONDS, $BOX, $HEADERS, $SELECTIONS, $MOLS, $ret);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($bgfFile,1);
&updateNumbonds($ATOMS, $BONDS);
$BOX = GetBox($ATOMS, undef, $HEADERS) if (HasCell($HEADERS));
&removeBonds($ATOMS, $BONDS, $dList, $BOX) if (defined($dList));
$MOLS = GetMols($ATOMS, $BONDS);
$ret = testField($ATOMS->{1}, $field);
die "ERROR: field \"$field\" is invalid\n" 
	if (! $ret);
if (defined($sfield)) {
	$ret = testField($ATOMS->{1}, $sfield);
	undef $sfield
		if (! $ret);
}
print "Done\nSelecting relevant atoms...";
$SELECTIONS = SelectAtoms($selection, $ATOMS, $BOX);
die "ERROR: No atoms matched selection\n" if (! keys %{ $SELECTIONS });
&AddMolsToSelection($SELECTIONS, $ATOMS) if ($molOpt);
print "Done\n";
&printAtoms($ATOMS, $SELECTIONS, $field, $sfield);

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
	my ($i,$nbonds);

	for $i (keys %{ $bonds }) {
		$nbonds = 0;
		$nbonds = scalar(@{ $bonds->{$i} }) if (exists($bonds->{$i}) and ref($bonds->{$i}) eq "ARRAY");
		$atoms->{$i}{NUMBONDS} = $nbonds;
	}
}

sub printAtoms {
	my ($atoms, $slist, $field, $sfield) = @_;
	my ($i, $j, $val);

	for $i (sort numerically keys %{ $slist }) {
		$val = "";
		for $j ($field,$sfield) {
			next if (! defined($j));
			if ($j =~ /MOLECULEID|MOLSIZE/) {
				$val .= "$j ${ $atoms->{$i}{$j} } ";
			} else {
				$val .= "$j $atoms->{$i}{$j} ";
			}
		}
		print "ATOM $i $val\n";
	}
}

sub testField {
	my ($atom, $fieldN) = @_;
	my ($i, $valid);

	for $i (keys %{ $atom }) {
		next if (!ref($atom->{$i}) eq "SCALAR");
		$valid = 1 if ($fieldN eq $i);
	}
	return $valid;
}

sub init {
	my (%OPTS, $atomSel, $delStr, $tmp);
	getopt('bafsmd',\%OPTS);
	($bgfFile, $atomSel, $field, $sfield, $molOpt, $delStr) = 
		($OPTS{b},$OPTS{a},$OPTS{f},$OPTS{s},$OPTS{m},$OPTS{d});
	for ($bgfFile, $atomSel,$field) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	$selection = BuildAtomSelectionString($atomSel);
	$molOpt = 0 if (! defined($molOpt) or $molOpt =~ /0|no/i);
	$molOpt = 1 if ($molOpt =~ /1|yes/i);
	die "ERROR: Expected word for field, Got \"$field\"\n" if ($field !~ /^(\w+)/);
	$field = uc $1;
	if(defined($sfield)) {
		if($sfield =~ /^(\w+)/) {
			$sfield = uc $1;
		} else {
			undef $sfield;
		}
	}
	if(defined($delStr) and $delStr =~ /\:\:/) {
		@{ $tmp } = split /::/,$delStr;
		$dList->{1} = BuildAtomSelectionString($tmp->[0]);
		$dList->{2} = BuildAtomSelectionString($tmp->[1]);
	}
	print "Done\n";
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b struct_file -a atom_selection -f field -s (secondary_field) -m (mol_opt) -d "(delete_bonds)"
Arguments:
  struct_file: name of struct_file
  field: Field to display
  secondary field: Optional. Print this field as well
  mol_opt: Optional. Add the remaning atoms in the molecule to the selection. 
           Prevents breaking up of molecules. Default is no
  delete_bonds: Optional. Temporarily delete bonds between atom. Expected 'atom1_selection::atom2_selection'		   
  atom selection options:
    any valid bgf field expression. E.g. resname eq 'WAT' will select
    all the "WAT" residues while index > 10 will select all indices > 10.
    combine multiple expressions to make complicated selections: e.g.
    (xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
