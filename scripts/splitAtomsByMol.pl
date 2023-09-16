#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;

use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(ParseStructFile addHeader createBGF);
use General qw(GetFileTypeStr);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString GetMols);
use Storable qw(dclone);

sub init;
sub splitAtoms;
sub showUsage;
sub removeBondTypes;

my ($struct_file, $saveFile, $SELECT, $atomSelect, $field, $dList);
my ($ATOMS, $BONDS, $HEADERS, $MOLS, $CONS, $offset);

$|++;
&init;
print "Parsing bgf file $struct_file...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($struct_file, 1);
$CONS = dclone($BONDS);
$SELECT = SelectAtoms($atomSelect, $ATOMS);
die "No selected atoms found in bgf file\n" if (! keys %{ $SELECT });
print "Done\nDetermining molecules based on connectivity...";
&removeBondTypes($ATOMS, $CONS, $dList) if (defined($dList));
$MOLS = &GetMols($ATOMS,$CONS,$SELECT);
print "Done\nSpliting atoms in molecules by ${field}...";
&splitAtoms($ATOMS, $MOLS, $SELECT, $offset);
print "Done\nCreating bgf file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub removeBondTypes {
	my ($atoms, $bonds, $list) = @_;
	my ($i, $j, $itype, $jtype);

	for $i (keys %{ $atoms }) {
		$itype = $atoms->{$i}{FFTYPE};
		for $j (0 .. $#{ $bonds->{$i} }) {
			$jtype = $atoms->{ $bonds->{$i}[$j] }{FFTYPE};
			next if !( exists($list->{$itype}{$jtype}));
			splice @{ $bonds->{$i} }, $j, 1;
		}
	}
}

sub splitAtoms {
	my ($atoms, $mols, $selection, $offset) = @_;
	my ($i);

	for $i (keys %{ $selection }) {
		if ($field eq "CHAIN") {
			$atoms->{$i}{$field} = chr((64 + ${ $atoms->{$i}{MOLECULEID} } + $offset));
			$atoms->{$i}{$field} = "X" if (${ $atoms->{$i}{MOLECULEID} } > 4);
		} else {
			$atoms->{$i}{$field} = ${ $atoms->{$i}{MOLECULEID} } + $offset;
		}
	}
}

sub init {
	my (%OPTS, $select, $usage, $dBondStr);

	getopt('bafsdo',\%OPTS);
	$usage = showUsage();
	die "$usage\n" if (! defined($OPTS{b}));
	print "Initializing...";
	($struct_file, $saveFile, $select, $field, $dBondStr, $offset) = 
		($OPTS{b}, $OPTS{s}, $OPTS{a}, $OPTS{f}, $OPTS{d}, $OPTS{o});
	if (! defined($saveFile)) {
		$struct_file =~ /^(\S+)/;
		$saveFile = basename($1);
		$saveFile =~ s/\.\w+//;
		$saveFile .= "_mod.bgf";
	}
	
	$select = "index > 0" if (! defined($select));
	$atomSelect = BuildAtomSelectionString($select);
	$field = "CHAIN" if (! defined($field) or $field !~ /^CHAIN|RESNUM|RESID|RES$/i);
	$field = uc($field);
	$field = "RESNUM" if ($field =~ /RES/);
	if (defined($dBondStr)) {
		while ($dBondStr =~ /(\S+) \- (\S+)/g) {
			$dList->{$1}{$2} = 1;
			$dList->{$2}{$1} = 1;
		}
	}
	$offset = 0 if (! defined($offset) or $offset =~ /\D/i);
	print "Done\n";
}

sub showUsage {
	my ($fTypeStr) = &GetFileTypeStr;
	my ($usage) = <<DATA;
usage: $0 -b structure_file -a (atom_selection) -f (field) -s (save_name) -d "delete bond between type1 - type2 ..."
Options:
	-b structure file (Required): location of structure file
$fTypeStr	
	-a atom selection (optional): any valid field (fftype, charge, resname etc) expression
		e.g. -a "resname ne 'WAT' and resnum > 10". Default: all atoms
	-f field (optional): field to use to split molecules. Either chain or resid (default)
	-s save name: name of new bgf file
DATA
	return $usage;
}

