#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester ShuffleArray);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(SelectAtoms BuildAtomSelectionString);

sub usage;

my ($bgfFile, $saveName, $sel1, $sel2, $num);
my ($ATOMS, $BONDS, $HEADERS, $SELECT, $FIELDS, $fStr);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
$fStr = testFields($ATOMS, $FIELDS);
print "Done\nSelecting relevant atoms...";
$SELECT->{1} = SelectAtoms($sel1, $ATOMS);
die "ERROR: Selected more atoms ($num) than exists in atom selection (" . scalar(keys %{ $SELECT->{1}}) . ")\n"
	if(scalar(keys %{ $SELECT->{1} }) < $num);
$SELECT->{2} = SelectAtoms($sel2, $ATOMS);
die "ERROR: Selected more atoms ($num) than exists in atom selection (" . scalar(keys %{ $SELECT->{2}}) . ")\n"
	if(scalar(keys %{ $SELECT->{2} }) < $num);
print "Done\nSwapping the $fStr of $num atoms...";
&swapFields($ATOMS, $SELECT, $FIELDS, $num);
print "Done\nCreating BGF file $saveName...";
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub swapFields {
	my ($atoms, $select, $fields, $num) = @_;
	my ($alist1, $alist2, $i, $j, $tmp);

	@{ $alist1 } = keys %{ $select->{1} };
	@{ $alist2 } = keys %{ $select->{2} };
	&ShuffleArray($alist1);
	&ShuffleArray($alist2);

	for $i (1 .. $num) {
		for $j (keys %{ $fields }) {
			$tmp = $atoms->{ $alist1->[$i-1] }{$j};
			$atoms->{ $alist1->[$i-1] }{$j} = $atoms->{ $alist2->[$i-1] }{$j};
			$atoms->{ $alist2->[$i-1] }{$j} = $tmp;
		}
	}

}

sub testFields {
	my ($atoms, $fields) = @_;
	my ($i, $retStr);

	$retStr = "";
	for $i (keys %{ $fields }) {
		if(!defined($atoms->{1}{$i})) {
			delete $fields->{$i};
		} else {
			$retStr .= "$i ";
		}
	}
	die "ERROR: No valid fields found!\n"
		if ($retStr eq "");
	return $retStr;
}

sub init {
	my (%OPTS);
	getopt('bsijnf',\%OPTS);
	($bgfFile, $saveName, $sel1, $sel2, $num, $fStr) = ($OPTS{b},$OPTS{s},$OPTS{i},$OPTS{j},$OPTS{n},$OPTS{f});
	for ($bgfFile, $sel1, $sel2, $num, $fStr) {
		&usage if (! defined($_));
	}
	print "Initializing...";
	FileTester($bgfFile);
	$sel1 = BuildAtomSelectionString($OPTS{i});
	$sel2 = BuildAtomSelectionString($OPTS{j});
	die "ERROR: Expected integer when searching number. Got \"$num\"\n"
		if ($num !~ /^\d+$/);
	while ($fStr =~ /(\S+)/g) {
		$FIELDS->{uc($1)} = 1;
	}
	if (! defined($saveName)) {
		$saveName = basename($bgfFile);
		$saveName =~ s/\.\w+$//;
		$saveName .= ".rand.bgf";
	}
}

sub usage {
	print STDOUT <<DATA;
usage: $0 -b bgf_file -s (save_name) -i atomselect1 -j atomselect2 -n number -f field(s)
Arguments:
  bgf_file: Required. name of bgf_file
  save_name: (Optional). name of file to save
  number: Required. Number of atoms positions to randomize
  field: Required. Field(s) to swap
  atomeslect: Required. any valid bgf field expression. E.g. resname eq 'WAT' will select
	all the "WAT" residues while index > 10 will select all indices > 10.
	combine multiple expressions to make complicated selections: e.g.
	(xcoord > 20.4 and moleculeid < 4) or sqrt((xcoord-23)**2+ycoord**2)>43.2
DATA
die "\n";

}
