#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(ParseStructFile addHeader createBGF);
use General qw(FileTester GetFileTypeStr);
use ManipAtoms qw(SplitAtomsByMol SelectAtoms BuildAtomSelectionString GetMols GroupAtomsByField);

my ($structFile, $saveFile, $SELECT, $selection, $field, $sort_field, $reverse_sort);
my ($ATOMS, $BONDS, $HEADERS, $MOLS);

$|++;
&init;
print "Parsing structure file $structFile...";
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($structFile, 1);
if(defined($selection)) {
  print "selecting atoms...";
  $SELECT = SelectAtoms($selection, $ATOMS);
  die "No selected atoms found in bgf file\n" if (! keys %{ $SELECT });
}
print "Done\nGrouping atoms based on $field...";
&GroupAtomsByField($ATOMS, $BONDS, $field, $SELECT, $sort_field, $reverse_sort);
print "Done\nCreating bgf file $saveFile...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
	my (%OPTS, $atomSelect, $usage);

	&getopt('bafsgr',\%OPTS);
	$usage = showUsage();
	die "$usage\n" if (! exists($OPTS{b}));

	print "Initializing...";
	($structFile, $saveFile, $atomSelect, $field, $sort_field, $reverse_sort) = ($OPTS{b}, $OPTS{s}, $OPTS{a}, $OPTS{f}, $OPTS{g}, $OPTS{r});
	if (! defined($saveFile)) {
		$structFile =~ /^\s*(\S+)/;
		$saveFile = basename($structFile);
		$saveFile =~ s/\.\w+//;
		$saveFile .= "_mod.bgf";
	}
	
	$field = "CHAIN" if (! defined($field) or $field !~ /^ATMNAME|LABEL|CHAIN|RESNUM|RESID|RES|MOL|MOLECULE|MOLID|MOLECULEID|MOLSIZE|NUMBONDS|FFTYPE|COORD$/i);
	$field = uc($field);
	$field = "RESNUM" if ($field =~ /RESID/);
	$field = "MOLECULEID" if ($field =~ /MOL/ && $field !~ /SIZE/);
	$field = "MOLSIZE" if ($field =~ /SIZE/);
	undef $sort_field if (defined($sort_field) and $sort_field !~ /^ATMNAME|LABEL|CHAIN|RESNUM|RESID|RES|MOL|MOLECULE|MOLID|MOLECULEID|MOLSIZE|NUMBONDS|FFTYPE|COORD/i);
	if (defined($sort_field)) {
		$sort_field = uc($sort_field);
		$sort_field = "RESNUM" if ($sort_field =~ /RESID/);
		$sort_field = "MOLECULEID" if ($sort_field =~ /MOL/ && $sort_field !~ /SIZE/);
		$sort_field = "MOLSIZE" if ($sort_field =~ /SIZE/);
		$reverse_sort = 0 if (! defined($reverse_sort) or $reverse_sort !~ /1|yes/i);
		$reverse_sort = 1 if ($reverse_sort =~ /1|yes/i);
	}
	if (defined($atomSelect)) {
		$selection = BuildAtomSelectionString($atomSelect);
	}
	print "Done\n";
}

sub showUsage {
	my $fTypeStr = GetFileTypeStr;
	my ($usage) = <<DATA;
usage: $0 -b struct_file -a (atom_selection) -f (field) -s (save_name) -g (sort_field) -r (reverse_sort)
Options:
		-b struct_file (Required): location of atom structure file
$fTypeStr		
		-a atom_selection (Optional): any valid field expression. E.g. "resname eq 'WAT' and zcoord > 50"
		-f field (Optional): chain|resid|resname|molecule(id)|molsize|numbonds. Default: chain
		-g sort_field: (Optional): field to sort the grouped atoms by. Default is index
		-r reverse_sort: (Optional): reverse sort the sort_field. Default is 0
		-s save_name (Optional): name of new bgf file
DATA

	return $usage;
}
