#!/usr/bin/perl
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo GetMOL2FileInfo GetPDBFileInfo GetBGFAtoms);
use Packages::General qw(FileTester GetSelections);
use Getopt::Std qw(getopt);

my ($readFunc, $dataFile, $SELECT, $fileType);
my ($ATOMS, $BONDS, $HEADERS);
my (%OPTS, $select, $selection);

$|++;
getopt('fts',\%OPTS);
my ($dataFile, $select, $fileType) = ($OPTS{f},$OPTS{s},$OPTS{t});
die "usage: $0 -f data file -t (file type [bgf|pdb|mol2]) -s [atom selection]\n" 
    if (! defined($dataFile));
print "Initializing...";
FileTester($dataFile);
$select = "*" if (! defined($select));
if (! defined($fileType)) {
    if ($dataFile =~ /.*\.(\w+)$/) {
	$fileType = $1;
	if ($fileType !~ /^(bgf|mol2|pdb)$/i) {
	    $fileType = "bgf";
	}
    } else {
	$fileType = "bgf";
    }
}
$fileType = uc($fileType);

if ($fileType =~ /bgf|mol2|pdb/i) {
    $readFunc = eval ('\&Get' . $fileType . 'FileInfo');
} else {
    die "ERROR: Expected bgf|mol2|pdb while parsing filetype. Got $fileType\n";
}

if ($select =~ /\s+/) {
    @{ $selection } = split /\s+/, $select;
} else {
    $selection->[0] = $select;
}

$SELECT = GetSelections($selection, 0);
print "Done\n";
print "Parsing $fileType file $dataFile...";
($ATOMS, $BONDS, $HEADERS) = $readFunc->($dataFile);
print "Done\nSelecting Relevant Atoms...";
($ATOMS, $BONDS) = GetBGFAtoms($SELECT, $ATOMS, $BONDS);
print "Done\nFound " . scalar(keys %{ $ATOMS }) . " atoms\n";
