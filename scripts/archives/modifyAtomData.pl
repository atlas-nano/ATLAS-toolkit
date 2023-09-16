#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use Packages::General qw(GetSelections FileTester);
use Packages::ManipAtoms qw(GetAtmList);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub showUsage;
sub init;
sub updateAtomFields;
sub parseFieldList;

my ($FILES, $FIELDS, $selection, $SELECT);
my ($ATOMS, $BONDS, $HEADERS, $SELECTATMS);

$|++;
&init;
print "Parsing Structure file $FILES->{STRUCTURE}...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($FILES->{STRUCTURE},1);
&parseFieldList($ATOMS, $FIELDS);
print "Done\nParsing atom/residue selection...";
$SELECT = GetSelections($selection, 0);
$SELECTATMS = GetAtmList($SELECT, $ATOMS);
die "ERROR: Not a valid atom selection!" if (! keys %{ $SELECTATMS });
print "Done\nUpdating fields..";
&updateAtomFields($ATOMS, $SELECTATMS, $FIELDS);
print "Done\nCreating $FILES->{SAVE}...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $FILES->{SAVE});
print "Done\n";

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
    my (%OPTS, $usage, $atomSelect, $fieldStr, $tmp, $rec);

    getopt('safw',\%OPTS);

    $usage = &showUsage;
    for ("s","a","f") {
	die "$usage\n" if (! defined($OPTS{$_}));
    }
    print "Initializing...";
    ($FILES->{STRUCTURE}, $atomSelect, $fieldStr, $FILES->{SAVE}) = 
	($OPTS{s}, $OPTS{a}, $OPTS{f}, $OPTS{w});
    FileTester($FILES->{STRUCTURE});
    if (! defined($FILES->{SAVE})) {
	$FILES->{SAVE} = $FILES->{STRUCTURE};
	$FILES->{SAVE} =~ s/\.\w+$/_mod\.bgf/;
    }
    if ($atomSelect =~ /\s+/) {
        @{ $selection } = split /\s+/, $atomSelect;
    } else {
        $selection->[0] = $atomSelect;
    }

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

    die "ERROR: Invalid string found while searching \"FIELD\".\n" .
	"Expected field:(+|-|.)new_val. Got \"${fieldStr}\"!\n" if (! defined($FIELDS));

    print "Done\n";
}

sub showUsage {
    return "usage: $0 -s structure file -a atom selection -f field -w [save name]\n" .
	"Options:\n\t-s structure file: Location of structure file\n" .
	"\t-a atom selection: : [:][I|N|T][a|r]x[-y:z]\n" .
	"\t-f field: field to adjust. Expected field:[+|-|.]new_val. Enclose multiple in quotes\n" .
	"\t-w [save name]: (Optional) Name to save the file as\n";
}
