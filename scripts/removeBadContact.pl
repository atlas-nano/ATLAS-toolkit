#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms);
use General qw(GetSelections FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetAtmList GetMols);

sub usage;
sub getCons;
sub addMolID;
sub invertSelection;

my ($bgfFile, $saveName, $selection);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF, $CONS, $tmp, $HEADERS, $BOX, $totAtoms, $totMols);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&GetMols($ATOMS, $BONDS);
print "Done\nParsing atom/residue selection...";
($totAtoms, $totMols, $SELECTIONS) = &addMolID($SELECTIONS, $ATOMS);
$SELECTIONS = invertSelection($SELECTIONS, $ATOMS);
die "ERROR: No valid atoms selected!\n" if (! keys %{ $SELECTIONS });
print "removing $totAtoms atoms ($totMols molecules)...Done\nSelecting relevant atoms...";
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\nCreating BGF file $saveName...";
addHeader($BGF,$HEADERS);
createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub invertSelection {
    my ($select, $atoms) = @_;
    my ($newSelect, $i);

    for $i (keys %{ $atoms }) {
	$newSelect->{$i} = $atoms->{$i}{MOLECULEID} if (! exists($select->{ATOMS}{$i}));
    }

    return $newSelect;
}

sub addMolID {
    my ($select, $atoms) = @_;
    my ($i, $molname, $j, $totA, $totM,$found); 

    for $i (keys %{ $select }) {
	if ($atoms->{$i}{RESNAME} =~ /POP|WAT/i) {
	    $totM++ if (! exists($found->{MOLECULES}{ $atoms->{$i}{MOLECULEID} }));
	    $found->{MOLECULES}{ $atoms->{$i}{MOLECULEID} } = 1;
	    for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
		$select->{$j} = $atoms->{$i}{MOLECULEID};
		$totA++ if (! exists($found->{ATOMS}{$i}));
		$found->{ATOMS}{$i} = 1;
	    }
	}
    }

    return($totA, $totM, $found);
}

sub init {
    my (%OPTS, $select);
    getopt('bos',\%OPTS);
    ($bgfFile, $saveName, $select) = ($OPTS{b},$OPTS{s},$OPTS{o});
    for ($bgfFile, $select) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    FileTester($select);
    open SELFILE, $select or die "ERROR: Cannot open $select: $!\n";
    while (<SELFILE>) {
	chomp;
	if ($_ =~ /^(\d+)\s+(\d+e?\+?\d*)/) {
	    $SELECTIONS->{$1} = $2 if($2>10);
	}
    }
    close SELFILE;
    die "ERROR: No valid data found while searching $select!" if (! $SELECTIONS);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -o bad_contact_file
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
DATA

die "\n";

}
