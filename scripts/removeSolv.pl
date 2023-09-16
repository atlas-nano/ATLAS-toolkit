#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF GetBGFAtoms sortByRes);
use General qw(GetSelections FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use ManipAtoms qw(GetAtmList);

sub usage;
sub getCons;
sub init;
sub removeSolvent;
sub randSort;

my ($bgfFile, $saveName, $selection, $num);
my ($ATOMS, $BONDS, $SELECTIONS, $BGF); 
my ($CONS, $tmp, $HEADERS, $BOX);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Done\n";

print "Parsing atom/residue selection...";
$SELECTIONS = GetSelections($selection, 0);
$SELECTIONS = GetAtmList($SELECTIONS, $ATOMS);
die "ERROR: No valid atoms selected!\n" if (! keys %{ $SELECTIONS });
print "Done\nRandomely removing $num solvent molecules...";
$SELECTIONS = &removeSolvent($ATOMS, $SELECTIONS, $num);
($BGF, $CONS, $tmp) = GetBGFAtoms($SELECTIONS, $ATOMS, $BONDS);
die "ERROR: No atoms matched selection\n" if (! keys %{ $CONS });
print "Done\n";

print "Creating BGF file $saveName...";
addHeader($BGF,$HEADERS);
createBGF($BGF, $CONS, $saveName);
print "Done\n";

sub randSort {
    my @a = splice @_;
    for my $i (0 .. $#a) {
        my $j = int rand @a;
        @a[$i, $j] = @a[$j, $i];
    }
    return @a;
}

sub removeSolvent {
    my ($atoms, $solventSelection, $numSolvent) = @_;
    my ($i, $SOLVENT, $RES, $allATOMS, @tmp, $tot, $j);

    for $i (keys %{ $solventSelection }) {
	$SOLVENT->{$i} = $atoms->{$i};
    }

    for $i (keys %{ $atoms }) {
	$allATOMS->{$i} = 1;
    }

    $RES = sortByRes($SOLVENT);
    @tmp = randSort(keys %{ $RES });
    $numSolvent = ($#tmp + 1) if (($numSolvent - 1) > $#tmp);
    for $i (0 .. ($numSolvent - 1)) {
	for $j (keys %{ $RES->{$tmp[$i]}{ATOMS} }) {
	    delete $allATOMS->{$j};
	}
    }

    return $allATOMS;
}
    
sub init {
    my (%OPTS, $select);
    getopt('basn',\%OPTS);
    ($bgfFile, $saveName, $select, $num) = ($OPTS{b},$OPTS{s},$OPTS{a}, $OPTS{n});
    for ($bgfFile, $select) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    die "ERROR: Expected integer for num for \"$num\"\n" if ($num !~ /^\d+$/);
    if ($select =~ /\s+/) {
        @{ $selection } = split /\s+/, $select;
    } else {
        $selection->[0] = $select;
    }
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -s save_name -a atom selection -n num solvent mols
Arguments:
  bgf_file: name of bgf_file
  save_name: name of file to save
  num solvent mols: number of solvent molecules to remove
  atom selection:
    [^][:][I|T|N][a|r]
    a   - atom
    r   - residue
    IaX - atom number X
    IrX - residue index X
    TaX - atom type X
    NaX - atom name X
    NrX - residue name X
    Use ":" to specify a range, eg. :Tr1-8 :Ia3-66
    Use "^" to exclude a selection. You can use multiple combinations
    range and exclusion enclosed in quotes, eg, "^:TrIP-IM ^:Ia23-45"
    to exclude residues of type IM and IP and atoms 23-45
DATA

die "\n";

}
