#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use warnings;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use General qw(FileTester GetSelections ShowSelectionInfo);

sub init;
sub getLMPatomIndex;
sub updateAtomIndex;
sub updateBonds;

my ($bgfFile, $saveFile, $lmpFile);
my ($ATOMS, $BONDS, $HEADERS, $INDICES);

$|++;
&init;
print "Parsing bgf file...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nGetting LAMMPS trajectory atom index from $lmpFile...";
$INDICES = getLMPatomIndex($lmpFile, $ATOMS);
print "Done\nUpdate Atom indicies...";
&updateAtomIndex($ATOMS, $INDICES);
print "Done\nCreating bgf file $saveFile...";
($ATOMS, $BONDS) = updateBonds($ATOMS);
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub updateBonds { 
    my ($atoms) = @_;
    my ($i, $j, $index, $newAtoms, $newBonds);

    $newBonds = ();

    for $i (keys %{ $atoms }) {
        $index = $atoms->{$i}{INDEX}; 
        %{ $newAtoms->{$index} } = %{ $atoms->{$i} };
        $newBonds->{$index} = ();
        next if (! exists($newAtoms->{$index}{BONDS}) || ! $newAtoms->{$index}{BONDS});
        for $j (0 .. $#{ $newAtoms->{$index}{BONDS} }) {
            $newBonds->{$index}[$j] = $newAtoms->{$index}{BONDS}[$j]{INDEX};
        }
        delete $newAtoms->{$index}{BONDS};
    }

    return ($newAtoms, $newBonds);
}

sub updateAtomIndex {
    my ($atoms, $indices) = @_;
    my ($i, $index, $j);

    for $i (keys %{ $atoms }) {
        next if (! exists($BONDS->{$i}) || ! $BONDS->{$i});
        $index = 0;
        for $j (@{ $BONDS->{$i} }) {
            $atoms->{$i}{BONDS}[$index] = \%{ $atoms->{$j} };
            $index++;
        }
    }

    $index = 0;
    for $i (sort {$a<=>$b} keys %{ $indices }) {
	$index++;
	$atoms->{ $indices->{$i} }{INDEX} = $index;
    }
}

sub getLMPatomIndex {
    my ($inFile, $atoms) = @_;
    my ($DATA, $i, $start, $index, $tot);

    $tot = scalar(keys %{ $atoms });
    $start = $index = 0;
    open INFILE, $inFile or die "ERROR: Cannot open $inFile: $!\n";
    while (<INFILE>) {
	chomp;
	last if ($start and $_ =~ /^ITEM/);
	if ($_ =~ /^ITEM: ATOMS/) {
	    $start = 1;
	} elsif ($start and $_ =~ /^(\d+)/) {
	    next if(! exists($atoms->{$1}));
	    $index++;
	    $DATA->{$index} = $1;
	}
    }
    close INFILE;
    die "ERROR: Lammpstrj is invalid\n" if (! $index);
    die "ERROR: Lammpstrj ($index) and bgf file ($tot) contains different number of atoms!\n"
	if ($index != $tot);
    return $DATA;
}

sub init {
    my (%OPTS);

    getopt('bls',\%OPTS);
    die "usage: $0 -b bgf_file -l lammps_trj -s (save_file)\n" 
	if (! exists($OPTS{b}) or ! exists($OPTS{l}));

    print "Initializing...";
    FileTester($OPTS{b});
    FileTester($OPTS{l});

    ($bgfFile, $lmpFile, $saveFile) = ($OPTS{b}, $OPTS{l}, $OPTS{s});
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_lmptrj_updated.bgf";
    }
    print "Done\n";
}
