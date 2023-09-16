#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetPDBFileInfo createPDB);
use General qw(FileTester);

sub init;
sub removeHydrogens;

my ($pdbFile, $saveFile, $writeBonds);
my ($ATOMS, $BONDS);

$|++;
&init;
print "Parsing PDBFile $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
undef $BONDS if (! $writeBonds);
print "Done\nRemoving hydrogens...";
&removeHydrogens($ATOMS);
print "Done\nSaving PDB file $saveFile...";
&createPDB($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub removeHydrogens {
    my ($atoms) = $_[0];
    my ($i);

    for $i (keys %{ $atoms }) {
	delete($atoms->{$i}) if ($atoms->{$i}{ATMNAME} =~ /^\s*\d*H/);
    }
}
    
sub init {
    my (%OPTS);
    getopt('psb',\%OPTS);
    die "usage: $0 -p pdb file -s [save name] -b [write bonds = yes]\n" if (! exists($OPTS{p}));
    print "Initializing...";
    ($pdbFile, $saveFile, $writeBonds) = ($OPTS{p}, $OPTS{s}, $OPTS{b});
    FileTester($pdbFile);
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "nopdb.pdb";
    }
    $writeBonds = 1 if (! defined($writeBonds) or $writeBonds !~ /0|no/i);
    $writeBonds = 0 if ($writeBonds =~ /0|no/i);
    print "Done\n";
}
