#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(GetPDBFileInfo createPDB);
use General qw(FileTester);

sub init;

my ($pdbFile, $terList, $saveFile);
my ($ATOMS, $BONDS, $JAGDATA);

$|++;
&init;
print "Parsing PDBFile $pdbFile...";
($ATOMS, $BONDS) = GetPDBFileInfo($pdbFile);
print "Done\nCreating PDB file $saveFile...";
&createPDB($ATOMS, $BONDS, $saveFile, $terList);
print "Done\n";


sub init {
    my (%OPTS, $terStr);
    getopt('pst',\%OPTS);
    
    die "usage: $0 -p pdbfile -t (ligand atoms) -s (save name)\n" if (! exists($OPTS{p}) or ! exists($OPTS{t}));
    print "Initializing...";
    ($pdbFile, $saveFile, $terStr) = ($OPTS{p}, $OPTS{s}, $OPTS{t});
    FileTester($pdbFile);
    
    if (! defined($saveFile)) {
	$saveFile = basename($pdbFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.pdb";
    }
    while ($terStr =~ /(\d+)/g) {
	$terList->{$1} = 1;
    }
    print "Done\n";
}
