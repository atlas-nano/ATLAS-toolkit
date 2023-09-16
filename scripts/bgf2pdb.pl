#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo createPDB);

# This program will open a bgf file and will write an msi file

sub init;

my ($bgfFile, $pdbFile);
my ($ATOMS, $BONDS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
print "Done\nCreating PDB file $pdbFile...";
createPDB($ATOMS, $BONDS, $pdbFile);
print "Done\n";

sub init {
    my (%OPTS);
    
    getopt('bp',\%OPTS);
    ($bgfFile, $pdbFile) = ($OPTS{b},$OPTS{p});
    die "usage: $0 -b bgf file -p [pdb file (optional)]\n" if (! defined($bgfFile));
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($pdbFile)) {
	$pdbFile = basename($bgfFile);
	$pdbFile =~ s/\.\w+$/\.pdb/;
    }
    print "Done\n";
}
