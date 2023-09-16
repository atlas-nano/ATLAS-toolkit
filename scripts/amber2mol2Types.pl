#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use FileFormats qw(ParseStructFile createMOL2 AMBER2MOL2Types);
use General qw(FileTester);

sub init;

my ($struct_file, $saveFile);

$|++;
&init;
my ($ATOMS, $BONDS) = ParseStructFile($struct_file, 0);
print "Done\nConverting atom types...";
&AMBER2MOL2Types($ATOMS);
print "Done\nCreating MOL2 file $saveFile...";
createMOL2($ATOMS, $BONDS, $saveFile, 1);
print "Done\n";

sub init {
    my (%OPTS);
    getopt('ms',\%OPTS);
    die "usage: $0 -m struct_file -s [save file]\n" if (! exists($OPTS{m}));
    print "Initializing...";
    ($struct_file, $saveFile) = ($OPTS{m}, $OPTS{s});
    FileTester($struct_file);
    if (! defined($saveFile)) {
	$saveFile = basename($struct_file);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_converted.mol2";
    }
    print "Done\n";
}
	
