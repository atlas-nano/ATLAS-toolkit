#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF insertHeaderRemark);
use General qw(FileTester IsInteger);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;
sub checkAtoms;
sub createBond;

my ($bgfFile, $saveName, $atom1, $atom2);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
&checkAtoms($ATOMS, $atom1, $atom2);
print "Done\nCreating Bond between atoms $atom1 and $atom2...";
&createBond($BONDS, $atom1, $atom2);
print "Done\nCreating BGF file $saveName...";
&insertHeaderRemark($HEADERS, "REMARK created bond between atoms $atom1 and $atom2");
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub createBond {
    my ($bonds, $i, $j) = @_;

    push @{ $bonds->{$i} }, $j;
    push @{ $bonds->{$j} }, $i;
}

sub checkAtoms {
    my ($atoms, $i, $j) = @_;

    for ($i, $j) {
	die "ERROR: Atom $_ is not in BGF file!\n"
	    if (! exists($atoms->{$i}));
    }
}

sub init {
    my (%OPTS, $select);
    getopt('bijs',\%OPTS);
    ($bgfFile, $atom1, $atom2, $saveName) = ($OPTS{b},$OPTS{i},$OPTS{j},$OPTS{s});
    for ($bgfFile, $atom1, $atom2) {
        &usage if (! defined($_));
    }
    print "Initializing...";
    FileTester($bgfFile);
    for ($atom1, $atom2) {
	die "ERROR: Expected integer for atom. Got \"$_\"!\n" if (! IsInteger($_));
    }
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$/_mod\.bgf/;
    }
}

sub usage {
    print STDOUT <<DATA;
usage: $0 -b bgf_file -i atom1 -j atom2 -s (save_name)
DATA

die "\n";

}
