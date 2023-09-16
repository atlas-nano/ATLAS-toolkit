#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(ParseStructFile createMOL2);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

my ($ATOMS, $BONDS, $HEADERS, $struct_file, $mol2File, $verbose, $forAmber);
$|++;
&init;
print "Parsing structure file $struct_file..." if ($verbose);
($ATOMS, $BONDS, $HEADERS) = ParseStructFile($struct_file, 1, $verbose);
print "Done\nCreating MOL2 file $mol2File..." if ($verbose);
&createMOL2($ATOMS, $BONDS, $mol2File, $forAmber);
print "Done\n" if ($verbose);

sub init {
    my (%OPTS);
    getopt('bamv',\%OPTS);
    ($struct_file, $mol2File, $forAmber, $verbose) = ($OPTS{b}, $OPTS{m}, $OPTS{a}, $OPTS{v});
    die "usage: $0 -b struct_file -m [mol2file] -a [mol24amber?=no] -v [verbose=yes]\n" if (! defined($struct_file));
    $verbose = 1 if (! defined($verbose) or $verbose !~ /^0|no$/i);
    $verbose = 0 if ($verbose =~ /^0|no$/i);
    print "Initializing..." if ($verbose);
    if (! defined($mol2File)) {
		$struct_file =~ /^\s*(\S+)/;
		$mol2File = basename($1);
		$mol2File =~ s/\.\w+$/\.mol2/;
    }
    $forAmber = 1 if (! defined($forAmber));
	if ($forAmber =~ /1|yes/i) {
		$forAmber = 1;
    } else {
		$forAmber = 0;
    }
    print "Done\n" if ($verbose);
}
