#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::FileFormats qw(GetMOL2FileInfo createBGF);
use Packages::General qw(FileTester); 
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub initialize;

my ($mol2File, $bgfFile) = @ARGV;
$|++;
&init;
print "Parsing MOL2 File $mol2File...";
my ($FDATA, $CONN) = GetMOL2FileInfo($mol2File, 0);
print "Done\nCreating BGF File $bgfFile...";
createBGF($FDATA, $CONN, $bgfFile);
print "Done\n";
 
sub init {
    my (%OPTS);
    getopt('mb',\%OPTS);
    die "usage: $0 -m mol2file -b [bgffile]\n" if (! defined($OPTS{m}));
    ($mol2File, $bgfFile) = ($OPTS{m},$OPTS{b});

    print "Initializing...";
    FileTester($mol2File);
    if (! defined($bgfFile)) {
	$bgfFile = basename($mol2File);
	$bgfFile =~ s/\..+$/\.bgf/;
    }
    print "Done\n";
}
