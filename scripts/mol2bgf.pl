#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetMOL2FileInfo createBGF addHeader);
use General qw(FileTester); 
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub initialize;

my ($mol2File, $bgfFile) = @ARGV;
$|++;
&init;
print "Parsing MOL2 File $mol2File...";
my ($FDATA, $CONN, $HEADERS) = GetMOL2FileInfo($mol2File, 0);
print "Done\nCreating BGF File $bgfFile...";
&addHeader($FDATA,$HEADERS);
&createBGF($FDATA, $CONN, $bgfFile);
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
	$bgfFile =~ s/\.mol2$/\.bgf/;
    }
    print "Done\n";
}
