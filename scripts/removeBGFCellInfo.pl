#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo addHeader createBGF createHeaders);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub usage;

my ($bgfFile, $saveName);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Getting atom information from $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile,1);
print "Creating BGF file $saveName...";
$HEADERS = createHeaders(undef, $saveName);
&addHeader($ATOMS,$HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub init {
    my (%OPTS, $atomSel);

    getopt('bs',\%OPTS);
    ($bgfFile, $saveName) = ($OPTS{b},$OPTS{s});
    die "usage: $0 -b bgf file -s (save name)\n" if (! exists($OPTS{b}));
    print "Initializing...";
    FileTester($bgfFile);
    if (! defined($saveName)) {
        $saveName = basename($bgfFile);
        $saveName =~ s/\.\w+$//;
	$saveName .= ".nocell.bgf";
    }
    print "Done\n";
}
