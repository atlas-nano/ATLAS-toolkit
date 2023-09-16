#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use Getopt::Std qw(getopt);
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo createBGF addHeader);

sub init;

$|++;

my ($bgfFile, $headerFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS);

&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, undef) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing Header BGF file $headerFile...";
(undef, undef, $HEADERS) = GetBGFFileInfo($headerFile, 1);
print "Done\nAdding header to $saveFile and saving...";
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
    my (%OPTS);

    getopt('bhs',\%OPTS);
    for ("b", "h") {
	die "usage: $0 -b bgf file -h header file -s [save file]\n" if (! exists($OPTS{$_}));
    }
    print "Initializing...";
    ($bgfFile, $headerFile, $saveFile) = ($OPTS{b}, $OPTS{h}, $OPTS{s});
    FileTester($bgfFile);
    FileTester($headerFile);

    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_header.bgf";
    }
    print "Done\n";
}
