#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createBGF addHeader);
use General qw(FileTester);
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub syncBonds;

my ($bgfFile, $saveFile);
my ($ATOMS, $BONDS, $HEADERS);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSyncing bonds...";
&syncBonds($BONDS);
print "Done\nCreate BGF file $saveFile...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub syncBonds {
    my ($bondList) = $_[0];
    my ($i, $j, $k, $isFound);

    for $i (keys %{ $bondList }) {
	for $j (@{ $bondList->{$i} }) {
	    $isFound = 0;
	    for $k (@{ $bondList->{$j} }) {
		if ($k == $i) {
		    $isFound = 1;
		    last;
		}
	    }
	    if (! $isFound) {
		push @{ $bondList->{$j} }, $i;
	    }
	}
    }
}

sub init {
    my (%OPTS);
    
    getopt('bs',\%OPTS);
    die "usage: $0 -b bgf file -s (save file)\n" if (! defined($OPTS{b}));

    print "Initializing...";
    ($bgfFile, $saveFile) = ($OPTS{b}, $OPTS{s});
    FileTester($bgfFile);
    if (! defined($saveFile)) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_mod.bgf";
    }
    print "Done\n";
}
