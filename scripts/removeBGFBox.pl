#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";

use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo createBGF addHeader);
use File::Basename qw(basename);
use Getopt::Std qw(getopt);

sub init;
		     
my ($bgfFile, $saveName, $i);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nSaving BGF file $saveName...";
$i = 0;
while($i < $#{ $HEADERS }) {
	if($HEADERS->[$i] =~ /^(PERIOD|AXES|SGNAME|CRYSTX|CELLS)/) {
		splice @{ $HEADERS }, $i, 1;
	} else {
		$i++;
	}
}
&addHeader($ATOMS, $HEADERS);
&createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";

sub init {
    my (%OPTS);
  
    getopt('bs',\%OPTS);

    die "usage: $0 -b bgf_file -s (savename)\n" if (! exists($OPTS{b}));
    print "Initializing...";
    ($bgfFile, $saveName) = ($OPTS{b}, $OPTS{s});
    FileTester($bgfFile);

    if (! $saveName) {
	$saveName = basename($bgfFile);
	$saveName =~ s/\.\w+$//;
	$saveName .= "_mod.bgf";
    }
    print "Done\n";
}

