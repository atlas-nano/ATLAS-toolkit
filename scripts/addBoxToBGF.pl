#!/usr/bin/perl -w
#
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(GetBGFFileInfo addHeader createBGF);
use BOX qw(GetBox);
use File::Basename;

die "usage: $0 bgfFile [save name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgfFile, $saveName) = @ARGV;

FileTester($bgfFile);

if (! $saveName) {
    $saveName = basename($bgfFile);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_mod.bgf";
}

my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
my ($BOX) = GetBox($ATOMS, undef, undef);

my ($bInfo) = "PERIOD 111\nAXES   ZYX\nSGNAME P 1                  1    1\n" .
    sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n", $BOX->{"X"}{"hi"} - $BOX->{"X"}{"lo"},
	    $BOX->{"Y"}{"hi"} - $BOX->{"Y"}{"lo"}, $BOX->{"Z"}{"hi"} - $BOX->{"Z"}{"lo"}, 
	    $BOX->{"X"}{"angle"}, $BOX->{"Y"}{"angle"}, $BOX->{"Z"}{"angle"}) . 
    "CELLS    -1    1   -1    1   -1    1";
push @{ $HEADERS }, $bInfo;
addHeader(\%{ $ATOMS }, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
		     
