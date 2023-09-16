#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use FileFormats qw(GetBGFFileInfo createHeaders addHeader createBGF);
use General qw(FileTester CombineMols);
use File::Basename qw(basename);

die "usage: $0 bgf_file1 bgf_file2 [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($mol1, $mol2, $saveName) = @ARGV;

FileTester($mol1);
FileTester($mol2);

if (! $saveName) {
    $saveName = basename($mol1);
    $saveName =~ s/\.\w+$//;
    $saveName .= "_combined.bgf";
}
my ($bond);
$|++;
print "Parsing bgf file $mol1...";
my ($ATOMS1, $BONDS,$HEADER) = GetBGFFileInfo($mol1, 1);
print "Done\nParsing bgf file $mol2...";
my ($ATOMS2, $BONDS2) = GetBGFFileInfo($mol2, 0);
print "Done\nCombining molecules...";
my ($ATOMS, $CONN) = CombineMols($ATOMS1, $ATOMS2, $BONDS, $BONDS2);
print "Done\nCreating $saveName...";
&addHeader(\%{ $ATOMS }, $HEADER);
createBGF($ATOMS, $CONN, $saveName);
print "Done\n";
