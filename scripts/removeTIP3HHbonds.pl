#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(FileTester);
use FileFormats qw(createBGF GetBGFFileInfo addHeader);
use File::Basename qw(basename);

die "usage: $0 bgfFile [saveName]\n"
    if (! @ARGV);

my ($bgfFile, $saveName) = @ARGV;
my ($i, $j, $bond);

FileTester($bgfFile);
if (! defined($saveName)) {
    $saveName = basename($bgfFile);
    $saveName =~ s/\.\w+$/_mod.bgf/;
}

$|++;
print "Parsing BGF file $bgfFile...";
my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\n";

print "Removing HW-HW Bonds...";
for $i (keys %{ $ATOMS }) {
    next if ($ATOMS->{$i}{"FFTYPE"} !~ /^H/);
    $j = 0;
    while ($j <= $#{ $BONDS->{$i} }) {
	$bond = $BONDS->{$i}[$j];
	if ($ATOMS->{$bond}{"FFTYPE"} =~ /^H/) {
	    splice @{ $BONDS->{$i} }, $j, 1;
	} else {
	    $j++;
	}
    }
}
print "Done\n";

print "Creating BGF file $saveName...";
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveName);
print "Done\n";
