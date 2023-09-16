#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "" . $ENV{HOME} . "/scripts";
}

use strict;
use Packages::General;
use Packages::FileFormats;
use File::Basename;

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
    next if ($ATOMS->{$i}{"FFTYPE"} ne "HW");
    $j = 0;
    while ($j <= $#{ $BONDS->{$i} }) {
	$bond = $BONDS->{$i}[$j];
	if ($ATOMS->{$bond}{"FFTYPE"} eq "HW") {
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
