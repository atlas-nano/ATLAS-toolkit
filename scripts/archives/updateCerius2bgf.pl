#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use Packages::FileFormats;
use Packages::General;
use strict;

die "usage: $0 original_bgf modified_bgf [save_name]\n"
    if (! @ARGV || $#ARGV < 1);

my ($file1, $file2, $saveName) = @ARGV;

FileTester($file1);
FileTester($file2);

if (! defined($saveName)) {
    $saveName = "out.bgf";
}

my ($ATOMS1, $BONDS1, $HEADER) = GetBGFFileInfo($file1, 1);
my ($ATOMS2, $BONDS2, $HEADER2) = GetBGFFileInfo($file2, 1);
my ($atom, $atmName, $oldAtmName);

for $atom (keys %{ $ATOMS2 }) {
    $atmName = $ATOMS2->{$atom}{"ATMNAME"};
    $oldAtmName = $ATOMS1->{$atom}{"ATMNAME"};
#    $oldAtmName =~ s/\d+//g;
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$ATOMS1->{$atom}{$_} = $ATOMS2->{$atom}{$_};
    }
}

addHeader(\%{ $ATOMS1 }, $HEADER);
createBGF($ATOMS1, $BONDS1, $saveName);
