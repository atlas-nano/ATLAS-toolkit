#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::FileFormats;
use Packages::General;

sub getDimension;
sub initialize;

die "usage: $0 bgffile distance dimension(x,y,z,xy,xz...) [savename]\n"
    if (! @ARGV or $#ARGV < 2);

initialize;

my ($bgf_file, $trans_dist, $dimension, $save_name) = @ARGV;
my ($cList, $counter, $coord);

my ($Atom_Info, $CONN, $HEADERS) = GetBGFFileInfo($bgf_file, 1);

$cList = getDimension($dimension);

for $counter (keys %{ $Atom_Info }) {
    for $coord (keys %{ $cList }) {
	$Atom_Info->{$counter}{$coord} += $trans_dist;
    }
}

print "Creating file $save_name...";
addHeader($Atom_Info, $HEADERS);
createBGF($Atom_Info, $CONN, $save_name);
print "Done\n";

sub initialize {
    my ($index);

    FileTester($ARGV[0]);
    for $index (0 .. $#ARGV) {
	$ARGV[$index] =  Trim($ARGV[$index]);
    }

    die "Error: Expected decimal(integer) for trans distance, got $ARGV[1]\n"
	if (! IsInteger($ARGV[1]) and ! IsDecimal($ARGV[1]));

    if (! $save_name ) {
	$save_name = "out.bgf";
    }

}

sub getDimension {
    my ($dimension) = $_[0];
    my (%COORD);

    while ($dimension =~ /(\w)/g) {
	if (lc($1) eq "x" or $1 eq "1") {
	    $COORD{"XCOORD"} = "";
	} elsif (lc($1) eq "y" or $1 eq "2") {
	    $COORD{"YCOORD"} = "";
	} elsif (lc($1) eq "z" or $1 eq "3") {
	    $COORD{"ZCOORD"} = "";
	}
    }

    return \%COORD;
}
