#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::FileFormats;
use Packages::General;

sub getDimension;
sub initialize;

die "usage: $0 bgffile angle dimension(x,y,z,xy,xz...) [savename]\n"
    if (! @ARGV or $#ARGV < 2);

my ($bgf_file, $rot_angle, $dimension, $save_name) = @ARGV;
my ($cList, $counter, @Angles);

$|++;

initialize;
my ($Atom_Info, $CONN, $HEADERS) = GetBGFFileInfo($bgf_file, 1);

$cList = getDimension($dimension);
@Angles = ($rot_angle, $rot_angle, $rot_angle);

for $counter (keys %{ $cList }) {
    Rotate($Atom_Info, \@Angles, $counter);
}
print "Creating BGF file $save_name...";
addHeader($Atom_Info, $HEADERS);
createBGF($Atom_Info, $CONN, $save_name);
print "Done\n";

sub initialize {
    my ($index);

    FileTester($ARGV[0]);
    for $index (0 .. $#ARGV) {
	$ARGV[$index] =  Trim($ARGV[$index]);
    }

    die "Error: Expected decimal(integer) for rotation angle, got $ARGV[1]\n"
	if (! IsInteger($ARGV[1]) and ! IsDecimal($ARGV[1]));

    if (! defined($save_name)) {
	$save_name = "out.bgf";
    }

}

sub getDimension {
    my ($dimension) = $_[0];
    my (%COORD, $counter);

    $counter = 0;
    while ($dimension =~ /(\w)/g) {
	if (lc($1) eq "x" or $1 eq "1") {
	    $COORD{0} = "";
	    $counter++;
	} elsif (lc($1) eq "y" or $1 eq "2") {
	    $COORD{1} = "";
	    $counter++;
	} elsif (lc($1) eq "z" or $1 eq "3") {
	    $COORD{2} = "";
	    $counter++;
	}
    }

    if ($counter == 3) {
	%COORD = ();
	$COORD{3} = "";
    }

    return \%COORD;
}
