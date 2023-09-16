#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo createBGF addHeader);
use Packages::LAMMPS qw (ReadDataFile);
use Packages::General qw (FileTester);
use File::Basename;

sub init;
sub updateCoordinates;

die "usage: $0 datafile bgffile [saveFile]\n"
    if (! @ARGV or $#ARGV < 1);

my ($dataFile, $bgfFile, $saveFile) = @ARGV;
my ($ATOMS, $BONDS, $HEADERS, $DATA);

$|++;

print "Initializing...";
&init;
print "Done\nParsing BGF file $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile, 1);
print "Done\nParsing LAMMPS data file $dataFile...";
$DATA = ReadDataFile($dataFile);
print "Done\nCreating updated BGF file $saveFile...";
updateCoordinates($ATOMS, $DATA->{ATOMS});
addHeader($ATOMS, $HEADERS);
createBGF($ATOMS, $BONDS, $saveFile);
print "Done\n";

sub init {
    FileTester($bgfFile);
    FileTester($dataFile);
    if (! $saveFile) {
	$saveFile = basename($bgfFile);
	$saveFile =~ s/\.\w+$//g;
	$saveFile .= "_updated.bgf";
    }
}

sub updateCoordinates {
    my ($ATOMS, $COORDS) = @_;
    my ($i, $dim, $j, @tmp, @tmp1);

    @tmp =  keys %{ $ATOMS };
    @tmp1 = keys %{ $COORDS };

    die "ERROR: Inconsistent atom count in BGF and DATA files... aborting\n"
	if ($#tmp != $#tmp1);

    for $i (@tmp) {
	$j = 3;
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $ATOMS->{$i}{$dim} = $COORDS->{$i}{$j};
	    $j++;
	}
    }
}
