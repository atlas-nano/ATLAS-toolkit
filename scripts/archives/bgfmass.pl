#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Packages::CERIUS2 qw(parseCerius2FF);

sub init;
sub calcMass;

die "usage: $0 bgf_file force_field\n"
    if (! @ARGV || $#ARGV < 1);

my ($bgfFile, $ff) = @ARGV;
my ($ATOMS);

$|++;
print "Initializing...";
&init;
print "Done\n";
calcMass($ATOMS);

sub init {
    my ($BONDS, $FFDATA);

    FileTester($bgfFile);
    FileTester($ff);

    print "Done\nParsing $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    print "Done\nParsing $ff...";
    $FFDATA = parseCerius2FF($ff);
    print "Done\nAdding Mass...";
    AddMass($ATOMS, $FFDATA);
}

sub calcMass {
    my ($atoms) = $_[0];
    my ($i, $mass);

    $mass = 0;
    for $i (keys %{ $atoms }) {
	$mass += $atoms->{$i}{MASS};
    }
    
    print "Total Mass: $mass amu\n";
}
