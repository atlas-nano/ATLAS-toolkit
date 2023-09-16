#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(CoM FileTester);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Packages::CERIUS2 qw(parseCerius2FF);
use Getopt::Std qw(getopt);

sub init;

my ($bgfFile, $cerius2FF);
my ($ATOMS, $BONDS, $PARMS, $CENTER);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
print "Done\nParsing CERIUS2 forcefield $cerius2FF...";
$PARMS = parseCerius2FF($cerius2FF);
AddMass($ATOMS, $PARMS);
print "Done\n";
$CENTER = CoM($ATOMS);
print "Center of Mass: ";
for ("XCOORD", "YCOORD", "ZCOORD") {
    printf "%8.3f", $CENTER->{$_};
}
print "\n";

sub init {
    my (%OPTS);
    
    getopt('bf',\%OPTS);
    ($bgfFile, $cerius2FF) = ($OPTS{b},$OPTS{f});
    
    for ($bgfFile, $cerius2FF) {
	die "usage: $0 -b bgf file -f cerius2 force field\n" if (! defined($_));
    }

    print "Initializing...";
    for ($bgfFile, $cerius2FF) {    
	FileTester($_);
    }
}
