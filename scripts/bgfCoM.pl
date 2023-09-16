#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(CoM FileTester LoadFFs ReadFFs);
use FileFormats qw(GetBGFFileInfo AddMass);
use CERIUS2 qw(parseCerius2FF);
use Getopt::Std qw(getopt);
use ManipAtoms qw(GetAtmList GetAtmData GetMols SelectAtoms BuildAtomSelectionString);

sub init;

my ($bgfFile, $FFILES, $selection);
my ($ATOMS, $BONDS, $PARMS, $CENTER, $SELECT);

$|++;
&init;
print "Parsing BGF file $bgfFile...";
($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
&GetMols($ATOMS,$BONDS);
$SELECT = SelectAtoms($selection, $ATOMS);
print "Done\n";
if (defined($FFILES)) {
    $PARMS = LoadFFs($FFILES,2);
    AddMass($ATOMS, $PARMS);
}
print "Center of Mass: ";
$CENTER = CoM(GetAtmData($ATOMS, $SELECT));
#$CENTER = CoM($ATOMS);
for ("XCOORD", "YCOORD", "ZCOORD") {
    printf "%15.5f", $CENTER->{$_};
}
print "\n";

sub init {
    my (%OPTS, $select, $fflist);
    
    getopt('bfa',\%OPTS);
    ($bgfFile, $fflist, $select) = ($OPTS{b},$OPTS{f}, $OPTS{a});
    
    die "usage: $0 -b bgf file -f cerius2 force field -a (atom selection. default all)\n" 
	if (! defined($bgfFile));

    print "Initializing...";
    FileTester($bgfFile);
    if (defined($fflist)) {
        ($FFILES, undef) = ReadFFs($fflist, undef);
 	}
    $select = "index>0" if (! defined($select));
    $selection = BuildAtomSelectionString($select);

}
