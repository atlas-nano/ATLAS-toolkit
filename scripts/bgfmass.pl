#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use General qw(ReadFFs FileTester LoadFFs GetFileTypeStr);
use FileFormats qw(ParseStructFile AddMass);
use Getopt::Std qw(getopt);
use ManipAtoms qw(GetMols SelectAtoms BuildAtomSelectionString);
use BOX qw(GetBox);

sub init;
sub calcMass;

my ($structFile, $ff, $selection);
my ($ATOMS, $BONDS, $SELECT, $FFILES, $FF, $ffList);

$|++;
&init;
print "Parsing $structFile...";
($ATOMS, $BONDS) = ParseStructFile($structFile, 0);
&GetMols($ATOMS, $BONDS);
$SELECT = SelectAtoms($selection, $ATOMS);
die "ERROR: No valid atoms found!\n" if (! keys %{ $SELECT });
print "Done\n";
($FFILES, undef) = ReadFFs($ffList);
$FF = LoadFFs($FFILES);
print "Adding Mass...";
&AddMass($ATOMS, $FF);
&calcMass($ATOMS, $SELECT);

sub init {
    my (%OPTS, $atomSelect, $fTypeStr);

    getopt('bfa',\%OPTS);
	$fTypeStr = GetFileTypeStr;
    for ("b", "f") {
		die "usage: $0 -b structure_file -f force field -a (atom selection = all)\n\tFile options:\n$fTypeStr"
			if (! exists($OPTS{$_}));
    }

    print "Initializing...";
    ($structFile, $ffList, $atomSelect) = ($OPTS{b}, $OPTS{f}, $OPTS{a});
    $atomSelect =  "index > 0" if (! defined($atomSelect));
    $selection = BuildAtomSelectionString($atomSelect);
    print "Done\n";

}

sub calcMass {
    my ($atoms, $atomSelect) = @_;
    my ($i, $mass);

    $mass = 0;
    for $i (keys %{ $atomSelect }) {
	$mass += $atoms->{$i}{MASS};
    }
    
    print "Total Mass: $mass amu\n";
}
