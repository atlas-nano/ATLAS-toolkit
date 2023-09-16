#!/usr/bin/perl -w
BEGIN {
    push @INC, "/ul/tpascal/scripts/";
}

use strict;
use Packages::FileFormats;
use Packages::CERIUS2;
use Packages::BOX;

sub initialize;

die "usage: $0 bgf_file par_file [save_name]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgf_file, $par_file, $save_name) = @ARGV;
initialize;

my ($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgf_file, 1);
my ($PARMS) = parseCerius2FF($par_file);
my ($BBOX) = GetBox($ATOMS, $PARMS, $HEADERS);

my ($counter, $isError, $errList, $radii);

for $counter (keys %{ $ATOMS }) {
    $radii = GetRadii(\%{ $ATOMS->{$counter } }, $PARMS);
if ($radii) {
    $ATOMS->{$counter}{"RADII"} = $radii;
    $ATOMS->{$counter}{"RESONANCE"} = 0;
    $ATOMS->{$counter}{"OCCUPANCY"} = 0;
} else {
    $isError = 1;
    $errList .= $ATOMS->{$counter}{"FFTYPE"} . "\n";
    $ATOMS->{$counter}{"RADII"} = 0;
    $ATOMS->{$counter}{"RESONANCE"} = 0;
    $ATOMS->{$counter}{"OCCUPANCY"} = 0;
}
}

addHeader($ATOMS, $HEADERS);
print "Creating $save_name...";
createBGF($ATOMS, $BONDS, $save_name, 1);
print "Done\n";

if ($isError) {
    print "WARNING: The following atoms didn't have a type in the forcefield. Radii is assigned to 0:\n$errList";
}

sub initialize {

    for ($bgf_file, $par_file) {
	die "ERROR: Cannot access file $_:\n"
	    if (! -e $_ or ! -r $_ or ! -T $_);
    }

    if (! $save_name) {
	$save_name = $bgf_file;
    }

}


sub getRadii {
    my ($atom, $PAR) = @_;
    my ($returnVal, $atmName);

    $atmName = $atom->{"FFTYPE"};
    if (! exists($PAR->{"VDW"}{$atmName})) {
	$atmName = $atom->{"ATMNAME"};
    }

    if (exists($PAR->{"VDW"}{$atmName})) {
	$returnVal = $PAR->{"VDW"}{$atmName}{"VALS"}[1]/2;
    }

    return $returnVal;
}
