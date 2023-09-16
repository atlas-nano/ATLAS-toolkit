#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(FileTester TrjSelections);
use Packages::FileFormats qw(GetBGFFileInfo AddMass);
use Packages::CERIUS2 qw(parseCerius2FF);
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsByteOffset);

sub init;
sub calcTemp;

die "usage: $0 lammps_vel bgf_file force_field \"trajectory selection\" [saveName]\n"
    if (! @ARGV || $#ARGV < 3);

my ($velFile, $bgfFile, $ff, $selection, $saveName) = @ARGV;
my ($ATOMS, $SELECT, $printStr);

$|++;
print "Initializing...";
&init;
print "Done\n";
GetLammpsByteOffset($SELECT, $velFile,  scalar keys %{ $ATOMS });
open TEMPFILE, "> $saveName" or die "ERROR: Cannot create file $saveName:$!\n";
$printStr = "Parsing LAMMPS Velocity file $velFile...";
ParseLAMMPSTrj($ATOMS, $velFile, $SELECT, "vel", \&calcTemp, $printStr, \*TEMPFILE);
close TEMPFILE;
print "Saved temperature profile to $saveName\n";

sub init {
    my ($BONDS, $FFDATA);

    FileTester($velFile);
    FileTester($bgfFile);
    FileTester($ff);
    
    print "Done\nParsing $bgfFile...";
    ($ATOMS, $BONDS) = GetBGFFileInfo($bgfFile, 0);
    print "Done\nParsing $ff...";
    $FFDATA = parseCerius2FF($ff);
    print "Done\nAdding Mass...";
    AddMass($ATOMS, $FFDATA);
    $saveName = "temp_profile.dat" if (! defined($saveName));
    $SELECT = TrjSelections($selection);
}

sub calcTemp {
    my ($DATA, $atomData, $fileHandle) = @_;
    my ($temp, $mvv2e, $boltz, $dof, $i, $atom, $mass, $tfactor, $v2, $field);

    $dof = $DATA->{"NUMBER OF ATOMS"}[0] * 3;
    $boltz = .83142;
    $mvv2e = 1;
    $tfactor = ($mvv2e/($dof * $boltz));
    $temp = 0.0;
    for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
	$atom = \%{ $DATA->{ATOMS}{$i} };
	$mass = $atomData->{$i}{MASS};
	$v2 = 0;
	for $field ("XVEL","YVEL","ZVEL") {
	    $v2 += ($atom->{$field} * 1000)**2;
	}
	$temp += $v2 * $mass;
    }

    $temp *= $tfactor;

    printf $fileHandle "%10d %12.3f\n", $DATA->{TIMESTEP}[0], $temp;
}
