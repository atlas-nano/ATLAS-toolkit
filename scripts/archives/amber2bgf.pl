#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use File::Basename;
use Packages::General;
use Packages::FileFormats;
use Packages::AMBER;
use Packages::LAMMPS;

# This program will read AMBER topology and coordinate files 
# and create a biograf file, with the correct charges
# and atoms labels as specified in the conversion file.
# usage: amber2bgf.pl topologyFile coordinateFile [output_bgf] [cnv_file]

sub initialize;

die "usage: amber2bgf.pl topologyFile coordinateFile [output_bgf] [cnv_file]\n"
    if (! @ARGV or $#ARGV < 1);

my ($topFile, $coordFile, $bgf_file, $cnv_file) = @ARGV;
my ($OPTS, $DATA, $RES, $timeStamp, $hasVels, $CONS, $HEADERS, $vels, $totAtms);

$|++;
initialize;
print "Obtaining Topology/Coordinate Information...";
$OPTS = getOpts;
($DATA, $totAtms) = getTopInfo($topFile, $OPTS);
if ($cnv_file) {
    $RES = loadCnvFile($cnv_file); #load the conversion file
    updateAtmLabels(\%{ $DATA->{"ATOMS"} }, $cnv_file, $RES);
}

($timeStamp, $hasVels) = parseCoordFile($DATA, $coordFile, $totAtms); # Add the coordinate information;
if (exists($DATA->{"BOX"})) {
    #imageAtoms(\%{ $DATA->{"ATOMS"} }, \%{ $DATA->{"BOX"} });
}
print "Done\nGenerating Connectivities...";
$CONS = getConn(\%{ $DATA->{"BONDLIST"} }, $DATA->{"ATOMS"});
print "Done\nCreating BGF file $bgf_file...";
if (keys %{ $DATA->{BOX} }) {
    $HEADERS = createHeaders(\%{ $DATA->{"BOX"} }, $bgf_file);
} else {
    $HEADERS = createHeaders(undef, $bgf_file);
}
addHeader(\%{ $DATA->{"ATOMS"} }, $HEADERS);
createBGF(\%{ $DATA->{"ATOMS"} }, $CONS, $bgf_file);
print "Done\n";
if ($hasVels) {
    $vels = GetVels(\%{ $DATA->{"ATOMS"} });
    open VELFILE, "> velocities.dat" or die "ERROR: Cannot write to file velocities.dat: $!\n";
    print VELFILE $vels;
    close VELFILE;
}
    
sub initialize {
    for (0 .. 1) {
	FileTester($ARGV[$_]);
    }

    if (! $bgf_file) {
	$bgf_file = basename($topFile);
	$bgf_file =~ s/\.\w+$//;
	$bgf_file .= ".bgf";
    }

}

