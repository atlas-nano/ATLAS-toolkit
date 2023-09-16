#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::FileFormats qw(GetBGFFileInfo);
use Packages::General qw(FileTester);

sub init;
sub writeGroupData;
sub readLAMMPSDataFile;
sub getValidTypes;

die "usage: $0 bgfFile dataFile [saveName]\n"
    if (! @ARGV or $#ARGV < 1);

my ($bgfFile, $dataFile, $saveName) = @ARGV;

my ($GROUPS, $PARMS, $ATOMS, $BONDS, $HEADERS, $TYPES, $watType);

$|++;

print "Initializing...";
&init;
print "Done\nParsing $bgfFile...";
($ATOMS, $BONDS, $HEADERS) = GetBGFFileInfo($bgfFile);
print "Done\nReading parms file...";
($PARMS, $watType) = getValidTypes($ATOMS);
print "Done\nGetting Atom types...";
($TYPES, $watType) = readLAMMPSDataFile($dataFile, $PARMS, $watType);
print "Done\nWriting group data to $saveName...";
writeGroupData($TYPES, $watType, $saveName);
print "Done\n";

sub init {
    FileTester($bgfFile);
    FileTester($dataFile);

    if (! defined($saveName)) {
	$saveName = "groupdata.dat";
    }
}

sub getValidTypes {
    my ($atoms) = $_[0];
    my (%TYPES, $atomC, $watFF);

    $watFF = "";
    for $atomC (keys %{ $atoms }) {
	if ($atoms->{$atomC}{ATMNAME} =~ /CYT|THY|GUA|ADE|H\d|D\d/) {
	    $TYPES{ $atoms->{$atomC}{FFTYPE} } = 1;
	} elsif ($atoms->{$atomC}{ATMNAME} =~ /WAT/) {
	    $watFF = $atoms->{$atomC}{FFTYPE};
	}
    }

    die "ERROR: BGF file does not contain any valid info!\n" if (! keys %TYPES or ! $watFF);

    return (\%TYPES, $watFF);
}

sub readLAMMPSDataFile {
    my ($inFile, $parms, $watFF) = @_;
    my ($isStart, %PARMS, $watType);

    $isStart = $watType = 0;
    open DATAFILE, $inFile or die "ERROR: Cannot read from LAMMPS datafile $inFile: $!\n";
    while (<DATAFILE>) {
	chomp;
	if ($_ =~ /^Masses/) {
	    $isStart = 1;
	} elsif ($_ =~ /Coeff|Atoms/) {
	    last;
	} elsif ($isStart and $_ =~ /\s*(\d+)\s+\d+\.\d+\s+\#\s(\S+)/) {
	    if (exists($parms->{$2})) {
		$PARMS{$1} = $2;
	    } elsif ($2 eq $watFF) {
		$watType = $1;
	    }
	}
    }
    close DATAFILE;

    die "ERROR: No valid data found in LAMMPS datafile $inFile!\n" if (! keys %PARMS or ! $watType);

    return (\%PARMS, $watType);
}

sub writeGroupData {
    my ($types, $WAT, $fileName) = @_;
    my ($i, $j, $groupName, $glist);

    open OUTDATA, "> $fileName" or die "ERROR: Cannot create file $fileName: $!\n";
    printf OUTDATA "%-15s nuc type", "group";
    for $i (keys %{ $types }) {
	print OUTDATA " $i";
    }
    print OUTDATA "\n";
    printf OUTDATA "%-15s solvent type $WAT\n", "group";
    printf OUTDATA "%-15s nucleic subtract all solvent\n", "group";
    printf OUTDATA "%-15s backbone subtract nucleic nuc\n", "group";
    printf OUTDATA "%-15s nucrigid nuc rigid molecule\n", "fix";
    printf OUTDATA "%-15s movable subtract all nuc\n", "group";
    printf OUTDATA "%-15s exclude molecule nuc\n", "neigh_modify";
    printf OUTDATA "%-15s nuc multi remove\n", "delete_bonds";
    
    close OUTDATA;
}
