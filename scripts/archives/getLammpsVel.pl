#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::LAMMPS qw(ParseLAMMPSTrj GetLammpsTrjType);
use Packages::General qw(FileTester TrjSelections);
use Packages::AMBER qw(CreateAmberTrj);
use File::Basename;
use Getopt::Std qw(getopt);

my ($velFile, $saveName, $saveType, $SELECT, $pStr);

sub init;
sub saveVel;

$|++;
&init;
open OUTDATA, "> $saveName" || die "ERROR: Cannot create $saveName: $!\n";
$pStr = "Parsing LAMMPS Velocity file $velFile...";
if ($saveType eq "lammps") {
    print OUTDATA "Velocities\n\n";
} else {
    print OUTDATA "AMBER velocity file created by $ENV{USER} at " . 
	scalar(localtime) . "\n";
}
&GetLammpsTrjType($SELECT, $velFile, undef);
&ParseLAMMPSTrj($saveType, $velFile, $SELECT, "vel", \&saveVel, $pStr, \*OUTDATA);
close OUTDATA or die "ERROR: Cannot close $saveName: $!\n";

sub saveVel {
    my ($DATA, $type, $fileHandle) = @_;
    my ($i, $j, $BOX);

    if ($type eq "lammps") {
	for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
	    printf $fileHandle "%8d", $i;
	    for $j ("XVEL", "YVEL", "ZVEL") {
		printf $fileHandle "%12.8f", $DATA->{ATOMS}{$i}{$j};
	    }
	    print $fileHandle "\n";
	}
    } else {
	for $i (1 .. $DATA->{"NUMBER OF ATOMS"}[0]) {
	    for $j ("X", "Y", "Z") {
		$DATA->{ATOMS}{$i}{$j . "COORD"} = $DATA->{ATOMS}{$i}{$j . "VEL"} * 
		    1000/20.455;
	    }
	}
	CreateAmberTrj($DATA->{ATOMS},$BOX, $fileHandle);
    }
}

sub init {
    my (%OPTS, $selections);
    
    getopt('vtsf',\%OPTS);
    ($velFile, $selections, $saveType, $saveName) = ($OPTS{v}, $OPTS{f}, $OPTS{t}, $OPTS{s});
    die "usage:$0 -v lammps velocity file -f trajectory frames (* for all)\n" . 
	"\t-t ouput type [=lammps|amber] (optional) " . "-s [saveName] (optional)\n" 
	if (! defined($velFile) || ! defined($selections));

    print "Initializing...";
    FileTester($velFile);
    if (! defined($saveName)) {
	$saveName = basename($velFile);
	$saveName =~ s/\.\w+$/_lammps\.vel/;
    }
    $saveType = "lammps" if (! defined($saveType) || $saveType !~ /^[lammps|amber]/i);
    $saveType = lc($saveType);
    $SELECT = TrjSelections($selections);
    print "Done\n";
}
