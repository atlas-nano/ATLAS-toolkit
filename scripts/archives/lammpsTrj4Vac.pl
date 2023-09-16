#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General qw(ParseSelect FileTester PrintProgress);
use Packages::LAMMPS qw(ParseLAMMPSLogFile ParseLammpsTrajectoryFile);

sub init;
sub parseLAMMPSDataFile;
sub saveData;
sub calcVolume;

die "usage: $0 velocity_file log_file saveName [\"frameStart frameEnd frameInterval\"]\n"
    if (! @ARGV || $#ARGV < 2);

my ($velFile, $logFile, $saveName, $selection) = @ARGV;
my ($DATA, $SELECT, $tmp, $TSList, $i, $thermo);

$|++;

$SELECT = ParseSelect($selection);
print "Initializing...";
&init;
print "Done\nParsing LAMMPS log file $logFile...";
$DATA = ParseLAMMPSLogFile($logFile);
print "Done\n";
open VACFILE, "> $saveName" or die "ERROR: Cannot create file $saveName:$!\n";
parseLAMMPSDataFile($velFile, \*VACFILE, $DATA, "Trajectory", \&saveData, $SELECT);
close VACFILE;
print "Created VAC compatible file $saveName\n";


sub init {
    FileTester($logFile);
    FileTester($velFile);
    die "ERROR: saveName cannot be the same as the velocity_file" if ($saveName eq $velFile);
    print "WARNING: file $saveName already exists. Overwriting..." if (-e $saveName);
    for ("Press","TotEng","PotEng","KinEng") {
	$thermo->{$_} = 0.0;
    }
    $thermo->{Temp} = 300.0;
}

sub parseLAMMPSDataFile {
    my ($lammpsFile, $fileHandle, $LOGDATA, $printStr, $doAnal) = @_;
    my ($inStr, $tStepData, @dim, $counter, $tot, $field, $filesize);
    my ($atomC, $currPos, $start, $strLen, $tStep, %header, $rec, $coords);

    if (keys %{ $SELECT }) {
	@dim = keys %{ $SELECT };
	$tot = $#dim + 1;
    } else {
	$tot = 0;
	$filesize = -s $lammpsFile;
    }

    $printStr = "Parsing LAMMPS $printStr file $lammpsFile...";
    $tStepData = ();
    $start = $strLen = 0; 
    $currPos = 1;
    
    @dim =  ("XVEL", "YVEL", "ZVEL");
    %header = ("TIMESTEP"=>1,"NUMBER OF ATOMS"=>1,"BOX BOUNDS"=>1,"ATOMS"=>1);
    open LAMMPS, $lammpsFile or die "ERROR: Cannot open LAMMPS data file $lammpsFile: $!\n";
    while (<LAMMPS>) {
        chomp;
        $inStr = $_;
	if ($inStr =~ /^ITEM: (.+)/) {
	    if (exists($header{$1})) {
		$field = $1;
		if ($field eq "TIMESTEP") {
		    if (keys %{ $tStepData } and (! keys %{ $SELECT } || exists($SELECT->{ $currPos })) ) {
			$doAnal->($tStepData, $LOGDATA, $fileHandle);
			if ($tot) {
			    $strLen = PrintProgress($currPos, $tot, $start, $printStr);
			    $currPos++;
			    last if ($currPos > $tot);
			} else {
			    $currPos = tell(LAMMPS);
			    $strLen = PrintProgress($currPos, $filesize, $start, $printStr);
			}
		    }
		    $tStepData = ();
		    if (! $start) {
			$start = time();
		    }
		}
	    } else {
		undef($field);
	    }
	} elsif (defined($field) && $_ =~ /^\s*(\-?\d+\.?\d*)(\s*\-?\d*.*)$/) {
	    if ($field eq "BOX BOUNDS") {
		$rec = (
			{
			    "lo" => $1,
			    "hi" => $2,
			}
			);
		push @{ $tStepData->{$field} }, $rec;
	    } elsif ($field eq "ATOMS") {
		$coords = $2;
		$atomC = $1;
		$counter = 0;
		while ($coords =~ /\s(\-?\d+\.?\d*e?\-?\d*)/g && ($counter <= $#dim)) {
		    $tStepData->{ATOMS}{$atomC}{ $dim[$counter] } = $1;
		    $counter++;
		}
	    } else {
		$tStepData->{$field} = $1;
	    }
	}
    }
		    
    close LAMMPS or die "ERROR: Cannot close $lammpsFile: $!\n";

    if (keys %{ $tStepData } and (! keys %{ $SELECT } || exists($SELECT->{ $currPos })) ) {
	$doAnal->($tStepData, $LOGDATA, $fileHandle);
    }
    printf "$printStr%-${strLen}s\n", "Done";
		    
}

sub saveData {
    my ($DATA, $LOGDATA, $outFile) = @_;
    my ($i, $tStep, $dim);
    $tStep = $DATA->{TIMESTEP};

    print $outFile "ITEM: TIMESTEP\n$tStep\n";
    print $outFile "ITEM: NUMBER OF ATOMS\n" . $DATA->{"NUMBER OF ATOMS"} . "\nITEM: THERMODYNAMICS\n";
    if (exists($LOGDATA->{$tStep})) {
        for $i (keys %{ $LOGDATA->{$tStep} }) {
	    $thermo->{$i} = $LOGDATA->{$tStep}{$i};
	}
    }
    $thermo->{Volume} = calcVolume($DATA->{"BOX BOUNDS"});
    for $i ("TotEng", "KinEng", "PotEng", "Temp","Volume","Press") {
	printf $outFile "%.3f ", $thermo->{$i};
    }
    printf $outFile "\nITEM: ATOMS\n", $i;
    for $i (1 .. $DATA->{"NUMBER OF ATOMS"}) {
	print $outFile "$i";
	for $dim ("XVEL","YVEL","ZVEL") {
	    if (! exists($DATA->{ATOMS}{$i}{$dim})) {
		print "Not found $i $dim tstep $tStep\n";
		print $outFile " 0.000";
	    } else {
		printf $outFile " %f", $DATA->{ATOMS}{$i}{$dim};
	    }
	}
	print $outFile "\n";
    }
}

sub calcVolume {
    my ($BOX) = $_[0];
    my ($volume, $i);

    for $i (@{ $BOX }) {
	if (! defined($volume)) {
	    $volume = ($i->{hi} - $i->{lo});
	} else {
	    $volume *= ($i->{hi} - $i->{lo});
	}
    }

    return $volume;
}
