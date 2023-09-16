#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts";
}

use strict;
use Packages::General;

sub init;
sub readLAMMPSDataFile;
sub readLAMMPSInputFile;
sub writePairData;
sub numerically;

die "usage: $0 lammps_inputfile lammps_datafile [save_name]\n"
    if (!@ARGV || $#ARGV < 1);

my ($inputFile, $dataFile, $saveFile) = @ARGV;
my ($mixType, $PAIRS);

$|++;

print "Initializing...";
&init;
print "Done\nParsing LAMMPS input file $inputFile...";
($mixType, $PAIRS) = readLAMMPSInputFile($inputFile);
print "Done\nParsing LAMMPS data file $dataFile...";
$PAIRS = readLAMMPSDataFile($dataFile, $PAIRS);
print "Done\nWriting n^2 pair mixing data to $saveFile...";
writePairData($PAIRS, $mixType, $saveFile);
print "Done\n";

sub init {
    FileTester($inputFile);
    FileTester($dataFile);

    $saveFile = "pairmix.dat" if (! $saveFile);
}

sub readLAMMPSInputFile {
    my ($inFile) = $_[0];
    my ($mixType, %PAIRS);

    open INPUTFILE, $inFile or die "ERROR: Cannot read from $inFile: $!\n";
    while(<INPUTFILE>) {
	chomp;
	if ($_ =~ /^pair_coeff\s+(\d+)\s+(\d+)\s+morse\/opt\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $PAIRS{$1}{$2} = (
			      {
				  d0    => $3,
				  alpha => $4,
				  r0    => $5,
			      }
			      );
	} elsif($_ =~ /^pair_modify\s+mix\s+(\S+)/) {
            $mixType = $1;
        } 
    }
    close INPUTFILE;
    die "ERROR: Cannot find mixing type in inputfile\n"
	if (! defined($mixType));
    die "ERRROR: Invalid mixing type in inputfile: $mixType\n"
	if (lc($mixType) !~ /geometric|arithmetic/);
    
    return ($mixType, \%PAIRS);
}

sub readLAMMPSDataFile {
    my ($datFile, $PAIRDATA) = @_;
    my ($startRead);
    
    $startRead = 0;
    open DATAFILE, $datFile or die "ERROR: Cannot open datafile $datFile: $!\n";
    while (<DATAFILE>) {
	chomp;
	if ($_ =~ /^Pair Coeffs/) {
	    $startRead = 1;
	} elsif ($_ =~ /Coeffs/) {
	    if ($startRead) {
		last;
	    } else {
		$startRead = 0;
	    }
	} elsif ($_ =~ /Atoms/) {
	    last;
	} elsif ($startRead and $_ =~ /morse/ and $_ =~ /^\s*(\d+)\s+\D*\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)/) {
	    $PAIRDATA->{$1}{$1} = (
				   {
				       d0    => $2,
				       alpha => $3,
				       r0    => $4,
				   }
				   );
	}
    }
    close DATAFILE;
    
    die "ERROR: No pair data found in either data or input files\n"
	if (! keys %{ $PAIRDATA });
    return $PAIRDATA;
}

sub writePairData {
    my ($PAIRDATA, $mixType, $saveName) = @_;
    my (@types, $i, $j);
    
    @types = sort numerically keys %{ $PAIRDATA };
    
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create data file $saveName: $!\n";
    for $i (@types) {
	next if (! exists($PAIRDATA->{$i}));
	for $j (@types) {
	    next if (! exists($PAIRDATA->{$j}));
	    next if ($i >= $j);
	    if (! exists($PAIRDATA->{$i}{$j})) {
		$PAIRDATA->{$i}{$j}{d0} = 
		    ($PAIRDATA->{$i}{$i}{d0} * $PAIRDATA->{$j}{$j}{d0});
		$PAIRDATA->{$i}{$j}{d0} = sqrt($PAIRDATA->{$i}{$j}{d0});
		$PAIRDATA->{$i}{$j}{alpha} = 
		    ($PAIRDATA->{$i}{$i}{alpha} + $PAIRDATA->{$j}{$j}{alpha})/2;
		if ($mixType eq "geometric") {
		    $PAIRDATA->{$i}{$j}{r0} =
			sqrt($PAIRDATA->{$i}{$i}{r0} * $PAIRDATA->{$j}{$j}{r0});
		} else {
		    $PAIRDATA->{$i}{$j}{r0} =
			($PAIRDATA->{$i}{$i}{r0} + $PAIRDATA->{$j}{$j}{r0})/2;
		}
		printf OUTDATA "%-15s %-4d %-4d %-18s%12.6f %12.6f %12.6f\n", "pair_coeff", $i, $j, "morse/opt",
		$PAIRDATA->{$i}{$j}{d0}, $PAIRDATA->{$i}{$j}{alpha}, $PAIRDATA->{$i}{$j}{r0};
	    }
	}
    }
    close OUTDATA;
}

sub numerically {
    ($a<=>$b);
}
