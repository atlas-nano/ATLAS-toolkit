#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use File::Basename qw(basename);
use General qw(FileTester STDev TrjSelections ReadParmFile);
use LAMMPS qw(ParseLAMMPSLogFile);
use constant kb => 1.3806504E-23; # boltzmann constant m^2 kg s^-2 K^-1
use Getopt::Std qw(getopt);
use Math::FFT;
use constant PI => atan2(1,1) * 4;

sub init;
sub writeData;
sub saveData;
sub correlation;
sub getFactor;
sub getData;
sub pruneData;
sub getNorm;

my ($logFile, $saveFile, $SELECT, $maxCorr, $cellVol, $sysTemp, $tStep, $normalize);
my ($FIELDS, $DATA, $autoCORR, $pts);
my ($TSTEP, $VOL,  $TEMP);

$|++;
&init;
print "Parsing LAMMPS log file $logFile...";
ParseLAMMPSLogFile($logFile, $SELECT, \&saveData, \*OUTDATA);
print "Done\nCalculating thermal conductivity correlation functions...";
($pts, $maxCorr) = &pruneData($DATA, $maxCorr);
$autoCORR = correlation($DATA, $pts);
print "Done\nWriting data to $saveFile...";
&writeData($autoCORR, $TSTEP, $maxCorr, $normalize, $saveFile);
print "Done\n";

sub writeData {
    my ($data, $tData, $len, $doNorm, $outName) = @_;
    my ($i, $j, @list, $tS, %sum, $factor);

    @list = keys %{ $data };
    #$factor = getFactor();
    open OUTDATA, "> $saveFile" or die "ERROR: Cannot create $saveFile: $!\n";
    printf OUTDATA "#%7s ", "TSTEP";
    for $i (@list) {
	for $j (@list) {
	    next if (! exists($data->{$i}{$j}));
	    printf OUTDATA "%12s %12s ", "${i}${j}", "int_${i}${j}";
	}
    }
    print OUTDATA "\n";
    for $tS (0 .. ($len - 1)) {
	printf OUTDATA "%-8.3f ", ($TSTEP->[$tS]*$tStep/1000);
	for $i (@list) {
	    for $j (@list) {
		next if (! exists($data->{$i}{$j}));
		$sum{$i}{$j} += $data->{$i}{$j}[$tS];
		printf OUTDATA "%12.8G %12.8G ", $data->{$i}{$j}[$tS], $sum{$i}{$j};
	    }
	}
	print OUTDATA "\n";
    }
    close OUTDATA or die "ERROR: Cannot close $saveFile: $!\n";

    open OUTDATA, "> test.dat" or die "ERROR: Cannot write to test.dat: $!\n";
    @list = keys %{ $DATA };
    printf OUTDATA "%-12s ", "#TSTEP";
    for $i (@list) {
	printf OUTDATA "%12s ", $i;
    }
    print OUTDATA "\n";
    for $i (0 .. $#{ $TSTEP }) {
	printf OUTDATA "%12d ", $TSTEP->[$i];
	for $j (@list) {
	    printf OUTDATA "%12.8f ", $DATA->{$j}[$i];
	}
	print OUTDATA "\n";
    }
    close OUTDATA;
}

sub getFactor {
    my ($factor);
    my ($avgTemp, $avgVol, $i);
    
    if (defined($sysTemp)) {
	$avgTemp = $sysTemp;
    } else {
	for $i (@{ $TEMP }) {
	    $avgTemp += $i;
	}
	$avgTemp /= scalar(@{ $TEMP });
    }

    if (defined($cellVol)) {
	$avgVol = $cellVol;
    } else {
	for $i (@{ $VOL }) {
	    $avgVol += $i;
	}
	$avgVol /= scalar(@{ $VOL });
    }

    $factor = 1/($avgVol*kb*$avgTemp**2); # 1/VkbT^2

    return $factor;
}

sub correlation {
    my ($data, $totPts) = @_;
    my ($i, $j, $k, @flist, $numSets); 
    my ($set1, $set2, $fft, $CORR, $norm);

    @flist = keys %{ $data };
    $numSets = scalar @flist;

    for $i (0 .. ($numSets - 1)) {
	$set1 = $data->{ $flist[$i] };
        $fft = new Math::FFT($set1);
	for $j ($i .. ($numSets - 1)) {
	    $set2 = $data->{ $flist[$j] };
	    $norm = getNorm($set1, $set2) if ($normalize);
	    $CORR->{$flist[$i]}{$flist[$j]} = $fft->correl($set2);
	    &normalize($CORR->{$flist[$i]}{$flist[$j]}, $norm) if ($normalize);
        }1
    }

    return $CORR;
}

sub normalize {
    my ($data, $norm) = @_;
    my ($i);

    for $i (0 .. $#{ $data }) {
	$data->[$i] /= $norm;
    }
}

sub getNorm {
    my ($data1, $data2) = @_;
    my ($sumAB,$sumA,$sumB,$norm,$i,$tot);
    # normalization constant = sum(A(i)B(i))/sqrt(sum(A(i)^2) sum(B(i)^2))

    $tot = $#{ $data1 };
    
    for $i (0 .. $tot) {
	$sumAB += $data1->[$i]*$data2->[$i];
	$sumA += $data1->[$i]**2;
	$sumB += $data2->[$i]**2;
    }

    $norm = $sumAB/sqrt($sumA*$sumB);
    
    return $norm;
}
    
sub pruneData {
    my ($data, $corrLen) = $_[0];
    my ($dataLen, $i, $rec, $j, @tmp);

    @tmp = keys %{ $data };
    $dataLen = scalar(@{ $data->{pop @tmp}  });
    #find the smallest power of 2 to include the data
    $i = 2;
    while ($i < $dataLen) {
        $i *=2;
    }
    $i /= 2;

    for $j (($i + 1) .. $dataLen) {
	for (keys %{ $data }) {
	    pop @{ $data->{$_} };
	}
	pop @{ $TSTEP } if (defined($TSTEP));
	pop @{ $VOL } if (defined($VOL));
    }

    $corrLen = $i/2 if (! defined($corrLen) or $corrLen > $i/2);
    return ($i, $corrLen);
}

sub saveData {
    my ($data, $ts, $FLEPTR) = @_;
    my ($i);
    
    push @{ $TSTEP }, $ts;
    die "ERROR: Need the system volume!\n" if (! defined($cellVol) and ! exists($data->{volume}));
    die "ERROR: Need the system temperature!\n" if (! defined($sysTemp) and ! exists($data->{temp}));
    push @{ $VOL }, $data->{volume} if (exists($data->{volume}));
    push @{ $TEMP }, $data->{temperature} if (exists($data->{temp}));
    for $i (keys %{ $data }) {
	next if (! exists($FIELDS->{$i}));
	push @{ $DATA->{$i} }, $data->{$i};
    }
}
    
sub init {
    my ($OPTS, $trjList, $fieldList, $fields, %tmp);

    $fields = { "jx" => 1, "jy" => 1,"jz" => 1};
    getopt('ps',\%tmp);
    die "usage: $0 -p parmfile -s (save file)\n" if (! exists($tmp{p}));
    $OPTS->{LOGFILE}{REQUIRED} = 1;
    &ReadParmFile($tmp{p}, $OPTS);
    
    print "Initializing...";
    ($logFile, $saveFile, $trjList, $fieldList, $maxCorr, $cellVol, $sysTemp, $tStep, $normalize) = 
	($OPTS->{LOGFILE}{VAL}, $tmp{s}, $OPTS->{TRAJSEL}{VAL}, $OPTS->{FIELDS}{VAL}, $OPTS->{CORRLEN}{VAL},
	 $OPTS->{VOL}{VAL}, $OPTS->{TEMP}{VAL}, $OPTS->{TSTEP}{VAL}, $OPTS->{NORMALIZE}{VAL});
    FileTester($logFile);
        
    if (defined($trjList) and $trjList ne "*") {
	$SELECT = TrjSelections($trjList);
	die "ERROR: No valid timestep selection found! Expected :Ita-b:c, got \"$trjList\"\n"
	    if (! $SELECT);
    } 
    
    if (! defined($saveFile)) {
	$saveFile = basename($logFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_lammps_log.dat";
    }
    if (defined($fieldList)) {
	while ($fieldList =~ /(\w+)/g) {
	    if (exists($fields->{lc $1})) {
		$FIELDS->{lc $1} = 1;
	    }
	}
    } else {
	$FIELDS = $fields;
    }
    die "ERROR: No valid fields found while searching \"$fieldList\"!\n" if (! $FIELDS);
    delete $FIELDS->{volume};

    if (defined($maxCorr)) {
	if ($maxCorr =~ /(\d+)/) {
	    $maxCorr = $1;
	} else {
	    undef($maxCorr);
	}
    }
    if (defined($cellVol)) {
	if ($cellVol =~ /(\d+\.?\d*)/) {
	    $cellVol = $1;
	} else {
	    undef($cellVol);
	}
    }

    if (defined($sysTemp)) {
	if ($sysTemp =~ /(\d+\.?\d*)/) {
	    $sysTemp = $1;
	} else {
	    undef($sysTemp);
	}
    }

    if (defined($tStep) and $tStep =~ /^(\d+\.?\d*)/) {
	$tStep = $1;
    } else {
	$tStep = 1; # 1 fs
    }

    if (defined($normalize) and $normalize =~ /^0|no/i) {
	$normalize = 0;
    } else {
	$normalize = 1;
    }
    print "Done\n";
}
