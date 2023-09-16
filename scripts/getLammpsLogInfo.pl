#!/usr/bin/perl -w

use strict;
use FindBin ();
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester STDev TrjSelections);
use LAMMPS qw(ParseLAMMPSLogFile);

sub init;
sub writeData;
sub numerically { ($a<=>$b); }
sub calcAvg;

my ($logFile, $saveFile, $SELECT, $FIELDS, $tstepOffset);
my ($DATA, $i, $doAvg, $STATS);

$|++;
&init;
open OUTDATA, "> $saveFile" or die "ERROR: Cannot create $saveFile: $!\n";
printf OUTDATA "#%7s ", "TSTEP";
for $i (@{ $FIELDS }) {
    printf OUTDATA "%12s ", $i;
}
print OUTDATA "\n";
print "Parsing LAMMPS log file $logFile...";
ParseLAMMPSLogFile($logFile, $SELECT, \&writeData, \*OUTDATA);
print "Done\n";
&calcAvg($DATA, $STATS, \*OUTDATA) if ($doAvg);
printf OUTDATA "#%7s ", "TSTEP";
for $i (@{ $FIELDS }) {
    printf OUTDATA "%12s ", $i;
}
printf OUTDATA "\n";
close OUTDATA or die "ERROR: Cannot close $saveFile: $!\n";


sub writeData {
    my ($data, $i, $FLEPTR) = @_;
    my ($j, $offset);
    
    printf $FLEPTR "%8d ", $i;
    for $j (@{ $FIELDS }) {
	if (! exists($data->{$j})) {
	    printf $FLEPTR "%12s ", " ";
	} else {
	    printf $FLEPTR "%12.4f ", $data->{$j};
	}
	if ($doAvg) {
	    next if (! exists($data->{$j}));
	    $STATS->{$j} .= $data->{$j} . " ";
	    $offset = sprintf("%.0f", $data->{$j});
	    $DATA->{$j}{$offset} = $i;
	}
    }
    print $FLEPTR "\n";
}
    
sub calcAvg {
    my ($data, $stats, $FLEPTR) = @_;
    my ($i, @tsteps, $avg, $stdev);
    my ($j, %BEST, $bestStr, $offset, $tmp);

    print "Calculating statistics...";
    printf $FLEPTR "#%7s ", "AVG";
    $tmp = sprintf("#%7s ", "STDEV");
    $bestStr = sprintf("#%7s ", "BEST");

    for $i (@{ $FIELDS }) {
	next if (! $stats->{$i});
	chomp $stats->{$i};
	($avg, $stdev, undef) = STDev($stats->{$i});
	printf $FLEPTR "%12.4f ", $avg;
	$tmp .= sprintf("%12.4f ", $stdev);
	for $j (keys %{ $data->{$i} }) {
	    $offset = abs(sprintf("%.0f", $avg) - $j);
	    $BEST{$i}{$offset} = $data->{$i}{$j};
	}
    }
    
    for $i (@{ $FIELDS }) {
	if (! keys %{ $BEST{$i} }) {
	    $bestStr .= sprintf("%12s ", " ");
	} else {
	    @tsteps = sort numerically keys %{ $BEST{$i} };
	    $bestStr .= sprintf("%12d ", $BEST{$i}{$tsteps[0]});
	}
    }
    print $FLEPTR "\n$tmp\n$bestStr\n";
    print "Done\n";
}

sub init {
    my (%OPTS, $fieldList, $fL, $trjList, $tmp);

    $fieldList = " TotEng KinEng Temp PotEng E_bond E_angle" . 
	" E_dihed E_impro E_vdwl E_coul E_long Press Volume v_comDist";

    getopt('lsftoa', \%OPTS);
    die "usage: $0 -l log file -f (fields - default all) -t (timestep selection) -s (save name) -a (compute averages = yes)\n"
	if (! defined($OPTS{l}));
    
    print "Initializing...";
    ($logFile, $saveFile, $fL, $trjList, $tstepOffset, $doAvg) = 
	($OPTS{l}, $OPTS{s}, $OPTS{f}, $OPTS{t}, $OPTS{o}, $OPTS{a});

    FileTester($logFile);
    
    $fL = $fieldList if (! defined($fL) or $fL eq "*");
    while ($fL =~ /(\w+)/g) {
	$tmp = $1;
	#if ($fieldList =~ / $tmp /i) {
	    $tmp = lc($tmp);
	    push @{ $FIELDS }, $tmp;
	#}
    }
    die "ERROR: No valid fields found! Got \"$fL\". Expected \"*\" or $fieldList\n"
	if (! $FIELDS );
    
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

    $tstepOffset = 0 if (! defined($tstepOffset) or $tstepOffset !~ /^\d+$/);

    if (defined($doAvg) and $doAvg =~ /^(0|no)/i) {
	$doAvg = 0;
    } else {
	$doAvg = 1;
    }
    print "Done\n";
}
