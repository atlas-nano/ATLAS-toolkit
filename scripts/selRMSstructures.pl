#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);
use File::Basename qw(basename);

sub init;
sub parseDatFile;
sub selectStructs;
sub writeStructData;
sub numerically { ($a<=>$b) }

my ($datFile, $num, $crms, $saveFile, $increments);
my ($DATA, $STRUCT);

$|++;
&init;
print "Parsing CRMS dat file $datFile...";
$DATA = parseDatFile($datFile);
print "Done\nSelecting top $num structures within $crms A...";
$STRUCT = selectStructs($DATA, $num, $crms, $increments);
print "Done\n";
&writeStructData($STRUCT, $saveFile);

sub writeStructData {
    my ($dat, $save) = @_;
    my ($snaps, $i);

    open OUTFILE, "> $save" or die "ERROR while writing to $save: $!\n";
    printf OUTFILE "%-8s%12s\n","#SNAP", "CRMS";
    print "Top " . scalar(keys %{ $dat }) . " data points:\n";
    @{ $snaps } = sort numerically keys %{ $dat };
    for $i (@{ $snaps }) {
	printf OUTFILE "%-8s%12.3f\n", $i, $dat->{$i};
	printf "%-8s%12.3f\n", $i, $dat->{$i}; 
    }
    close OUTFILE;
}

sub selectStructs {
    my ($dat, $num, $crms, $incr) = @_;
    my (@tmp, $i, $BEST, $inc, $count, $USED, $snap, $isvalid, $j);

    if (scalar(keys %{ $dat }) <= $num) {
	print "wanted $num snapshot but only " . scalar(keys %{ $dat }) . " read in...";
	for $i (keys %{ $dat }) {
	    $BEST->{$dat->{$i}[0]} = $i;
	}
	return $BEST;
    }
    
    for $i (keys %{ $dat }) {
	delete $dat->{$i} if ($i > $crms);
    }
    die "ERROR: No snapshot matched criteria of CRMS <= $crms\n" if (! keys %{ $dat });

    if (scalar(keys %{ $dat }) <= $num) {
	print "wanted $num snapshot but only " . scalar(keys %{ $dat }) . " matched criteria...";
	for $i (keys %{ $dat }) {
	    $BEST->{$dat->{$i}[0]} = $i;
	}
	return $BEST;
    }

    @tmp = sort numerically keys %{ $dat };
    $inc = int(scalar(@tmp)/$num) + 1; # increments
    $i = 0;
    while ($i <= $#tmp and $i < $num) {
	$isvalid = 0;
	for $j (0 .. $#{  $dat->{$tmp[$i]} }) {
	    $snap = $dat->{$tmp[$i]}[$j];
	    if (! exists($USED->{$snap}) ) {
		for (($snap - $incr) .. ($snap + $incr)) {
		    $USED->{$_} = 1;
		}
		$isvalid = 1;
		last;
	    }
	}
	if (! $isvalid) {
	    splice @{ @tmp }, $i, 1;
	    next;
	}
	$BEST->{$snap} = $tmp[$i];
	$count++;
	$i++;
	#$i += $inc;
    }
    $BEST->{ $dat->{ $tmp[$#tmp] }[0] } = $tmp[$#tmp] if ($count < $num);
    print "found best " . scalar(keys %{ $BEST }) . " points that matched criteria...";
    return $BEST;
}

sub parseDatFile {
    my ($dat) = $_[0];
    my ($DATA);

    open DATAFILE, $dat or die "ERROR while reading $dat: $!\n";
    while (<DATAFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+)\s*(\d+\.?\d*)\s*$/) {
	    push @{ $DATA->{$2} }, $1;
	}
    }
    close DATAFILE;
    die "ERROR: $dat does not contain any valid data!\n" if (! $DATA);
   
    return $DATA;
}

sub init {
    my (%OPTS);
    getopt('dcnsi',\%OPTS);
    die "usage: $0 -d data file -n (num structures = 10) -c (max crms = 2.5) -i (increments = 5) -s (save file)\n"
	if (! exists($OPTS{d}));
    print "Initializing...";
    ($datFile, $num, $crms, $saveFile, $increments) = ($OPTS{d}, $OPTS{n}, $OPTS{c}, $OPTS{s}, $OPTS{i});
    $num = 10 if (! defined($num));
    $crms = 2.5 if (! defined($crms));
    $increments = 5 if (! defined($increments));
    die "ERROR: Cannot access $datFile: $!\n" if (! -e $datFile or ! -r $datFile or ! -T $datFile);
    die "ERROR: Expected positive integer for num structures. Got \"$num\"\n"
	if ($num !~ /^\s*\d+\s*$/);
    die "ERROR: Expected positive integer for increments. Got \"$increments\"\n"
	if ($increments !~ /^\s*\d+\s*$/);   
    die "ERROR: Expected positive integer/decimal for max crms. Got \"$crms\"\n"
	if ($crms !~ /^\s*\d+\.?\d+\s*$/);
    if (! defined($saveFile)) {
	$saveFile = basename($datFile);
	$saveFile =~ s/\.\w+$//;
	$saveFile .= "_bestCRMS.dat";
    }
    print "Done\n";
}
