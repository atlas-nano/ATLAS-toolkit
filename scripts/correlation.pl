#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use Getopt::Std qw(getopt);
use File::Basename qw(basename);
use General qw(FileTester);
use Math::FFT;
use constant PI => atan2(1,1) * 4;

sub init;
sub calcCORR;
sub getData;
sub padData;
sub getSetData;
sub writeData;

my ($dataFile, $fPrefix); 
my ($fDATA, $pts, $nset, $CORR);

$|++;
($dataFile, $fPrefix) = &init;
print "Parsing data file $dataFile...";
$fDATA = getData($dataFile);
($pts, $nset) = &padData($fDATA);
print "Done\nCalculating correlations for $nset sets...";
$CORR = calcCORR($fDATA, $fPrefix, $pts, $nset);
print "Done\nWriting data to ${fPrefix}_correlation.dat...";
&writeData($CORR, $fPrefix, $nset, $pts, $fDATA);
print "Done\n";

sub writeData {
    my ($corr, $prefix, $num, $pts, $data) = @_;
    my ($i, $j, $k, $saveName, $corrLen, $time);

    $saveName = "${prefix}_correlation.dat";
    $corrLen = int($pts/2) -1;
    open OUTDATA, "> $saveName" or die "Cannot write to $saveName: $!\n";
    for $i (1 .. $num) {
	for $j ($i .. $num) {
	    print OUTDATA "#Set $i autocorrelation\n" if ($i == $j);
	    print OUTDATA "#Sets $i - $j cross correlation\n" if ($i != $j);
	    for $k (0 .. (2*$corrLen - 1)) {
		$time = $data->[$k]{time};
		printf OUTDATA "%8.3f %12.3f\n", $time, $corr->{$i}{$j}->[$k*2]/($corr->{$i}{$j}->[0]);
	    }
	    print OUTDATA "\n";
	}
    }
    close OUTDATA;

}

sub calcCORR {
    my ($data, $prefix, $dataLen, $numSets) = @_;
    my ($set1, $set2, $corrLen, $i, $j, $fft, $CORR); 
    my ($pwr1, $pwr2, $k, $dtmp, $invData, $setcorr, $tot);

    $tot = scalar(@{ $data }) *2;
    for $i (1 .. $numSets) {
	$set1 = getSetData($data, ($i - 1));
	$fft = new Math::FFT($set1);
	$pwr1 = $fft->cdft();
	for $j ($i .. $numSets) {
	    $invData = ();
	    $set2 = getSetData($data, ($j - 1));
	    undef($fft);
	    $fft = new Math::FFT($set2);
	    $pwr2 = $fft->cdft();
	    $k = 0;
	    while ($k < $tot) {
		$dtmp            = $pwr1->[$k  ] * $pwr2->[$k  ] + $pwr1->[$k+1] * $pwr2->[$k+1];
		$invData->[$k+1] = $pwr1->[$k  ] * $pwr2->[$k+1] - $pwr2->[$k  ] * $pwr1->[$k+1];
		$invData->[$k  ] = $dtmp;
		$k += 2;
	    }
	    $fft = new Math::FFT($invData);
	    $setcorr = $fft->invcdft($invData); #autocorrelation function
	    $CORR->{$i}{$j} = $setcorr;
	    undef($setcorr);
	}
	
    }
    return $CORR;
}

sub getSetData {
    my ($data, $setnum) = @_;
    my ($i, $setData);

    for $i (@{ $data }) {
	push @{ $setData }, $i->{VALS}[$setnum];
	push @{ $setData }, 0;
    }

    return $setData;
}

sub padData {
    my ($data) = $_[0];
    my ($dataLen, $i, $numPts, $rec, $j);

    $dataLen = scalar(@{ $data });
    $numPts = scalar(@{ $data->[0]{VALS} });
    #find the smallest power of 2 to include the data, pad with zeros
    $i = 2;
    while ($i < $dataLen) {
	$i *=2;
    }
    for (($dataLen +1) .. $i) {
	$rec = ();
	for $j (1 .. $numPts) {
	    push @{ $rec->{VALS} }, 0;
	}
	push @{ $data }, $rec;
    }
    return ($dataLen, $numPts);
}

sub getData {
    my ($dFile) = $_[0];
    my (@FDATA, $rec, $time, $data);

    open DATAFILE, $dFile or die "ERROR: Cannot open $dFile: $!\n";
    while (<DATAFILE>) {
	chomp;
	$rec = ();
	if ($_ =~ /^\s*(\d+\.?\d*)\s+(\-?\d+\.\d+.*)/) {
	    ($time, $data) = ($1, $2);
	    while ($data =~ /(\-?\d+\.\d+)/g) {
		push @{ $rec->{VALS} }, $1;
	    }
	    $rec->{time} = $time;
	    push @FDATA, $rec;
	}
    }
    close DATAFILE;

    die "ERROR: $dFile does not contain any valid data!\n" if (! @FDATA);

    return \@FDATA;
}

sub init {
    my (%OPTS, $dFile, $prefix);

    getopt('ds',\%OPTS);

    die "usage: $0 -d data file -s [save prefix]\n" if (! defined($OPTS{d}));

    print "Initializing...";
    ($dFile, $prefix) = ($OPTS{d}, $OPTS{s});
    FileTester($dFile);

    if (! defined($prefix)) {
	$prefix = basename($dFile);
	$prefix =~ s/\.\w+$//;
    }

    print "Done\n";

    return ($dFile, $prefix);
}
