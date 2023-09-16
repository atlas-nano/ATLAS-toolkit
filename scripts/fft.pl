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
sub calcFFT;
sub getData;

my ($dataFile, $fPrefix, $CORR);

$|++;
($dataFile, $fPrefix) = &init;
print "Parsing data file $dataFile...";
$CORR = getData($dataFile);
print "Done\n";
&calcFFT($CORR, $fPrefix);

sub getData {
    my ($dFile) = $_[0];
    my (@FDATA);

    open DATAFILE, $dFile or die "ERROR: Cannot open $dFile: $!\n";
    while (<DATAFILE>) {
	chomp;
	if ($_ =~ /^\d+\.?\d* (\-?\d+\.\d+)/) {
	    push @FDATA, $1;
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
sub calcFFT {
    my ($correlation, $prefix) = @_;
    my ($fft, $i, $rFFT, $cFFT, $spectrum);

    print "Writing data...dipole fourier transform        \r";
    $fft = new Math::FFT($correlation);
    $spectrum = $fft->rdft();
    open rFFT, "> ${prefix}_rfft.dat" or die "ERROR: Cannot create ${prefix}_rfft.dat: $!\n"; # real part of fft
    open cFFT, "> ${prefix}_cfft.dat" or die "ERROR: Cannot create ${prefix}_cfft.dat: $!\n"; # complex part of fft
    for $i (0 .. ($#{ $spectrum } - 1)/2) {
        $rFFT->[$i] = $spectrum->[2*$i];
        print rFFT "$i $rFFT->[$i]\n";
        next if ($i == 0);
        $cFFT->[$i] = $spectrum->[2*$i + 1];
        print cFFT "$i $cFFT->[$i]\n";
    }
    close rFFT;
    close cFFT;
    $spectrum = ();
    $correlation = ();

    print "Writing data...frequency dependent dipole\r";
    open DIPOLEFREQ, "> ${prefix}_dipole_freq.dat" or die "ERROR: Cannot create ${prefix}_dipole_freq.dat: $!\n";
    for $i (1 .. $#{ $cFFT }) {
        print DIPOLEFREQ "$i " . (1 - $i * $cFFT->[$i]) . "\n";
    }
    close DIPOLEFREQ;
    $cFFT = ();
    $rFFT = ();

    print "Writing data...Done                        \n";
}
