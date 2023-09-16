#!/usr/bin/perl -w
use strict;
use warnings;

my ($counter, $total);

if (@ARGV) {
    open INFILE, $ARGV[0] or die "Cannot open $ARGV[0]: $!\n";
    $counter = 0;
    $total = 0.0;
    while (<INFILE>) {
	if ($_ =~ /\d+\s*(-?\d+\.\d+)/) {
	    chomp;
	    $total += $1;
	    $counter++;
	}
    }
    close INFILE;

    if ($counter > 0) {
	printf "Total Energy: %10.2f\nTotal Bases: %5d\nAverage Energy per base: %10.2f\n", $total, $counter, ($total/$counter);
    } else {
	print "Invalid file: $ARGV[0]\n";
    }
} else {
    die "usage: avg_energy.pl datfile\n";
}
