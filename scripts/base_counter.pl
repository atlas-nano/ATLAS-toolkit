#!/usr/bin/perl -w
use strict;
use warnings;

if (! @ARGV) {
    die "$0 sequencefile\n";
}

my ($sequencefile) = $ARGV[0];
my ($cg_counter, $at_counter);

-e $sequencefile or die "Cannot locate $sequencefile: $!\n";

open SEQFILE, $sequencefile or die "Cannot open $sequencefile: $!\n";
while (<SEQFILE>) {
    if ($_ =~ /[CG|GC]/i) {
	$cg_counter++;
    } elsif ($_ =~ /[AT|TA]/i) {
	$at_counter++;
    }
}

close SEQFILE;

if ($at_counter or $cg_counter) {
    print "$at_counter AT units\n$cg_counter CG units\n";
} else {
    print "$sequencefile is not a valid sequence file\n";
}
