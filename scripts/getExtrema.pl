#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);

my (%OPTS, $min, $i, $max, $datafile, %DATA);

getopt('d',\%OPTS);

die "usage: $0 -d data file\n" if (! exists($OPTS{d}));

$datafile = $OPTS{d};

open DATAFILE, $datafile || die "ERROR: Cannot open $datafile: $!\n";
while (<DATAFILE>) {
    chomp;
    if ($_ =~ /^\s*(\d+\.\d+)\s+(\-?\d+\.\d+)/) {
	$DATA{$1} = $2;
    }
}
close DATAFILE;

die "ERROR: Invalid data file $datafile\n" if (! %DATA);

for $i (keys %DATA) {
    if (! defined($max) or $DATA{$i} > $max->{Y}) {
	$max->{Y} = $DATA{$i};
	$max->{X} = $i;
     }
     if (! defined($min) or $DATA{$i} < $min->{Y}) {
	$min->{Y} = $DATA{$i};
	$min->{X} = $i;
      }
}

print "MIN $min->{X} $min->{Y}\nMAX $max->{X} $max->{Y}\n";
