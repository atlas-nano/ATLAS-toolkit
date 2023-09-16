#!/usr/bin/perl -w

use strict;
use Getopt::Std qw(getopt);

my (%OPTS, $RES, %short, $i, $j);

getopt('f',\%OPTS);
die "usage: $0 -f fasta string\n" if (! exists($OPTS{f}));
if (-e $OPTS{f}) {
    open FASTAFILE, $OPTS{f} or die "ERROR: Cannot read $OPTS{f}: $!\n";
    while (<FASTAFILE>) {
	chomp;
	while ($_ =~ /(\w)/g) {
	    push @{ $RES }, $1;
	}
    }
    close FASTAFILE;
    die "ERROR: $OPTS{f} is not a valid fastafile!\n" if (! $RES);
} else {
    while ($OPTS{f} =~ /(\w)/g) {
	push @{ $RES }, $1;
    }
}

%short = (
	  "ALA" => "A",
	  "CYS" => "C",
	  "ASP" => "D",
	  "GLU" => "E",
	  "PHE" => "F",
	  "GLY" => "G",
	  "HIS" => "H",
	  "ILE" => "I",
	  "LYS" => "K",
	  "LEU" => "L",
	  "MET" => "M",
	  "ASN" => "N",
	  "PRO" => "P",
	  "GLN" => "Q",
	  "ARG" => "R",
	  "SER" => "S",
	  "THR" => "T",
	  "VAL" => "V",
	  "TRP" => "W",
	  "TYR" => "Y" );
for $i (@{ $RES }) {
    for $j (keys %short) {
	if (lc($i) eq lc($short{$j})) {
	    print "$j\n";
	}
    }
}
