#!/usr/bin/perl -w

use strict;

if (!@ARGV) {
    die "usage: stripNaH2O.pl pdbfile [savefile]\n";
}

my (@outfile, $i, $instr);

my ($in_name, $out_name) = @ARGV;

-e $in_name or die "Cannot find $in_name: $!\n";

$out_name = $in_name if (! $out_name);

open INFILE, $in_name or die "Cannot open $in_name: $!\n";
while (<INFILE>) {
    chomp;
    if ($_ =~ /Na|WAT|Mg/i) {
	last;
    } else {
	push @outfile, $_;
    }
}

close INFILE;

open OUTFILE, "> $out_name" or die "Cannot write changes to $out_name: $!\n";
for $i (@outfile) {
    print OUTFILE "$i\n";
}

close OUTFILE;
